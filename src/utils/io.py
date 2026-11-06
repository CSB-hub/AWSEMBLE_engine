"""
파일 입출력 유틸리티 (FASTA, GenBank)

이 모듈은 DNA 서열 파일을 읽어서 pydna Dseqrecord로 변환하고 관리합니다.

클래스 구성:
    - DNAFileReader: FASTA/GenBank 파일 읽기 (배치 처리 포함)
    - DNAFileWriter: FASTA/GenBank 파일 쓰기
    - DNAValidator: 조립된 서열의 정확성 및 순서 검증
    - DatabaseImport: DB 서열 가져오기 및 스키마 변환 유틸리티

작성일: 2025-10-30
수정일: 2026-03-03 (클래스 기반 전면 리팩토링)
"""

from pathlib import Path
from typing import Union, List, Dict, Optional, Tuple
import os
import random
import mysql.connector
from mysql.connector import Error

# BioPython imports
from Bio import SeqIO

# pydna imports
from pydna.dseqrecord import Dseqrecord

# awsemble 전용 annotation(part_type 등)은 GenBank 표준이 아니라서 Biopython이 파일에
# 기록하지도, 읽어서 복원하지도 않는다. Galaxy job들은 서로 **파일로만** 데이터를 주고받으므로
# (B Import → C Arrange/Build → D Insert가 각각 별도 프로세스), 그냥 두면 매번 유실된다.
# (part_type 유실 → module feature가 전부 misc_feature로 떨어지던 버그가 이 경우)
# → 왕복 보존되는 COMMENT 필드에 `awsemble_<키>:<값>` 마커 한 줄씩 실어 보내고, 읽을 때 복원한다.
# 마커 형식은 genbank_part_search/bulk_part_register가 이미 쓰는 관례와 동일하게 맞춘다.
AWSEMBLE_MARKER_PREFIX = "awsemble_"

# 왕복 보존할 awsemble annotation 키:
#   part_type    : 부품 종류(promoter/rbs/cds/…) → module part feature의 GenBank type 결정
#   insert_whole : "이 insert는 통째로 vector에 넣는다"(biopart module/multimodule) 표시.
#                  Gibson/GG의 CDS 자동감지(_find_target_feature)가 module 안의 CDS 부품을
#                  삽입 대상으로 오인해 그 구간만 증폭하는 것을 막는다.
ROUNDTRIP_ANNOTATION_KEYS = ("part_type", "insert_whole")

# 부품 gbk가 feature note에 같은 마커를 다는 기존 관례(genbank_part_search_helper 참고)
PART_TYPE_MARKER = AWSEMBLE_MARKER_PREFIX + "part_type:"


def _comment_to_text(comment) -> str:
    """Biopython의 annotations["comment"]는 str 또는 list[str] 둘 다 올 수 있다."""
    if not comment:
        return ""
    if isinstance(comment, (list, tuple)):
        return "\n".join(str(c) for c in comment)
    return str(comment)

# TODO: 아래 코드들은 DB기준으로 변경이 될때 내가 작업하지 않았기 때문에 1차적인 코드설명이 필요함.
# ============================================================================ 
# 데이터베이스 관리 및 입출력 클래스
# ============================================================================ 

class DatabaseImport:
    """
    데이터베이스 스키마와 상호작용하고 서열 데이터를 관리하는 클래스.
    """

    def __init__(self):
        """환경 변수를 확인하고 데이터베이스 접속 정보를 설정합니다."""
        # .env 파일 로드 (awsemble/.env 위치 우선 검색, 시스템 변수 덮어쓰기 허용)
        from dotenv import load_dotenv
        env_path = Path(__file__).resolve().parent.parent.parent / ".env"
        load_dotenv(dotenv_path=env_path, override=True)

        self.host = os.getenv("HOST")
        self.port = os.getenv("PORT")
        self.user = os.getenv("USER")
        self.password = os.getenv("PASSWORD")

        if not all([self.host, self.port, self.user, self.password]):
            raise ConnectionError(
                "HOST, PORT, USER and PASSWORD must all be set. "
                f"(.env path searched: {env_path})"
            )
        
        self.conn = None

    def _connect(self, schema_name: str) -> bool:
        """지정된 스키마의 데이터베이스에 연결합니다."""
        try:
            self.conn = mysql.connector.connect(
                host=self.host,
                port=int(self.port),
                user=self.user,
                password=self.password,
                database=schema_name
            )
            return self.conn.is_connected()
        except Error as e:
            print(f"Error connecting to database '{schema_name}': {e}")
            self.conn = None
            return False

    def _disconnect(self):
        """데이터베이스 연결을 닫습니다."""
        if self.conn and self.conn.is_connected():
            self.conn.close()

    def convert_schema(self, schema_name: str, output_dir: Union[str, Path]):
        """스키마의 모든 테이블을 각각의 FASTA 파일로 변환합니다."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if not self._connect(schema_name):
            return

        try:
            cursor = self.conn.cursor()
            cursor.execute("SHOW TABLES;")
            tables = cursor.fetchall()
            print(f"Found {len(tables)} tables in '{schema_name}'. Starting conversion...")

            for (table_name,) in tables:
                cursor.execute(f"SHOW COLUMNS FROM `{table_name}`")
                columns = [col[0].lower() for col in cursor.fetchall()]

                name_col = next((c for c in ['name', 'id', 'part_name', 'identifier', 'test_id'] if c in columns), None)
                seq_col = next((c for c in ['sequence', 'seq', 'dna_sequence', 'seq_data'] if c in columns), None)

                if not name_col or not seq_col:
                    continue
                
                cursor.execute(f"SELECT `{name_col}`, `{seq_col}` FROM `{table_name}`")
                rows = cursor.fetchall()
                dseqs = [Dseqrecord(seq, name=str(name)) for name, seq in rows if name and seq]

                if dseqs:
                    output_file = output_dir / f"{table_name}.fasta"
                    DNAFileWriter.write_fasta(dseqs, output_file)
            cursor.close()
        except Error as e:
            print(f"An error occurred: {e}")
        finally:
            self._disconnect()

    def _execute(self, schema_name: str, query: str, params: tuple = (), fetch: str = "all"):
        """공통 DB 실행 헬퍼 (L-1): connect→cursor→execute→fetch→disconnect 중복 일원화.

        Args:
            schema_name: 접속 스키마
            query: SQL (placeholder %s)
            params: 바인딩 파라미터 튜플
            fetch: 'all'→fetchall() / 'one'→fetchone() / None→결과 없음

        Returns:
            fetch='all' → list(실패 시 []), 'one' → row 또는 None, None → None
        """
        empty = [] if fetch == "all" else None
        if not self._connect(schema_name):
            return empty
        try:
            cursor = self.conn.cursor()
            cursor.execute(query, params)
            if fetch == "all":
                result = cursor.fetchall()
            elif fetch == "one":
                result = cursor.fetchone()
            else:
                result = None
            cursor.close()
            return result
        except Error as e:
            print(f"An error occurred: {e}")
            return empty
        finally:
            self._disconnect()

    def get_part_by_id(self, schema_name: str, table_name: str, part_id: str) -> Optional[Dseqrecord]:
        """데이터베이스에서 ID로 특정 부품을 가져옵니다."""
        result = self._execute(
            schema_name, f"SELECT sequence FROM `{table_name}` WHERE `ID` = %s", (part_id,), fetch="one")
        if not result:
            return None
        dseq = Dseqrecord(result[0], name=part_id)
        # table_name이 곧 part_type(promoter/rbs/cds/terminator …) — module feature 기입(P6)에서
        # part 본체 feature를 GenBank 표준 type으로 달기 위해 부품에 type을 실어 전파한다.
        dseq.annotations["part_type"] = table_name
        return dseq

    def get_all_test_ids_from_construct(self) -> List[str]:
        """Institution이 'igem' 또는 'KRIBB'인 Test_ID를 가져옵니다."""
        rows = self._execute(
            "Test", "SELECT Test_ID FROM test_construct WHERE Institution IN ('igem', 'KRIBB')")
        return [r[0] for r in rows]

    def get_ordered_parts_by_test_id(self, test_id: str) -> List[Tuple[str, str]]:
        """Test_ID에 해당하는 순서대로 정렬된 부품 목록을 가져옵니다."""
        return self._execute(
            "Test", "SELECT part_id, part_type FROM test_part WHERE Test_ID = %s ORDER BY part_order ASC", (test_id,))

    def get_all_overhang_ids(self) -> List[str]:
        """'Bioparts' 스키마의 'overhang' 테이블에서 Institution이 'Awsemble'인 모든 오버행 ID를 가져옵니다."""
        rows = self._execute(
            "Bioparts", "SELECT ID FROM `overhang` WHERE `institution` = 'Awsemble' ORDER BY ID ASC")
        return [r[0] for r in rows]

    def get_test_ids_by_institution(self, institution: str) -> List[str]:
        """`Test` 스키마에서 특정 Institution에 해당하는 Test_ID를 가져옵니다."""
        rows = self._execute(
            "Test", "SELECT Test_ID FROM test_construct WHERE Institution = %s", (institution,))
        return [r[0] for r in rows]

    def get_construct_sequence_by_test_id(self, test_id: str) -> Optional[str]:
        """`Test` 스키마에서 정답 서열을 가져옵니다."""
        result = self._execute(
            "Test", "SELECT sequence FROM test_construct WHERE Test_ID = %s", (test_id,), fetch="one")
        return result[0] if result else None

    def bioparts_to_fasta(self, output_dir: Union[str, Path]):
        """'Bioparts' 스키마의 모든 테이블을 FASTA 파일로 저장합니다."""
        self.convert_schema("Bioparts", output_dir)

    def schema_test_to_fasta(self, output_dir: Union[str, Path]):
        """'Test' 스키마의 모든 테이블을 FASTA 파일로 저장합니다."""
        self.convert_schema("Test", output_dir)


# ============================================================================ 
# 파일 읽기 클래스
# ============================================================================ 

# TODO: Dead code인지 체크하고 어느 단계에서 사용되는지 확인.
class DNAFileReader:
    """
    FASTA 및 GenBank 파일을 읽는 유틸리티 클래스.
    """

    @staticmethod
    def read_fasta(file_path: Union[str, Path]) -> List[Dseqrecord]:
        """FASTA 파일 읽기 → Dseqrecord 리스트"""
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        try:
            records = list(SeqIO.parse(file_path, "fasta"))
            if not records:
                raise ValueError(f"The FASTA file has no sequence: {file_path}")

            dseqrecords = []
            for record in records:
                dseq = Dseqrecord(str(record.seq), circular=False)
                dseq.name = record.id
                dseq.description = record.description
                dseqrecords.append(dseq)
            return dseqrecords
        except Exception as e:
            raise ValueError(f"FASTA parse error ({file_path}): {e}")

    @staticmethod
    def read_genbank(file_path: Union[str, Path]) -> List[Dseqrecord]:
        """GenBank 파일 읽기 → Dseqrecord 리스트 (Features 및 Circular 여부 보존)"""
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        try:
            records = list(SeqIO.parse(file_path, "genbank"))
            if not records:
                raise ValueError(f"The GenBank file has no sequence: {file_path}")

            return [
                DNAFileReader._restore_annotations(Dseqrecord.from_SeqRecord(record))
                for record in records
            ]
        except Exception as e:
            raise ValueError(f"GenBank parse error ({file_path}): {e}")

    @staticmethod
    def _restore_annotations(dseq: Dseqrecord) -> Dseqrecord:
        """gbk 왕복에서 유실된 awsemble annotation을 복원한다 (제자리 수정 후 반환).

        - 1순위: write_genbank가 COMMENT에 남긴 `awsemble_<키>:<값>` 마커
        - 2순위(part_type 한정): 같은 마커를 feature note에 다는 부품 gbk
          (genbank_part_search_helper가 만드는 형식 — bulk_part_register_helper 참고)

        마커가 없으면(외부에서 받은 일반 gbk 등) 아무것도 하지 않는다.
        """
        for line in _comment_to_text(dseq.annotations.get("comment")).splitlines():
            line = line.strip()
            if not line.startswith(AWSEMBLE_MARKER_PREFIX) or ":" not in line:
                continue
            key, _, value = line[len(AWSEMBLE_MARKER_PREFIX):].partition(":")
            if key in ROUNDTRIP_ANNOTATION_KEYS and not dseq.annotations.get(key):
                dseq.annotations[key] = value.strip()

        if not dseq.annotations.get("part_type"):
            for feat in dseq.features:
                for note in feat.qualifiers.get("note", []):
                    note = str(note).strip()
                    if note.startswith(PART_TYPE_MARKER):
                        dseq.annotations["part_type"] = note[len(PART_TYPE_MARKER):].strip()
                        return dseq

        return dseq

    @classmethod
    def read_genbank_dir(cls, directory: Union[str, Path]) -> Dict[str, Dseqrecord]:
        """디렉토리 내 모든 GenBank 파일 배치 읽기 (.gb 확장자)"""
        directory = Path(directory)
        if not directory.is_dir():
            raise ValueError(f"Not a directory: {directory}")

        result = {}
        for gb_file in sorted(directory.glob("*.gb")):
            try:
                dseqs = cls.read_genbank(gb_file)
                if dseqs:
                    result[gb_file.name] = dseqs[0]
            except Exception as e:
                print(f"Could not read file ({gb_file.name}): {e}")
        return result


# ============================================================================ 
# 파일 쓰기 클래스
# ============================================================================ 

# TODO: 유저에게 마지막 output 줄때 사용?
class DNAFileWriter:
    """
    DNA 서열을 파일로 저장하는 유틸리티 클래스.
    """

    @staticmethod
    def write_genbank(dseq: Dseqrecord, file_path: Union[str, Path]) -> Path:
        """Dseqrecord를 GenBank 파일로 저장 (Features 포함)

        pydna `dseq.write()`는 대상 파일이 이미 존재하면 먼저 그 파일을 읽어서
        서열을 비교하는데(다르면 _OLD_ 백업), Galaxy처럼 출력 파일을 빈 파일로
        미리 만들어두는 환경에선 그 빈 파일을 파싱하다
        "ValueError: not enough values to unpack"로 깨진다. 조립 결과는 항상
        새로 쓰면 되므로 비교/백업 로직 없이 직접 포맷해서 쓴다.
        (프로덕션은 새 경로에 써서 원래 문제가 없지만, galaxy와 동일 코드로 통일해
         src→galaxy 동기화 시 재발을 막는다.)
        """
        file_path = Path(file_path)
        file_path.parent.mkdir(parents=True, exist_ok=True)

        # awsemble annotation은 표준 필드가 아니라 그냥 쓰면 유실된다 → COMMENT 마커로 실어
        # 보낸다(ROUNDTRIP_ANNOTATION_KEYS 주석 참고). 호출자의 record는 원복해 안 건드린다.
        markers = [
            f"{AWSEMBLE_MARKER_PREFIX}{k}:{dseq.annotations[k]}"
            for k in ROUNDTRIP_ANNOTATION_KEYS
            if dseq.annotations.get(k)
        ]
        prev_comment = dseq.annotations.get("comment")
        if markers:
            kept = [
                ln for ln in _comment_to_text(prev_comment).splitlines()
                if not ln.strip().startswith(AWSEMBLE_MARKER_PREFIX)  # 기존 마커는 갱신(중복 방지)
            ]
            dseq.annotations["comment"] = "\n".join(kept + markers)

        try:
            with open(file_path, "w", encoding="utf8") as f:
                f.write(dseq.format("gb"))
        finally:
            if markers:
                if prev_comment is None:
                    dseq.annotations.pop("comment", None)
                else:
                    dseq.annotations["comment"] = prev_comment
        return file_path

    @staticmethod
    def write_fasta(dseqs: List[Dseqrecord], file_path: Union[str, Path]) -> Path:
        """Dseqrecord 리스트를 단일 FASTA 파일로 저장"""
        file_path = Path(file_path)
        file_path.parent.mkdir(parents=True, exist_ok=True)

        with open(file_path, 'w') as f:
            for dseq in dseqs:
                f.write(f">{dseq.name if dseq.name else 'sequence'}\n")
                seq_str = str(dseq.seq)
                for i in range(0, len(seq_str), 60):
                    f.write(seq_str[i:i+60] + "\n")
        return file_path


# ============================================================================ 
# 검증 클래스
# ============================================================================ 

class DNAValidator:
    """
    조립된 서열의 구성 및 정확성을 검증하는 클래스.
    """

    def __init__(self, db_importer: Optional[DatabaseImport] = None):
        self.db = db_importer or DatabaseImport()

    def check_parts_correct(self, gbk_path: Union[str, Path], test_id: str) -> bool:
        """결과 GenBank 파일이 DB의 부품 순서 및 서열과 일치하는지 정밀 검증"""
        try:
            results = DNAFileReader.read_genbank(gbk_path)
        except Exception:
            print(f"[{test_id}] Could not read the result file.")
            return False
        
        full_seq = str(results[0].seq).upper()
        ordered_parts = self.db.get_ordered_parts_by_test_id(test_id)
        
        if not ordered_parts:
            print(f"[{test_id}] Could not find the part composition in the DB.")
            return False

        current_pos = 0
        success = True
        print(f"\n[{test_id}] Starting detailed validation...")
        
        for i, (part_id, part_type) in enumerate(ordered_parts, 1):
            part_record = self.db.get_part_by_id("Bioparts", part_type, part_id)
            if not part_record:
                print(f"  [{i}] {part_id}: sequence missing in the DB")
                success = False
                continue
                
            ref_seq = str(part_record.seq).upper()
            # 이전 부품 이후의 위치에서 서열 검색
            match_pos = full_seq.find(ref_seq, max(0, current_pos - 100))
            
            if match_pos == -1:
                print(f"  [{i}] {part_id}: sequence or order mismatch")
                success = False
                continue
                
            print(f"  [{i}] {part_id}: Pass (at {match_pos})")
            current_pos = match_pos + len(ref_seq)

        if success:
            print(f"[{test_id}] All checks passed (order and sequence match)")
        else:
            print(f"[{test_id}] Some checks failed")
            
        return success


# ============================================================================ 
# 모듈 메타데이터
# ============================================================================ 

__all__ = [
    'DNAFileReader',
    'DNAFileWriter',
    'DNAValidator',
    'DatabaseImport',
]
