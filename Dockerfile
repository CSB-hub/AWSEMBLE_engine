# awsemble 테스트/실행 환경
# - spark(docker 기반) 머신에서 회귀테스트를 돌리기 위한 이미지
# - 비밀번호(.env)는 이미지에 굽지 않고 런타임에 마운트한다.
#
# 빌드:
#   docker build -t awsemble:dev awsemble/
# 회귀테스트(전수):
#   docker run --rm --network host -v "$PWD/awsemble/.env:/app/.env:ro" awsemble:dev
# 단일 Test_ID:
#   docker run --rm --network host -v "$PWD/awsemble/.env:/app/.env:ro" awsemble:dev \
#       python tests/igem/test_single.py <Test_ID>
#
# DB 접속: --network host 로 호스트 네트워크 공유 → .env 의 HOST(10.131.81.3:3306) 에 직접 연결.

FROM python:3.11-slim

# pydna/biopython 등이 빌드 시 필요로 할 수 있는 최소 런타임 라이브러리
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# 의존성 먼저 설치 (레이어 캐시)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 소스/테스트 복사 (.env 는 .dockerignore 로 제외 → 런타임 마운트)
COPY __init__.py .
COPY src/ ./src/
COPY tests/ ./tests/

# 테스트가 `from utils.io ...`, `from workflow.dbtl ...` 로 import → src 를 경로에 추가
ENV PYTHONPATH=/app/src
ENV PYTHONUNBUFFERED=1

# 기본 동작: iGEM 전수 회귀테스트 (docker run 인자로 override 가능)
CMD ["python", "tests/igem/test_all.py"]
