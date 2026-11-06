
import sys, os
sys.path.append(os.path.abspath("Awsemble/awsemble_dev/BioCAD_private/awsemble/src"))
from utils.io import DatabaseImport

db = DatabaseImport()
db._connect("Bioparts")
cursor = db.conn.cursor()
cursor.execute("SHOW TABLES")
tables = [r[0] for r in cursor.fetchall()]

target = "BBa_K115001"
found = False
for t in tables:
    try:
        cursor.execute(f"SELECT ID FROM `{t}` WHERE `ID` = %s", (target,))
        res = cursor.fetchone()
        if res:
            print(f"Found exactly '{target}' in table: {t}")
            found = True
        else:
            # Try partial match just in case
            cursor.execute(f"SELECT ID FROM `{t}` WHERE `ID` LIKE %s", (f"%{target}%",))
            res = cursor.fetchall()
            if res:
                print(f"Found partial match for '{target}' in table: {t} -> {res}")
                found = True
    except Exception:
        continue

if not found:
    print(f"'{target}' not found in any table of Bioparts schema.")
db._disconnect()
