#!/usr/bin/env python3
import os, psycopg2

# ── pick up TCP details from the environment ───────────────────────────
DB_HOST     = os.environ["DB_HOST"]
DB_PORT     = os.environ["DB_PORT"]
DB_NAME     = os.environ["DB_NAME"]
DB_USER     = os.environ["DB_USER"]
DB_PASSWORD = os.environ["DB_PASSWORD"]
# ────────────────────────────────────────────────────────────────────────

# Connect via TCP
conn = psycopg2.connect(
    host     = DB_HOST,
    port     = DB_PORT,
    dbname   = DB_NAME,
    user     = DB_USER,
    password = DB_PASSWORD
)
# ────────────────────────────────────────────────────────────────────────

# autocommit for CREATE INDEX CONCURRENTLY
conn.autocommit = True

with conn.cursor() as cur:
    print("Creating idx_union_id_basic …")

    cur.execute("""
        VACUUM ANALYZE "BasicInfo", "SampleInfo", "CelltypeInfo";
    """)

    cur.execute("""
        CREATE INDEX CONCURRENTLY IF NOT EXISTS idx_union_id_basic
          ON "BasicInfo" ("Union ID");
    """)

    print("Creating idx_union_id_celltype …")
    cur.execute("""
        CREATE INDEX CONCURRENTLY IF NOT EXISTS idx_union_id_celltype
          ON "CelltypeInfo" ("Union ID");
    """)

    print("Creating idx_union_id_sample …")
    cur.execute("""
        CREATE INDEX CONCURRENTLY IF NOT EXISTS idx_union_id_sample
          ON "SampleInfo" ("Union ID");
    """)

conn.close()
print("All indexes created.")