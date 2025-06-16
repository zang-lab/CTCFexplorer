#!/usr/bin/env python3
import os, psycopg2

# ── pick up TCP details from the environment ───────────────────────────
PGHOST     = os.environ["PGHOST"]
PGPORT     = os.environ["PGPORT"]
PGDATABASE = os.environ["PGDATABASE"]
PGUSER     = os.environ["PGUSER"]
# ────────────────────────────────────────────────────────────────────────

# Connect via TCP
conn = psycopg2.connect(
    host     = PGHOST,
    port     = PGPORT,
    dbname   = PGDATABASE,
    user     = PGUSER,
)
# ────────────────────────────────────────────────────────────────────────

# autocommit for CREATE INDEX CONCURRENTLY
conn.autocommit = True

with conn.cursor() as cur:
    print("Creating idx_union_id_basic …")
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