#!/usr/bin/env python3
import os
from pathlib import Path

import pandas as pd
import psycopg2
from psycopg2.extras import execute_values

# ── CONFIGURATION ────────────────────────────────────────────────────────
SUPP_DIR = Path(__file__).parent.resolve()

DB_HOST     = os.environ["PGHOST"]
DB_PORT     = os.environ["PGPORT"]
DB_NAME     = os.environ["PGDATABASE"]
DB_USER     = os.environ["PGUSER"]
# ────────────────────────────────────────────────────────────────────────

def get_connection():
    return psycopg2.connect(
        host=DB_HOST,
        port=DB_PORT,
        dbname=DB_NAME,
        user=DB_USER
    )

def create_ctcf_labels_table(conn):
    csv_path = SUPP_DIR / "sample_summary_20250528.csv"
    if not csv_path.exists():
        print(f"[Skip] CTCFLabels: {csv_path} not found.")
        return

    df = pd.read_csv(csv_path)
    cols = list(df.columns)

    # Build and run DDL
    cols_def = ", ".join(f'"{c}" TEXT' for c in cols)
    ddl = f'CREATE TABLE IF NOT EXISTS "CTCFLabels" ({cols_def});'
    with conn.cursor() as cur:
        cur.execute(ddl)
        cur.execute('DELETE FROM "CTCFLabels";')

        # Prepare INSERT
        col_list = ", ".join(f'"{c}"' for c in cols)
        insert_sql = f'INSERT INTO "CTCFLabels" ({col_list}) VALUES %s'
        execute_values(cur, insert_sql, list(df.itertuples(index=False, name=None)))
        print(f"[Done] CTCFLabels: loaded {len(df):,} rows.")

def create_hg38gene_table(conn):
    csv_path = SUPP_DIR / "hg38_gene_annotation_GRCh38.p14_v2.2.csv"
    if not csv_path.exists():
        print(f"[Skip] hg38gene: {csv_path} not found.")
        return

    df = pd.read_csv(
        csv_path,
        header=None,
        names=["Chromosome", "Start", "End", "Gene", "Strand"]
    )

    ddl = """
    CREATE TABLE IF NOT EXISTS "hg38gene" (
      "Chromosome" TEXT,
      "Start"      INTEGER,
      "End"        INTEGER,
      "Gene"       TEXT,
      "Strand"     TEXT
    );
    """
    with conn.cursor() as cur:
        cur.execute(ddl)
        cur.execute('DELETE FROM "hg38gene";')

        insert_sql = """
        INSERT INTO "hg38gene"
          ("Chromosome","Start","End","Gene","Strand")
        VALUES %s
        """
        execute_values(cur, insert_sql, list(df.itertuples(index=False, name=None)))
        print(f"[Done] hg38gene: loaded {len(df):,} rows.")

def print_table_counts(conn):
    with conn.cursor() as cur:
        cur.execute("""
            SELECT table_name
              FROM information_schema.tables
             WHERE table_schema = 'public'
               AND table_type   = 'BASE TABLE'
             ORDER BY table_name;
        """)
        tables = [r[0] for r in cur.fetchall()]

        print("\n Table Record Counts:")
        for tbl in tables:
            cur.execute(f'SELECT COUNT(*) FROM "{tbl}";')
            cnt = cur.fetchone()[0]
            print(f" - {tbl}: {cnt:,} rows")

def main():
    print(f"Connecting to {DB_NAME}@{DB_HOST}:{DB_PORT} as {DB_USER}")
    conn = get_connection()
    conn.autocommit = True

    create_ctcf_labels_table(conn)
    create_hg38gene_table(conn)
    print_table_counts(conn)

    conn.close()
    print("\n All supplementary tables loaded.")

if __name__ == "__main__":
    main()