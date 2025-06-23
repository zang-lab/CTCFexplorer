#!/bin/bash
set -euo pipefail

# ── Load required environment variables ─────────────────
: "${PGHOST:?PGHOST is not set}"
: "${PGPORT:?PGPORT is not set}"
: "${PGDATABASE:?PGDATABASE is not set}"
: "${PGUSER:?PGUSER is not set}"
: "${LOGFILE:?LOGFILE is not set}"
# ────────────────────────────────────────────────────────

export PGHOST PGPORT PGDATABASE PGUSER LOGFILE

echo "Using PGHOST: $PGHOST"
echo "Connecting as user: $PGUSER to database: $PGDATABASE on port: $PGPORT"

# ── Ensure the PostgreSQL server is running ─────────────
if ! pg_ctl -D "$PGHOST" status >/dev/null 2>&1; then
    echo "Starting PostgreSQL server..."
    pg_ctl -D "$PGHOST" \
           -o "-p $PGPORT -k $PGHOST" \
           -l "$LOGFILE" start
    sleep 2
fi
# ────────────────────────────────────────────────────────

DUMP_DIR=app/pg_dumps
SCHEMA_SQL=app/schema.sql

# Check schema file exists
if [[ ! -f "$SCHEMA_SQL" ]]; then
    echo "ERROR: schema file not found at $SCHEMA_SQL"
    exit 1
fi

# Check for dump files
NUM_DUMPS=$(ls -1 "$DUMP_DIR"/dump_*.sql 2>/dev/null | wc -l)
if [[ $NUM_DUMPS -eq 0 ]]; then
    echo "ERROR: no dump_*.sql files found in $DUMP_DIR"
    exit 1
fi

# Load schema first
echo "Loading schema from $SCHEMA_SQL ..."
psql -v ON_ERROR_STOP=1 -f "$SCHEMA_SQL" || {
    echo "ERROR: Failed to load schema."
    exit 1
}

# Load dumps in sorted order
echo "Replaying $NUM_DUMPS dump files ..."
for f in $(ls "$DUMP_DIR"/dump_*.sql | sort -V); do
    echo ">> $(basename "$f")"
    psql -v ON_ERROR_STOP=1 -f "$f" || {
        echo "ERROR in $f"
        exit 1
    }
done

echo "Merge job complete."

# Create index
echo "Running index creation …"
python create_union_id_indexes.py

echo "All tasks finished."

# Create supplementary tables
echo "Running index creation …"
python create_supplementary_tables.py

echo "All tasks finished."