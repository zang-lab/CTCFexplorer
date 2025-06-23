#!/bin/bash

# ── Manual configuration ────────────────────────────────
DATA_DIR="/home/vdj6tq/pg17.2_data"
PGHOST="$DATA_DIR"
PGPORT=55432
PGDATABASE="CTCFexplorer"
PGUSER="vdj6tq"
LOGFILE="$DATA_DIR/postgres.log"
# ─────────────────────────────────────────────────────────

export PGHOST PGPORT PGDATABASE PGUSER

echo "Using data directory: $DATA_DIR"
echo "Connecting as user: $PGUSER to database: $PGDATABASE on port: $PGPORT"

# ── Ensure the master server is running ─────────────────
if ! pg_ctl -D "$DATA_DIR" status >/dev/null 2>&1; then
    echo "Starting PostgreSQL server..."
    pg_ctl -D "$DATA_DIR" \
           -o "-p $PGPORT -k $PGHOST" \
           -l "$LOGFILE" start
    sleep 2
fi
# ───────────────────────────────────────────────────────────────────────

DUMP_DIR=app/pg_dumps
SCHEMA_SQL=app/schema.sql

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

# create index
echo "Running index creation …"
python create_union_id_indexes.py

echo "All tasks finished."

# create supplementary
echo "Running index creation …"
python create_supplementary_tables.py

echo "All tasks finished."