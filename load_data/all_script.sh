#!/bin/bash
set -euo pipefail

# ── Load required environment variables ─────────────────
: "${DB_HOST:?DB_HOST is not set}"
: "${DB_PORT:?DB_PORT is not set}"
: "${DB_NAME:?DB_NAME is not set}"
: "${DB_USER:?DB_USER is not set}"
: "${DB_PASSWORD:?DB_PASSWORD is not set}"
: "${LOGFILE:?LOGFILE is not set}"
# ────────────────────────────────────────────────────────

export DB_HOST DB_PORT DB_NAME DB_USER DB_PASSWORD LOGFILE

echo "Using DB_HOST: $DB_HOST"
echo "Connecting as user: $DB_USER to database: $DB_NAME on port: $DB_PORT"

DUMP_DIR=/app/pg_dumps
SCHEMA_SQL=/app/schema.sql

# Check schema file exists
if [[ ! -f "$SCHEMA_SQL" ]]; then
    echo "ERROR: Schema file not found at $SCHEMA_SQL"
    exit 1
fi

# Check for dump files
NUM_DUMPS=$(ls -1 "$DUMP_DIR"/dump_*.sql 2>/dev/null | wc -l)
if [[ $NUM_DUMPS -eq 0 ]]; then
    echo "ERROR: No dump_*.sql files found in $DUMP_DIR"
    exit 1
fi

# Load schema
echo "Loading schema from $SCHEMA_SQL ..."
PGPASSWORD="$DB_PASSWORD" psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -v ON_ERROR_STOP=1 -f "$SCHEMA_SQL" || {
    echo "ERROR: Failed to load schema."
    exit 1
}

# Load dumps
echo "Replaying $NUM_DUMPS dump files ..."
for f in $(ls "$DUMP_DIR"/dump_*.sql | sort -V); do
    echo ">> $(basename "$f")"
    PGPASSWORD="$DB_PASSWORD" psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -v ON_ERROR_STOP=1 -f "$f" || {
        echo "ERROR in $f"
        exit 1
    }
done

echo "Merge job complete."

# Create indexes
echo "Running index creation ..."
python create_union_id_indexes.py

# Create supplementary tables
echo "Running supplementary table creation ..."
python create_supplementary_tables.py

echo "All tasks finished successfully."