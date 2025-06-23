#!/bin/bash
#SBATCH --job-name=ctcf_merge
#SBATCH --partition=standard
#SBATCH --account=cphg_cz3d
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00

set -euo pipefail    # abort on any error, undefined var, or failed pipeline

# ── Connection to master PostgreSQL ──────────────────────────────
PGDATA=$(ps -ef | grep '[p]ostgres' | grep ' -D ' | awk '{for (i=1;i<=NF;i++) if ($i=="-D") print $(i+1)}' | head -n1)

if [ -z "$PGDATA" ]; then
    echo "Error: Could not auto-detect PGDATA. Is PostgreSQL running?"
    exit 1
fi

export PGHOST=$PGDATA
export PGPORT=5432
export PGDATABASE=${USER}
export PGUSER=${USER}

echo "Detected PGDATA: $PGDATA"
echo "Connecting as user: $PGUSER to database: $PGDATABASE on port: $PGPORT"
# ───────────────────────────────────────────────────────────────────────

# ── Ensure the master server is running (start it if it’s not) ────────
if ! pg_ctl -D "$PGDATA" status >/dev/null 2>&1; then
    echo "Starting master PostgreSQL …"
    pg_ctl -D "$PGDATA" \
           -o "-p $PGPORT -k $PGHOST" \
           -l "$PGDATA/postgres.log" start
    sleep 2
fi
# ───────────────────────────────────────────────────────────────────────

DUMP_DIR=/project/zanglab_project/vdj6tq/project_data/CTCFexplorer/hs/pg_dumps
NUM_DUMPS=$(ls -1 "$DUMP_DIR"/dump_*.sql 2>/dev/null | wc -l)

if [[ $NUM_DUMPS -eq 0 ]]; then
    echo "ERROR: no dump_*.sql files found in $DUMP_DIR"
    exit 1
fi

echo "Replaying $NUM_DUMPS dump files …"

# replay in sorted order to be deterministic
for f in $(ls "$DUMP_DIR"/dump_*.sql | sort -V); do
    echo ">> $(basename "$f")"
    psql -v ON_ERROR_STOP=1 -f "$f"
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