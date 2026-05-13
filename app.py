from flask import Flask, render_template, request, Response, abort, send_file, redirect, url_for
from io import StringIO
import csv
import os
import re
from pathlib import Path

import psycopg2
import psycopg2.extras
from psycopg2 import sql


VALID_SPECIES = {"human", "mouse"}
MAX_UNION_ID_BY_SPECIES = {"human": 531851, "mouse": 297825}

# hg38 chromosome sizes
HG38_CHROM_SIZES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415, "chrM": 16569,
}

# mm10 chromosome sizes
MM10_CHROM_SIZES = {
    "chr1": 195471971, "chr2": 182113224, "chr3": 160039680, "chr4": 156508116,
    "chr5": 151834684, "chr6": 149736546, "chr7": 145441459, "chr8": 129401213,
    "chr9": 124595110, "chr10": 130694993, "chr11": 122082543, "chr12": 120129022,
    "chr13": 120421639, "chr14": 124902244, "chr15": 104043685, "chr16": 98207768,
    "chr17": 94987271, "chr18": 90702639, "chr19": 61431566, "chrX": 171031299,
    "chrY": 91744698, "chrM": 16299,
}

CHROM_SIZES_BY_SPECIES = {
    "human": HG38_CHROM_SIZES,
    "mouse": MM10_CHROM_SIZES,
}

DEFAULT_GENE_TABLES = {
    "human": "public.hg38gene",
    "mouse": "public.mm10gene",
}

DEFAULT_LABEL_TABLES = {
    "human": "public.human_sample_summary",
    "mouse": "public.mouse_sample_summary",
}

GENE_TABLES = {
    "human": os.environ.get("HUMAN_GENE_TABLE", DEFAULT_GENE_TABLES["human"]),
    "mouse": os.environ.get("MOUSE_GENE_TABLE", DEFAULT_GENE_TABLES["mouse"]),
}

LABEL_TABLES = {
    "human": os.environ.get("HUMAN_LABEL_TABLE", DEFAULT_LABEL_TABLES["human"]),
    "mouse": os.environ.get("MOUSE_LABEL_TABLE", DEFAULT_LABEL_TABLES["mouse"]),
}

BASE_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = Path(os.environ.get("TEMPLATES_DIR", str(BASE_DIR / "templates")))
STATIC_DIR = Path(os.environ.get("STATIC_DIR", str(BASE_DIR / "static")))
GENOME_DIR = Path(os.environ.get("GENOME_DIR", str(STATIC_DIR / "data" / "genome")))

app = Flask(__name__, template_folder=str(TEMPLATES_DIR), static_folder=str(STATIC_DIR))


def get_db_connection():
    conn = psycopg2.connect(
        dbname=os.environ.get("DB_NAME", "CTCFDB_PostgreSQL"),
        user=os.environ.get("DB_USER", "postgres"),
        password=os.environ.get("DB_PASSWORD", ""),
        host=os.environ.get("DB_HOST", "localhost"),
        port=os.environ.get("DB_PORT", 5432),
    )
    conn.autocommit = True
    return conn


def normalize_species(species):
    value = (species or "").strip().lower()
    if value not in VALID_SPECIES:
        abort(400, description="Invalid species. Use 'human' or 'mouse'.")
    return value


def parse_fq_table(name, default_schema="public"):
    if not name:
        return None
    cleaned = name.strip()
    if "." in cleaned:
        schema_name, table_name = cleaned.split(".", 1)
        return schema_name.strip('"'), table_name.strip('"')
    return default_schema, cleaned.strip('"')


def table_exists(conn, schema_name, table_name):
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT EXISTS (
              SELECT 1
              FROM information_schema.tables
              WHERE table_schema = %s AND table_name = %s
            )
            """,
            (schema_name, table_name),
        )
        return bool(cur.fetchone()[0])


def query_with_dict_cursor(conn, query, params=()):
    with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
        cur.execute(query, params)
        return cur.fetchall()


def parse_loci(text):
    match = re.match(r"^(chr[0-9XYM]+):(\d+)-(\d+)$", (text or "").strip())
    if not match:
        return None
    chrom, start, end = match.groups()
    return chrom, int(start), int(end)


def valid_loci(species, chrom, start, end):
    sizes = CHROM_SIZES_BY_SPECIES[species]
    if chrom not in sizes:
        return False, f"Invalid chromosome '{chrom}' for {species}."
    if start < 1 or end > sizes[chrom] or start >= end:
        return False, f"Invalid loci range for {chrom}. Valid range is 1-{sizes[chrom]}, with start < end."
    return True, ""


def overlaps_region(loci, search_chr, search_start, search_end):
    match = parse_loci(loci)
    if not match:
        return False
    loci_chr, loci_start, loci_end = match
    return loci_chr == search_chr and not (search_end < loci_start or search_start > loci_end)


def split_celltypes(field_value):
    if not field_value:
        return []
    return [token.strip().lower() for token in re.split(r"[;,]\s*|\s+", str(field_value)) if token.strip()]


def feature_contains_gene(feature_text, gene_symbol):
    if not feature_text:
        return False
    pattern = r"(?<![A-Za-z0-9_-])" + re.escape(gene_symbol) + r"(?![A-Za-z0-9_-])"
    return re.search(pattern, str(feature_text), flags=re.IGNORECASE) is not None


@app.route("/")
def index():
    return render_template("index.html", species=request.args.get("species", ""))


@app.route("/search_loci", methods=["POST"])
def search_loci():
    species = normalize_species(request.form.get("species"))
    loci_input = request.form.get("loci", "").strip()
    parsed = parse_loci(loci_input)
    if not parsed:
        return render_template("loci_error.html", species=species, error="Use format chr#:start-end, e.g. chr1:10000-20000.")

    search_chr, search_start, search_end = parsed
    valid, error = valid_loci(species, search_chr, search_start, search_end)
    if not valid:
        return render_template("loci_error.html", species=species, error=error)

    conn = get_db_connection()
    try:
        query = sql.SQL(
            """
            SELECT "Union ID" AS union_id, "Loci", "Occupancy score" AS occupancy_score,
                   "Occupancy frequency" AS occupancy_frequency, "Motif"
            FROM {}."BasicInfo"
            WHERE "Loci" LIKE %s
            """
        ).format(sql.Identifier(species))
        rows = query_with_dict_cursor(conn, query, (f"{search_chr}:%",))
    finally:
        conn.close()

    results = [dict(row) for row in rows if overlaps_region(row.get("Loci"), search_chr, search_start, search_end)]
    results.sort(key=lambda row: int(row.get("union_id")) if str(row.get("union_id", "")).isdigit() else float("inf"))
    search_region = {"chr": search_chr, "start": search_start, "end": search_end}
    if not results:
        return render_template(
            "loci_not_found.html",
            species=species,
            search_loci=loci_input,
            search_region=search_region,
            error=f"No union binding site found for {loci_input}.",
        )
    return render_template("results_Loci.html", species=species, results=results, search_region=search_region, query=loci_input)


@app.route("/search_union", methods=["POST"])
def search_union():
    species = normalize_species(request.form.get("species"))
    union_text = request.form.get("union", "").strip()
    if not union_text.isdigit():
        return render_template("union_not_found.html", species=species, search_loci=union_text, error="Union ID must be numeric.")

    union_id = int(union_text)
    max_union_id = MAX_UNION_ID_BY_SPECIES[species]
    if union_id < 1 or union_id > max_union_id:
        return render_template(
            "union_not_found.html",
            species=species,
            search_loci=union_id,
            error=(
                f"Union ID {union_id} is out of range. "
                "Valid Union ID range is 1 to 531851 for human and 1 to 297825 for mouse."
            ),
        )

    conn = get_db_connection()
    try:
        with conn.cursor() as cur:
            query = sql.SQL('SELECT 1 FROM {}."BasicInfo" WHERE "Union ID" = %s LIMIT 1').format(sql.Identifier(species))
            cur.execute(query, (union_id,))
            exists = cur.fetchone() is not None
    finally:
        conn.close()

    if exists:
        return redirect(url_for("union_info", species=species, union_id=union_id))
    return render_template("union_not_found.html", species=species, search_loci=union_id, error=f"Union ID {union_id} was not found in {species}.")


@app.route("/union_info/<species>/<int:union_id>")
def union_info(species, union_id):
    species = normalize_species(species)
    conn = get_db_connection()
    try:
        basic_q = sql.SQL(
            '''
            SELECT "Union ID", "Loci", "Motif", "Genomic feature", "Cell type gain", "Cell type lost",
                   "Constitutive", "Occupancy score", "Occupancy frequency"
            FROM {}."BasicInfo"
            WHERE "Union ID" = %s
            '''
        ).format(sql.Identifier(species))
        cell_q = sql.SQL(
            '''
            SELECT "Celltype", "Sample size", "Occupancy frequency in cell type dataset",
                   "Average RPKM (cell type)", "Average RPKM (others)", "-log10(FDR)"
            FROM {}."CelltypeInfo"
            WHERE "Union ID" = %s
            ORDER BY "Occupancy frequency in cell type dataset" DESC
            '''
        ).format(sql.Identifier(species))
        sample_q = sql.SQL(
            '''
            SELECT "GSM", "Occupancy", "RPKM"
            FROM {}."SampleInfo"
            WHERE "Union ID" = %s
            ORDER BY "Occupancy" DESC
            '''
        ).format(sql.Identifier(species))

        basic = query_with_dict_cursor(conn, basic_q, (union_id,))
        celltype = query_with_dict_cursor(conn, cell_q, (union_id,))
        sample = query_with_dict_cursor(conn, sample_q, (union_id,))
    finally:
        conn.close()

    if not basic:
        return render_template("union_not_found.html", species=species, search_loci=union_id, error=f"Union ID {union_id} was not found in {species}.")

    return render_template(
        "union_info.html",
        species=species,
        union_id=union_id,
        basic=[dict(row) for row in basic],
        celltype=[dict(row) for row in celltype],
        sample=[dict(row) for row in sample],
        results=[dict(row) for row in basic],
    )


@app.route("/search_celltype", methods=["GET", "POST"])
def search_celltype():
    species = normalize_species((request.form.get("species") if request.method == "POST" else request.args.get("species")) or "human")
    celltype = (request.form.get("celltype") if request.method == "POST" else request.args.get("celltype", "")).strip()
    if not celltype:
        return render_template("celltype_not_found.html", species=species, celltype=celltype)

    normalized_input = celltype.lower()
    token_pattern = r"(^|[ ,;])" + re.escape(normalized_input) + r"([ ,;]|$)"

    normalized_celltype = celltype
    conn = get_db_connection()
    try:
        # Resolve display case using the species metadata table when available.
        label_cfg = parse_fq_table(LABEL_TABLES[species])
        if label_cfg and table_exists(conn, label_cfg[0], label_cfg[1]):
            canonical_q = sql.SQL(
                '''
                SELECT "Label"
                FROM {}.{}
                WHERE LOWER("Label") = LOWER(%s)
                LIMIT 1
                '''
            ).format(sql.Identifier(label_cfg[0]), sql.Identifier(label_cfg[1]))
            canonical_rows = query_with_dict_cursor(conn, canonical_q, (celltype,))
            if canonical_rows:
                normalized_celltype = canonical_rows[0]["Label"]

        gain_q = sql.SQL('SELECT * FROM {}."BasicInfo" WHERE COALESCE("Cell type gain", \'\') ~* %s').format(sql.Identifier(species))
        loss_q = sql.SQL('SELECT * FROM {}."BasicInfo" WHERE COALESCE("Cell type lost", \'\') ~* %s').format(sql.Identifier(species))
        gain_rows = query_with_dict_cursor(conn, gain_q, (token_pattern,))
        loss_rows = query_with_dict_cursor(conn, loss_q, (token_pattern,))

        gsm_rows = []
        if label_cfg and table_exists(conn, label_cfg[0], label_cfg[1]):
            label_q = sql.SQL('SELECT * FROM {}.{} WHERE LOWER("Label") = LOWER(%s)').format(
                sql.Identifier(label_cfg[0]), sql.Identifier(label_cfg[1])
            )
            gsm_rows = query_with_dict_cursor(conn, label_q, (normalized_celltype,))
    finally:
        conn.close()

    gain_results = [dict(row) for row in gain_rows]
    loss_results = [dict(row) for row in loss_rows]

    if not gain_results and not loss_results and not gsm_rows:
        return render_template("celltype_not_found.html", species=species, celltype=celltype)

    return render_template(
        "results_Celltype.html",
        species=species,
        celltype=normalized_celltype,
        gsm_results=[dict(row) for row in gsm_rows],
        gain_results=gain_results,
        loss_results=loss_results,
    )


@app.route("/search_gsm", methods=["POST"])
def search_gsm():
    species = normalize_species(request.form.get("species"))
    gsm = request.form.get("gsm", "").strip()

    if not re.match(r"^GSM\d+$", gsm):
        return render_template("gsm_not_found.html", species=species, error="Invalid GSM format. Ensure the input starts with 'GSM' followed by numbers.", is_valid_gsm=False)

    conn = get_db_connection()
    metadata = []
    try:
        label_cfg = parse_fq_table(LABEL_TABLES[species])
        if label_cfg and table_exists(conn, label_cfg[0], label_cfg[1]):
            label_q = sql.SQL('SELECT * FROM {}.{} WHERE "GSM" = %s').format(
                sql.Identifier(label_cfg[0]), sql.Identifier(label_cfg[1])
            )
            metadata = query_with_dict_cursor(conn, label_q, (gsm,))
    finally:
        conn.close()

    if not metadata:
        return render_template("gsm_not_found.html", species=species, error="", is_valid_gsm=True, gsm=gsm)

    gsm_payload = {"GSM": gsm, "GSE": "", "In Situ Site": "", "Label": "", "PubMed ID": "", "Release Date": "", "Num Peaks": "", "FRiP": "", "motif_ov_ratio": ""}
    if metadata:
        gsm_payload.update(dict(metadata[0]))

    return render_template(
        "gsm_info.html",
        species=species,
        gsm=gsm_payload,
    )


@app.route("/search_gene", methods=["POST"])
def search_gene():
    species = normalize_species(request.form.get("species"))
    gene = request.form.get("gene", "").strip()
    if not gene:
        return render_template("gene_not_found.html", species=species, search_loci=gene, error="Please provide a gene symbol.")

    try:
        extension_window = int(request.form.get("window", 5000))
    except ValueError:
        extension_window = 5000

    conn = get_db_connection()
    try:
        gene_table_cfg = parse_fq_table(GENE_TABLES[species])
        if gene_table_cfg and table_exists(conn, gene_table_cfg[0], gene_table_cfg[1]):
            query = sql.SQL(
                '''
                SELECT "Chromosome", "Start", "End"
                FROM {}.{}
                WHERE LOWER("Gene") = LOWER(%s)
                LIMIT 1
                '''
            ).format(sql.Identifier(gene_table_cfg[0]), sql.Identifier(gene_table_cfg[1]))
            gene_rows = query_with_dict_cursor(conn, query, (gene,))

            if not gene_rows:
                return render_template("gene_not_found.html", species=species, search_loci=gene, error=f"Gene '{gene}' was not found in the configured gene table.")

            gene_row = gene_rows[0]
            search_chr = gene_row["Chromosome"]
            gene_start = int(gene_row["Start"])
            gene_end = int(gene_row["End"])

            sizes = CHROM_SIZES_BY_SPECIES[species]
            if search_chr not in sizes:
                return render_template("loci_not_found.html", species=species, search_loci=gene, error=f"Chromosome '{search_chr}' in gene table is invalid for {species}.", search_region=None)

            search_start = max(1, gene_start - extension_window)
            search_end = min(sizes[search_chr], gene_end + extension_window)

            basic_q = sql.SQL(
                '''
                SELECT "Union ID" AS union_id, "Loci", "Occupancy score" AS occupancy_score,
                       "Occupancy frequency" AS occupancy_frequency, "Motif", "Genomic feature"
                FROM {}."BasicInfo"
                WHERE "Loci" LIKE %s
                '''
            ).format(sql.Identifier(species))
            rows = query_with_dict_cursor(conn, basic_q, (f"{search_chr}:%",))
            results = [dict(row) for row in rows if overlaps_region(row.get("Loci"), search_chr, search_start, search_end)]

            return render_template(
                "results_Gene.html",
                species=species,
                search_gene=gene,
                results=results,
                no_overlaps=(len(results) == 0),
                search_region={"chr": search_chr, "start": search_start, "end": search_end},
                extension_window=extension_window,
            )

        # Fallback mode for incomplete deployments: match by genomic feature text only.
        basic_q = sql.SQL(
            '''
            SELECT "Union ID" AS union_id, "Loci", "Occupancy score" AS occupancy_score,
                   "Occupancy frequency" AS occupancy_frequency, "Motif", "Genomic feature"
            FROM {}."BasicInfo"
            WHERE "Genomic feature" ILIKE %s
            '''
        ).format(sql.Identifier(species))
        rows = query_with_dict_cursor(conn, basic_q, (f"%{gene}%",))
        results = [dict(row) for row in rows if feature_contains_gene(row.get("Genomic feature"), gene)]
    finally:
        conn.close()

    fallback_note = (
        f"Gene coordinate table is not configured for {species}. "
        "Results below are feature-text matches from BasicInfo only."
    )
    if results:
        parsed = parse_loci(results[0].get("Loci", ""))
        if parsed:
            search_region = {"chr": parsed[0], "start": parsed[1], "end": parsed[2]}
        else:
            search_region = {"chr": "chr1", "start": 1, "end": 10000}
    else:
        search_region = {"chr": "chr1", "start": 1, "end": 10000}
    return render_template(
        "results_Gene.html",
        species=species,
        search_gene=gene,
        results=results,
        no_overlaps=(len(results) == 0),
        search_region=search_region,
        extension_window=extension_window,
        fallback_note=fallback_note,
    )


@app.route("/download_table/<table_name>/<identifier>")
def download_table_legacy(table_name, identifier):
    species = normalize_species(request.args.get("species", "human"))
    return download_table(species, table_name, identifier)


@app.route("/download_table/<species>/<table_name>/<identifier>")
def download_table(species, table_name, identifier):
    species = normalize_species(species)

    queries = {
        "basic_info": sql.SQL('SELECT * FROM {}."BasicInfo" WHERE "Union ID" = %s').format(sql.Identifier(species)),
        "celltype_info": sql.SQL(
            '''
            SELECT "Celltype", "Sample size", "Occupancy frequency in cell type dataset",
                   "Average RPKM (cell type)", "Average RPKM (others)", "-log10(FDR)"
            FROM {}."CelltypeInfo"
            WHERE "Union ID" = %s
            ORDER BY "Occupancy frequency in cell type dataset" DESC
            '''
        ).format(sql.Identifier(species)),
        "sample_info": sql.SQL(
            '''
            SELECT "GSM", "Occupancy", "RPKM"
            FROM {}."SampleInfo"
            WHERE "Union ID" = %s
            ORDER BY "Occupancy" DESC
            '''
        ).format(sql.Identifier(species)),
        "gain_results": sql.SQL('SELECT * FROM {}."BasicInfo" WHERE "Cell type gain" ILIKE %s').format(sql.Identifier(species)),
        "loss_results": sql.SQL('SELECT * FROM {}."BasicInfo" WHERE "Cell type lost" ILIKE %s').format(sql.Identifier(species)),
    }

    label_cfg = parse_fq_table(LABEL_TABLES[species])
    if table_name == "gsm_results":
        if not label_cfg:
            abort(404, description=f"No label metadata table configured for {species}.")
        queries["gsm_results"] = sql.SQL('SELECT * FROM {}.{} WHERE LOWER("Label") = LOWER(%s)').format(
            sql.Identifier(label_cfg[0]), sql.Identifier(label_cfg[1])
        )

    if table_name not in queries:
        abort(400, description=f"Invalid table name: {table_name}")

    if table_name in {"basic_info", "celltype_info", "sample_info"}:
        if not identifier.isdigit():
            abort(400, description="Identifier for union-based downloads must be numeric.")
        param = int(identifier)
    elif table_name == "gsm_results":
        param = identifier
    else:
        param = f"%{identifier}%"

    conn = get_db_connection()
    try:
        with conn.cursor() as cur:
            cur.execute(queries[table_name], (param,))
            rows = cur.fetchall()
            headers = [desc[0] for desc in cur.description]
    finally:
        conn.close()

    if not rows:
        abort(404, description=f"No data found for {identifier} in {species}/{table_name}.")

    output = StringIO()
    writer = csv.writer(output)
    writer.writerow(headers)
    writer.writerows(rows)

    filename = f"{species}_{table_name}_{identifier}.csv"
    return Response(output.getvalue(), mimetype="text/csv", headers={"Content-Disposition": f"attachment; filename={filename}"})


@app.route("/download_file/<folder>/<filename>")
def download_file(folder, filename):
    file_path = os.path.join(app.static_folder, "data", folder, filename)
    if not os.path.exists(file_path):
        abort(404, description=f"File {filename} not found in {folder}.")
    return send_file(file_path, as_attachment=True, download_name=filename)


@app.route("/genome_file/<filename>")
def genome_file(filename):
    if filename not in {"hg38.fa", "hg38.fa.fai", "mm10.fa", "mm10.fa.fai"}:
        abort(404, description="Genome file not found.")
    file_path = GENOME_DIR / filename
    if not file_path.exists():
        abort(404, description=f"Genome file not found: {filename}")
    return send_file(file_path, as_attachment=False, download_name=filename)


@app.route("/annotation_file/<filename>")
def annotation_file(filename):
    allowed = {
        "hg38.refGene.transcripts.bed",
        "mm10.ncbiRefSeq.transcripts.bed",
    }
    if filename not in allowed:
        abort(404, description="Annotation file not found.")
    file_path = GENOME_DIR / filename
    if not file_path.exists():
        abort(404, description=f"Annotation file not found: {filename}")
    return send_file(file_path, as_attachment=False, conditional=True)


@app.route("/available_celltypes")
def available_celltypes():
    # Keep supporting species query parameter for backward compatibility,
    # but this page now shows both human and mouse sections together.
    _ = request.args.get("species", "human")
    conn = get_db_connection()
    try:
        def get_celltypes_for_species(species_name):
            label_cfg = parse_fq_table(LABEL_TABLES[species_name])
            if label_cfg and table_exists(conn, label_cfg[0], label_cfg[1]):
                query = sql.SQL(
                    '''
                    SELECT label
                    FROM (
                        SELECT DISTINCT "Label" AS label
                        FROM {}.{}
                        WHERE "Label" IS NOT NULL AND BTRIM("Label") <> ''
                    ) t
                    ORDER BY LOWER(label)
                    '''
                ).format(sql.Identifier(label_cfg[0]), sql.Identifier(label_cfg[1]))
                with conn.cursor() as cur:
                    cur.execute(query)
                    return [row[0] for row in cur.fetchall()]

            fallback_query = sql.SQL(
                '''
                SELECT DISTINCT "Celltype"
                FROM {}."CelltypeInfo"
                WHERE "Celltype" IS NOT NULL AND BTRIM("Celltype") <> ''
                ORDER BY LOWER("Celltype")
                '''
            ).format(sql.Identifier(species_name))
            with conn.cursor() as cur:
                cur.execute(fallback_query)
                return [row[0] for row in cur.fetchall()]

        human_celltypes = get_celltypes_for_species("human")
        mouse_celltypes = get_celltypes_for_species("mouse")
    finally:
        conn.close()

    return render_template(
        "available_celltypes.html",
        human_celltypes=human_celltypes,
        mouse_celltypes=mouse_celltypes,
    )


@app.route("/contact")
def contact():
    return render_template("contact.html")


@app.route("/help")
def help_page():
    return render_template("help.html")


@app.route("/download")
def download_page():
    return render_template("download.html")


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=os.environ.get("FLASK_ENV") == "development")
