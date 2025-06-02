from flask import Flask, render_template, request, Response, abort, send_file, redirect, url_for
import re
import os
import csv
from io import StringIO
import psycopg2
import psycopg2.extras


# Define paths
TEMPLATES_PATH = "templates"
STATIC_PATH = "static"

# Dictionary of valid hg38 chromosome sizes
HG38_CHROM_SIZES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415
}

# Initialize Flask app with the templates folder path
app = Flask(__name__, template_folder=TEMPLATES_PATH, static_folder=STATIC_PATH)

def get_db_connection():
    conn = psycopg2.connect(
        dbname=os.environ.get("DB_NAME", "CTCFDB_PostgreSQL"),
        user=os.environ.get("DB_USER", "postgres"),
        password=os.environ.get("DB_PASSWORD", ""),  # must be set in deployment
        host=os.environ.get("DB_HOST", "localhost"),
        port=os.environ.get("DB_PORT", 5432)
    )
    conn.autocommit = True
    return conn


# Route for the home page
@app.route('/')
def index():
    return render_template('index.html')

# Route for searching loci
@app.route('/search_loci', methods=['POST'])
def search_loci():
    search_loci = request.form['loci']

    # Validate loci format
    match = re.match(r"(chr[\dXY]+):(\d+)-(\d+)", search_loci)
    if not match:
        return render_template(
            'loci_error.html',
            error="Invalid loci format. Use the format chr#:start-end (e.g., chr1:10000-20000)."
        )

    search_chr, search_start, search_end = match.groups()
    search_start, search_end = int(search_start), int(search_end)

    # Check if chromosome is valid for hg38
    if search_chr not in HG38_CHROM_SIZES:
        return render_template(
            'loci_error.html',
            error=f"Invalid chromosome: {search_chr}. Please use a valid hg38 chromosome."
        )

    # Validate start and end positions
    chr_length = HG38_CHROM_SIZES[search_chr]
    if search_start < 1 or search_end > chr_length or search_start >= search_end:
        return render_template(
            'loci_error.html',
            error=f"Invalid loci range for {search_chr}. Ensure the range is within 1-{chr_length} and start < end."
        )

    conn = get_db_connection()
    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    # Query the database
    cursor.execute("""
        SELECT "Union ID" AS union_id, "Loci", "Occupancy score" AS occupancy_score, "Occupancy frequency" AS occupancy_frequency, "Motif"
        FROM "BasicInfo"
        WHERE "Loci" LIKE %s
    """, (f"{search_chr}:%",))
    results = cursor.fetchall()
    conn.close()

    # Debugging: Log the results
    print("Raw Results from Database:", results)

    # Check for overlaps
    overlapping_results = []
    for row in results:
        print("Processing Row:", row)  # Debug: Inspect each row
        loci_match = re.match(r"(chr[\dXY]+):(\d+)-(\d+)", row['Loci'])
        if loci_match:
            loci_chr, loci_start, loci_end = loci_match.groups()
            loci_start, loci_end = int(loci_start), int(loci_end)
            if loci_chr == search_chr and not (search_end < loci_start or search_start > loci_end):
                overlapping_results.append(dict(row))

    # Debug: Inspect overlapping results
    print("Overlapping Results:", overlapping_results)

    # Define search_region
    search_region = {"chr": search_chr, "start": search_start, "end": search_end}

    if overlapping_results:
        return render_template('results_Loci.html', results=overlapping_results, search_region=search_region)
    else:
        return render_template(
            'loci_not_found.html',
            search_loci=search_loci,
            search_region=search_region,  # Pass the search region to display IGV
            error=f"No union binding site found for {search_loci}. You may try increasing the search range or check for potential errors in the input format."
        )

# Route to search by GSM
@app.route('/search_gsm', methods=['POST'])
def search_gsm():
    gsm = request.form['gsm']

    # Validate GSM format
    if not re.match(r'^GSM\d+$', gsm):
        return render_template(
            'gsm_not_found.html',
            error="Invalid GSM format. Ensure the input starts with 'GSM' followed by numbers.",
            is_valid_gsm=False
        )

    # Search GSM in the database
    conn = get_db_connection()
    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    cursor.execute('SELECT * FROM "CTCFLabels" WHERE "GSM" = %s', (gsm,))
    result = cursor.fetchone()
    conn.close()

    if not result:
        return render_template(
            'gsm_not_found.html',
            error="",
            is_valid_gsm=True,
            gsm=gsm
        )

    return render_template(
        'gsm_info.html',
        gsm=result  # result is already a dict-like object
    )

# Route for searching by cell type
@app.route('/search_celltype', methods=['GET', 'POST'])
def search_celltype():
    if request.method == 'POST':
        # Handle POST request from a search form
        celltype = request.form['celltype'].strip()
    else:
        # Handle GET request when a celltype link is clicked
        celltype = request.args.get('celltype', '').strip()

    # Normalize the cell type input and query case-insensitively
    normalized_input = re.sub(r'[-_/\\#]', '', celltype).lower()

    conn = get_db_connection()
    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    # Query for normalized cell type
    cursor.execute("""
        SELECT DISTINCT "Label"
        FROM "CTCFLabels"
        WHERE LOWER(REPLACE(REPLACE(REPLACE(REPLACE(REPLACE("Label", '-', ''), '_', ''), '/', ''), '\\', ''), '#', '')) = %s
    """, (normalized_input,))
    normalized_result = cursor.fetchone()

    if not normalized_result:
        conn.close()
        return render_template(
            'celltype_not_found.html',
            celltype=celltype
        )

    normalized_celltype = normalized_result["Label"]

    # Query for GSMs and union bindings
    cursor.execute("""
        SELECT *
        FROM "CTCFLabels"
        WHERE "Label" = %s
    """, (normalized_celltype,))
    gsm_results = cursor.fetchall()

    cursor.execute("""
        SELECT *
        FROM "BasicInfo"
        WHERE "Cell type gain" LIKE %s
    """, (f'%{normalized_celltype}%',))
    gain_results = cursor.fetchall()

    cursor.execute("""
        SELECT *
        FROM "BasicInfo"
        WHERE "Cell type lost" LIKE %s
    """, (f'%{normalized_celltype}%',))
    loss_results = cursor.fetchall()

    conn.close()

    return render_template(
        'results_Celltype.html',
        celltype=normalized_celltype,
        gsm_results=gsm_results,
        gain_results=gain_results,
        loss_results=loss_results
    )

# Route to view details for a specific Union ID
@app.route('/union_info/<union_id>')
def union_info(union_id):
    conn = get_db_connection()
    try:
        cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

        cursor.execute(
            '''
            SELECT "Union ID", "Loci", "Motif", "Genomic feature", "Cell type gain", "Cell type lost",
                   "Constitutive", "Occupancy score", "Occupancy frequency"
            FROM "BasicInfo"
            WHERE "Union ID" = %s
            ''',
            (union_id,)
        )
        basic = cursor.fetchall()

        cursor.execute(
            '''
            SELECT "Celltype", "Sample size", "Occupancy frequency in cell type dataset",
                   "Average RPKM (cell type)", "Average RPKM (others)", "-log10(FDR)"
            FROM "CelltypeInfo"
            WHERE "Union ID" = %s
            ORDER BY "Occupancy frequency in cell type dataset" DESC
            ''',
            (union_id,)
        )
        celltype = cursor.fetchall()

        cursor.execute(
            '''
            SELECT "GSM", "Occupancy", "RPKM"
            FROM "SampleInfo"
            WHERE "Union ID" = %s
            ORDER BY "Occupancy" DESC
            ''',
            (union_id,)
        )
        sample = cursor.fetchall()

        cursor.close()
    finally:
        conn.close()

    # Render the template with sorted data
    return render_template(
        'union_info.html',
        union_id=union_id,
        basic=[dict(row) for row in basic],
        celltype=[dict(row) for row in celltype],
        sample=[dict(row) for row in sample]
    )

@app.route('/search_gene', methods=['POST'])
@app.route('/search_gene', methods=['POST'])
def search_gene():
    search_gene = request.form['gene'].strip()
    try:
        extension_window = int(request.form.get('window', 5000))  # Default to 5000 if no input
    except ValueError:
        extension_window = 5000

    # Connect to the database
    conn = get_db_connection()
    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)

    # Query for gene coordinates
    cursor.execute("""
        SELECT "Chromosome", "Start", "End"
        FROM "hg38gene"
        WHERE LOWER("Gene") = LOWER(%s)
    """, (search_gene,))
    gene_result = cursor.fetchone()

    if not gene_result:
        conn.close()
        return render_template(
            'gene_not_found.html',
            search_loci=search_gene,
            error=f"Gene '{search_gene}' not found in the database. Please input a valid human gene."
        )

    # Unpack and validate gene coordinates
    search_chr = gene_result['Chromosome']
    gene_start = int(gene_result['Start'])
    gene_end = int(gene_result['End'])

    search_start = max(1, gene_start - extension_window)
    search_end = gene_end + extension_window

    if search_chr not in HG38_CHROM_SIZES:
        conn.close()
        return render_template(
            'loci_not_found.html',
            search_loci=f"{search_chr}:{gene_start}-{gene_end}",
            error=f"Invalid chromosome: {search_chr}. Please use a valid hg38 chromosome."
        )

    chr_length = HG38_CHROM_SIZES[search_chr]
    if search_end > chr_length:
        search_end = chr_length

    # Query for loci overlaps
    cursor.execute("""
        SELECT "Union ID" AS union_id, "Loci", "Occupancy score" AS occupancy_score, 
               "Occupancy frequency" AS occupancy_frequency, "Motif"
        FROM "BasicInfo"
        WHERE "Loci" LIKE %s
    """, (f"{search_chr}:%",))
    results = cursor.fetchall()
    conn.close()

    # Check for overlaps
    overlapping_results = []
    for row in results:
        loci = row.get('Loci', '')
        if isinstance(loci, str):
            loci_match = re.match(r"(chr[\dXY]+):(\d+)-(\d+)", loci)
        else:
            loci_match = None

        if not loci_match:
            print(f"Skipping invalid or missing Loci: {loci}")
            continue

        loci_chr, loci_start, loci_end = loci_match.groups()
        loci_start, loci_end = int(loci_start), int(loci_end)

        if loci_chr == search_chr and not (search_end < loci_start or search_start > loci_end):
            overlapping_results.append(dict(row))

    # Prepare output region
    search_region = {"chr": search_chr, "start": search_start, "end": search_end}

    return render_template(
        'results_Gene.html',
        results=overlapping_results,
        search_region=search_region,
        search_gene=search_gene,
        extension_window=extension_window,
        no_overlaps=(len(overlapping_results) == 0)
    )

@app.route('/search_union', methods=['POST'])
def search_union():
    search_union = request.form['union'].strip()

    # Ensure input is numeric
    if not search_union.isdigit():
        return render_template(
            'union_not_found.html',
            search_loci=search_union,
            error="Union ID must be a number within the range of 1 to 531851."
        )

    search_union = int(search_union)

    # Check if Union ID is within the valid range
    if not (1 <= search_union <= 531851):
        return render_template(
            'union_not_found.html',
            search_loci=search_union,
            error=f"Union ID '{search_union}' is out of range. "
                  f"Valid Union IDs range from 1 to 531851."
        )

    # Connect to the database and check if Union ID exists
    conn = get_db_connection()
    try:
        cursor = conn.cursor()
        cursor.execute('SELECT "Union ID" FROM "BasicInfo" WHERE "Union ID" = %s LIMIT 1', (search_union,))
        result = cursor.fetchone()
        cursor.close()
    finally:
        conn.close()

    # If Union ID exists, redirect to details page
    if result:
        return redirect(url_for('union_info', union_id=search_union))

    # If not found, show exclusion message
    return render_template(
        'union_not_found.html',
        search_loci=search_union,
        error=f"Union ID '{search_union}' is excluded from the database as it is not a high-confidence binding site."
    )

@app.route('/download_table/<table_name>/<identifier>')
def download_table(table_name, identifier):
    # Map table names to database queries
    queries = {
        "basic_info": 'SELECT * FROM "BasicInfo" WHERE "Union ID" = %s',
        "celltype_info": '''
            SELECT "Celltype", 
                   "Sample size", 
                   "Occupancy frequency in cell type dataset", 
                   "Average RPKM (cell type)", 
                   "Average RPKM (others)", 
                   "-log10(FDR)"
            FROM "CelltypeInfo"
            WHERE "Union ID" = %s
            ORDER BY "Occupancy frequency in cell type dataset" DESC
        ''',
        "sample_info": '''
            SELECT GSM, Occupancy, RPKM
            FROM "SampleInfo"
            WHERE "Union ID" = %s
            ORDER BY Occupancy DESC
        ''',
        "gsm_results": """
            SELECT "GSM", 
                   "GSE", 
                   "Label" AS "Cell Type", 
                   "PubMed ID", 
                   "Release Date", 
                   "In Situ Site", 
                   "FRiP", 
                   "Num Peaks" AS "Number of Peaks", 
                   "motif_ov_ratio" AS "Motif Overlap Ratio"
            FROM "CTCFLabels"
            WHERE "Label" = %s
        """,
        "gain_results": """
            SELECT *
            FROM "BasicInfo"
            WHERE "Cell type gain" LIKE %s
        """,
        "loss_results": """
            SELECT *
            FROM "BasicInfo"
            WHERE "Cell type lost" LIKE %s
        """
    }

    # Check if the table name is valid
    if table_name not in queries:
        abort(400, description=f"Invalid table name: {table_name}")

    # Determine the query parameter
    param = identifier
    if table_name in ["gain_results", "loss_results"]:
        param = f"%{identifier}%"  # Add wildcards for LIKE queries

    # Query the database
    conn = get_db_connection()
    cursor = conn.cursor()
    cursor.execute(queries[table_name], (param,))
    rows = cursor.fetchall()

    # Fetch column headers
    column_headers = [description[0] for description in cursor.description]
    conn.close()

    # If no data found, return a 404 error
    if not rows:
        abort(404, description=f"No data found for identifier: {identifier} in table: {table_name}")

    # Create a CSV file in memory
    def generate_csv():
        output = StringIO()
        writer = csv.writer(output)
        writer.writerow(column_headers)  # Write headers
        writer.writerows(rows)          # Write data rows
        output.seek(0)
        return output.getvalue()

    # Generate a dynamic filename
    if table_name in ["basic_info", "celltype_info", "sample_info"]:
        filename = f"{table_name}_union_{identifier}.csv"
    else:
        filename = f"{table_name}_{identifier}.csv"

    # Return the CSV as a downloadable response
    return Response(
        generate_csv(),
        mimetype='text/csv',
        headers={"Content-Disposition": f"attachment; filename={filename}"}
    )

@app.route('/download_file/<folder>/<filename>')
def download_file(folder, filename):
    # Define the path to the file in the static/data directory
    file_path = os.path.join(STATIC_PATH, 'data', folder, filename)

    # Check if the file exists
    if not os.path.exists(file_path):
        abort(404, description=f"File {filename} not found in {folder} folder.")

    # Send the file for download
    return send_file(
        file_path,
        as_attachment=True,
        download_name=filename
    )

@app.route('/contact')
def contact():
    return render_template('contact.html')

@app.route('/help')
def help():
    return render_template('help.html')

@app.route('/download')
def download():
    return render_template('download.html')

@app.route('/available_celltypes')
def available_celltypes():
    conn = get_db_connection()
    cursor = conn.cursor()

    # Query distinct cell types from the database
    cursor.execute("""
        SELECT DISTINCT "Label"
        FROM "CTCFLabels"
        ORDER BY LOWER("Label") ASC
    """)
    celltypes = [row[0] for row in cursor.fetchall()]
    conn.close()

    return render_template('available_celltypes.html', celltypes=celltypes)


# Run the Flask app
if __name__ == '__main__':
    app.run(debug=True)
