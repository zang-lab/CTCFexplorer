DROP SCHEMA IF EXISTS human CASCADE;
DROP SCHEMA IF EXISTS mouse CASCADE;
CREATE SCHEMA human;
CREATE SCHEMA mouse;

DROP TABLE IF EXISTS public.hg38gene;
DROP TABLE IF EXISTS public.mm10gene;
DROP TABLE IF EXISTS public.human_sample_summary;
DROP TABLE IF EXISTS public.mouse_sample_summary;

CREATE TABLE human."BasicInfo" (
    "Union ID" INTEGER PRIMARY KEY,
    "Loci" TEXT,
    "Motif" TEXT,
    "Genomic feature" TEXT,
    "Cell type gain" TEXT,
    "Cell type lost" TEXT,
    "Constitutive" TEXT,
    "Occupancy score" DOUBLE PRECISION,
    "Occupancy frequency" DOUBLE PRECISION
);

CREATE TABLE human."CelltypeInfo" (
    "Union ID" INTEGER,
    "Celltype" TEXT,
    "Sample size" INTEGER,
    "Occupancy frequency in cell type dataset" DOUBLE PRECISION,
    "Average RPKM (cell type)" DOUBLE PRECISION,
    "Average RPKM (others)" DOUBLE PRECISION,
    "-log10(FDR)" DOUBLE PRECISION
);

CREATE TABLE human."SampleInfo" (
    "Union ID" INTEGER,
    "GSM" TEXT,
    "Occupancy" DOUBLE PRECISION,
    "RPKM" DOUBLE PRECISION
);

CREATE TABLE mouse."BasicInfo" (LIKE human."BasicInfo" INCLUDING ALL);
CREATE TABLE mouse."CelltypeInfo" (LIKE human."CelltypeInfo" INCLUDING ALL);
CREATE TABLE mouse."SampleInfo" (LIKE human."SampleInfo" INCLUDING ALL);

CREATE TABLE public.hg38gene (
    "Chromosome" TEXT,
    "Start" INTEGER,
    "End" INTEGER,
    "Gene" TEXT,
    "Strand" TEXT
);

CREATE TABLE public.mm10gene (LIKE public.hg38gene INCLUDING ALL);

CREATE TABLE public.human_sample_summary (
    "GSM" TEXT,
    "GSE" TEXT,
    "In Situ Site" TEXT,
    "Label" TEXT,
    "PubMed ID" TEXT,
    "Release Date" TEXT,
    "Num Peaks" TEXT,
    "FRiP" TEXT,
    "motif_ov_ratio" TEXT
);

CREATE TABLE public.mouse_sample_summary (LIKE public.human_sample_summary INCLUDING ALL);
