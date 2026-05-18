TRUNCATE TABLE human."BasicInfo", human."CelltypeInfo", human."SampleInfo", mouse."BasicInfo", mouse."CelltypeInfo", mouse."SampleInfo", public.hg38gene, public.mm10gene, public.human_sample_summary, public.mouse_sample_summary;

INSERT INTO human."BasicInfo"
("Union ID", "Loci", "Motif", "Genomic feature", "Cell type gain", "Cell type lost", "Constitutive", "Occupancy score", "Occupancy frequency")
VALUES
(101, 'chr1:100-200', 'CTCF', 'MYC promoter', 'Bcell', '', 'no', 8.5, 0.75);

INSERT INTO human."CelltypeInfo"
("Union ID", "Celltype", "Sample size", "Occupancy frequency in cell type dataset", "Average RPKM (cell type)", "Average RPKM (others)", "-log10(FDR)")
VALUES
(101, 'Bcell', 5, 0.8, 12.5, 3.1, 4.2);

INSERT INTO human."SampleInfo"
("Union ID", "GSM", "Occupancy", "RPKM")
VALUES
(101, 'GSM100001', 1, 10.5);

INSERT INTO mouse."BasicInfo"
("Union ID", "Loci", "Motif", "Genomic feature", "Cell type gain", "Cell type lost", "Constitutive", "Occupancy score", "Occupancy frequency")
VALUES
(201, 'chr1:300-420', 'CTCF', 'Myc promoter', 'Neuron', '', 'no', 7.2, 0.60);

INSERT INTO mouse."CelltypeInfo"
("Union ID", "Celltype", "Sample size", "Occupancy frequency in cell type dataset", "Average RPKM (cell type)", "Average RPKM (others)", "-log10(FDR)")
VALUES
(201, 'Neuron', 4, 0.7, 9.5, 2.0, 3.5);

INSERT INTO mouse."SampleInfo"
("Union ID", "GSM", "Occupancy", "RPKM")
VALUES
(201, 'GSM200001', 1, 8.2);

INSERT INTO public.hg38gene VALUES ('chr1', 120, 180, 'MYC', '+');
INSERT INTO public.mm10gene VALUES ('chr1', 320, 380, 'Myc', '+');

INSERT INTO public.human_sample_summary
("GSM", "GSE", "In Situ Site", "Label", "PubMed ID", "Release Date", "Num Peaks", "FRiP", "motif_ov_ratio")
VALUES
('GSM100001', 'GSETEST1', 'NA', 'Bcell', 'PMID1', '2026-01-01', '1', '0.5', '0.8');

INSERT INTO public.mouse_sample_summary
("GSM", "GSE", "In Situ Site", "Label", "PubMed ID", "Release Date", "Num Peaks", "FRiP", "motif_ov_ratio")
VALUES
('GSM200001', 'GSETEST2', 'NA', 'Neuron', 'PMID2', '2026-01-01', '1', '0.5', '0.8');
