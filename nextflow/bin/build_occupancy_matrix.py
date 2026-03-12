#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
CHROM_RANK = {chrom: idx for idx, chrom in enumerate(CHROM_ORDER)}


def read_union_bed(path):
    regions = {}
    union_rows = []
    with open(path) as handle:
        for line in handle:
            chrom, start, end, union_id, score, strand = line.rstrip("\n").split("\t")[:6]
            start = int(start)
            end = int(end)
            regions.setdefault(chrom, []).append((start, end, union_id))
            union_rows.append((chrom, start, end, union_id, score, strand))
    return regions, union_rows


def read_peak_windows(path, width, fe_cutoff):
    half = width // 2
    by_chrom = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom = fields[0]
            if chrom not in CHROM_RANK:
                continue
            start = int(fields[1])
            enrichment = float(fields[6])
            summit_offset = int(fields[9])
            if summit_offset < 0 or enrichment < fe_cutoff:
                continue
            summit = start + summit_offset
            by_chrom.setdefault(chrom, []).append((max(0, summit - half), summit + half))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


def overlaps(union_regions, sample_regions):
    result = {}
    for chrom, union_intervals in union_regions.items():
        peaks = sample_regions.get(chrom, [])
        j = 0
        for start, end, union_id in union_intervals:
            while j < len(peaks) and peaks[j][1] < start:
                j += 1
            k = j
            found = 0
            while k < len(peaks) and peaks[k][0] <= end:
                if peaks[k][1] >= start:
                    found = 1
                    break
                k += 1
            result[union_id] = found
    return result


def sample_name(path_str):
    name = Path(path_str).name
    return name.replace("_peaks.narrowPeak", "")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--union", required=True)
    parser.add_argument("--peaks", nargs="+", required=True)
    parser.add_argument("--output-matrix", required=True)
    parser.add_argument("--output-scores", required=True)
    parser.add_argument("--fe-cutoff", type=float, default=4.0)
    parser.add_argument("--sample-width", type=int, default=50)
    args = parser.parse_args()

    union_regions, union_rows = read_union_bed(args.union)
    matrix = pd.DataFrame(index=[row[3] for row in union_rows])

    for peak_file in args.peaks:
        matrix[sample_name(peak_file)] = pd.Series(
            overlaps(union_regions, read_peak_windows(peak_file, args.sample_width, args.fe_cutoff))
        )

    matrix = matrix.fillna(0).astype(int)
    matrix.index.name = "Union_ID"
    matrix.to_csv(args.output_matrix, sep="\t")

    scores = []
    occupancy_sum = matrix.sum(axis=1)
    for chrom, start, end, union_id, score, strand in union_rows:
        scores.append(
            {
                "chr": chrom,
                "start": start,
                "end": end,
                "Union_ID": union_id,
                "score": score,
                "strand": strand,
                "occupancy_score": int(occupancy_sum.loc[union_id]),
            }
        )
    pd.DataFrame(scores).to_csv(args.output_scores, sep="\t", index=False)


if __name__ == "__main__":
    main()
