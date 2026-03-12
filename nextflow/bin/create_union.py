#!/usr/bin/env python3

import argparse
from pathlib import Path


CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
CHROM_RANK = {chrom: idx for idx, chrom in enumerate(CHROM_ORDER)}


def iter_summit_windows(path, width, fe_cutoff):
    half = width // 2
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
            yield chrom, max(0, summit - half), summit + half


def merge_intervals(intervals):
    merged = []
    for chrom in CHROM_ORDER:
        chrom_intervals = sorted(intervals.get(chrom, []))
        if not chrom_intervals:
            continue
        cur_start, cur_end = chrom_intervals[0]
        for start, end in chrom_intervals[1:]:
            if start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                merged.append((chrom, cur_start, cur_end))
                cur_start, cur_end = start, end
        merged.append((chrom, cur_start, cur_end))
    return merged


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--width", type=int, default=150)
    parser.add_argument("--fe-cutoff", type=float, default=4.0)
    args = parser.parse_args()

    intervals = {chrom: [] for chrom in CHROM_ORDER}
    for peak_file in args.inputs:
        for chrom, start, end in iter_summit_windows(peak_file, args.width, args.fe_cutoff):
            intervals[chrom].append((start, end))

    merged = merge_intervals(intervals)
    with open(args.output, "w") as out:
        for idx, (chrom, start, end) in enumerate(merged, start=1):
            out.write(f"{chrom}\t{start}\t{end}\tunion_{idx:07d}\t0\t+\n")


if __name__ == "__main__":
    main()
