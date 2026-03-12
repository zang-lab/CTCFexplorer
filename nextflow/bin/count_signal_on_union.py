#!/usr/bin/env python3

import argparse

import pandas as pd
import pysam


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--union", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    total_mapped = bam.mapped
    rows = []

    with open(args.union) as handle:
        for line in handle:
            chrom, start, end, union_id, score, strand = line.rstrip("\n").split("\t")[:6]
            start = int(start)
            end = int(end)
            region_len = max(1, end - start)
            count = bam.count(contig=chrom, start=start, end=end, read_callback="all")
            rpkm = 0.0
            if total_mapped > 0:
                rpkm = (count * 1_000_000_000.0) / (total_mapped * region_len)
            rows.append((union_id, round(rpkm, 4)))

    bam.close()

    df = pd.DataFrame(rows, columns=["Union_ID", args.sample]).set_index("Union_ID")
    df.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()
