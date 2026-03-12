#!/usr/bin/env python3

import argparse
import gzip
import hashlib
import random
from pathlib import Path


def read_fasta(path):
    seq_parts = []
    chrom = None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                chrom = line[1:].split()[0]
            else:
                seq_parts.append(line.upper())
    if chrom is None:
        raise ValueError(f"No FASTA header found in {path}")
    return chrom, "".join(seq_parts)


def sample_profile(sample):
    profiles = [
        [1500, 5000, 8500],
        [1500, 5000, 12000],
        [1500, 8500, 12000],
    ]
    digest = hashlib.sha256(sample.encode()).digest()[0]
    return profiles[digest % len(profiles)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--read-count", type=int, required=True)
    parser.add_argument("--read-length", type=int, default=75)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    chrom, seq = read_fasta(args.reference)
    seq_len = len(seq)
    peaks = sample_profile(args.sample)
    weights = [0.5, 0.3, 0.2]
    counts = [int(args.read_count * w) for w in weights]
    counts[-1] = args.read_count - sum(counts[:-1])

    seed = int(hashlib.sha256(args.sample.encode()).hexdigest()[:16], 16)
    rng = random.Random(seed)
    quality = "I" * args.read_length

    with gzip.open(args.output, "wt") as out:
        read_index = 1
        for peak_center, peak_count in zip(peaks, counts):
            start_min = max(0, peak_center - 350)
            start_max = min(seq_len - args.read_length, peak_center + 350)
            if start_min >= start_max:
                raise ValueError(f"Peak center {peak_center} is too close to edge for read length {args.read_length}")

            for _ in range(peak_count):
                read_start = rng.randint(start_min, start_max)
                read_seq = seq[read_start:read_start + args.read_length]
                out.write(f"@{args.sample}.{read_index} {chrom}:{read_start + 1}-{read_start + args.read_length}\n")
                out.write(read_seq + "\n+\n" + quality + "\n")
                read_index += 1

    print(f"sample={args.sample}")
    print(f"reads={args.read_count}")
    print(f"reference={Path(args.reference).resolve()}")
    print(f"output={Path(args.output).resolve()}")


if __name__ == "__main__":
    main()
