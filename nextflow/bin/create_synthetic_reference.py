#!/usr/bin/env python3

import argparse
import random


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--length", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20260311)
    parser.add_argument("--chrom", default="chr1")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    bases = [rng.choice("ACGT") for _ in range(args.length)]

    with open(args.output, "w") as out:
        out.write(f">{args.chrom}\n")
        for i in range(0, len(bases), 80):
            out.write("".join(bases[i:i + 80]) + "\n")


if __name__ == "__main__":
    main()
