#!/usr/bin/env python3

import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    frames = [pd.read_csv(path, sep="\t", index_col=0) for path in args.inputs]
    merged = pd.concat(frames, axis=1).fillna(0.0)
    merged = merged.reindex(sorted(merged.columns), axis=1)
    merged.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()
