#!/usr/bin/env python3

import argparse
import math

import numpy as np
import pandas as pd


def quantile_normalize(df):
    arr = df.to_numpy(dtype=float)
    sort_idx = np.argsort(arr, axis=0)
    sorted_arr = np.sort(arr, axis=0)
    rank_mean = sorted_arr.mean(axis=1)
    norm = np.zeros_like(arr)
    for col in range(arr.shape[1]):
        norm[sort_idx[:, col], col] = rank_mean
    return pd.DataFrame(norm, index=df.index, columns=df.columns)


def write_bed(df, path):
    df.loc[:, ["chr", "start", "end", "Union_ID"]].to_csv(path, sep="\t", index=False, header=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--union", required=True)
    parser.add_argument("--occupancy-matrix", required=True)
    parser.add_argument("--occupancy-scores", required=True)
    parser.add_argument("--signal-matrix", required=True)
    parser.add_argument("--sample-count", required=True, type=int)
    parser.add_argument("--high-confidence-min-samples", required=True, type=int)
    parser.add_argument("--constitutive-min-fraction", required=True, type=float)
    args = parser.parse_args()

    occupancy_matrix = pd.read_csv(args.occupancy_matrix, sep="\t", index_col=0)
    occupancy_scores = pd.read_csv(args.occupancy_scores, sep="\t")
    signal_matrix = pd.read_csv(args.signal_matrix, sep="\t", index_col=0)

    constitutive_min_samples = int(math.ceil(args.sample_count * args.constitutive_min_fraction))

    high_confidence_ids = occupancy_scores.loc[
        occupancy_scores["occupancy_score"] >= args.high_confidence_min_samples, "Union_ID"
    ]
    constitutive_ids = occupancy_scores.loc[
        occupancy_scores["occupancy_score"] >= constitutive_min_samples, "Union_ID"
    ]

    high_confidence_df = occupancy_scores[occupancy_scores["Union_ID"].isin(high_confidence_ids)].copy()
    constitutive_df = occupancy_scores[occupancy_scores["Union_ID"].isin(constitutive_ids)].copy()

    write_bed(high_confidence_df, "high_confidence.bed")
    write_bed(constitutive_df, "constitutive.bed")

    occupancy_matrix.loc[high_confidence_ids].to_csv("union_occupancy_high_confidence.tsv", sep="\t")
    occupancy_matrix.loc[constitutive_ids].to_csv("union_occupancy_constitutive.tsv", sep="\t")

    signal_matrix.loc[high_confidence_ids].fillna(0.0).to_csv("union_RPKM_high_confidence.tsv", sep="\t")
    signal_matrix.loc[constitutive_ids].fillna(0.0).to_csv("union_RPKM_constitutive.tsv", sep="\t")

    qn = quantile_normalize(signal_matrix.loc[high_confidence_ids].fillna(0.0)).round(4)
    qn.to_csv("union_RPKM_high_confidence_QN.tsv", sep="\t")


if __name__ == "__main__":
    main()
