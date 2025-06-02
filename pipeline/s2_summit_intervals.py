import pandas as pd

# Load the file (no header)
df = pd.read_csv("all_peaks.bed", sep="\t", header=None)

# Assign column names
df.columns = [
    "chrom", "start", "end", "name", "score", "strand",
    "fold_enrichment", "pvalue", "qvalue", "summit_offset"
]

# Filter by fold enrichment > 4
df = df[df["fold_enrichment"] > 4]

# Compute summit position
df["summit"] = df["start"] + df["summit_offset"]

# Sort by chromosome and summit
df = df.sort_values(by=["chrom", "summit"])

# Compute intervals between summits per chromosome
intervals = []
for chrom, group in df.groupby("chrom"):
    summits = group["summit"].values
    for i in range(len(summits) - 1):
        interval = int(summits[i+1] - summits[i])
        intervals.append(interval)

# Save output: one interval per line, no header
with open("summit_intervals.txt", "w") as f:
    for interval in intervals:
        f.write(f"{interval}\n")