#!/usr/bin/env python3
import pandas as pd
from plotnine import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--tsv", required=True, help="Input TSV")
parser.add_argument("--output", required=True, help="Output name without extension")
args = parser.parse_args()

# Read data
data = pd.read_csv(args.tsv, sep="\t")

# Columns to normalize
numeric_cols = ["entropy", "avg_coverage", "supporting_kmers", "#_transcripts_in_cluster"]

# numeric_cols <- c("entropy", "avg_coverage", "supporting_kmers", "`#_transcripts_in_cluster`")

# Normalization function
def normalize(x):
    return (x - x.mean())

# Create normalized columns
for col in numeric_cols:
    data[f"normalized_{col}"] = normalize(data[col])

# Reshape to long format
heatmap_data = (
    data.melt(
        id_vars="SL_ID",
        value_vars=[f"normalized_{c}" for c in numeric_cols],
        var_name="Metric",
        value_name="Normalized_Value"
    )
)

# Clean metric names
heatmap_data["Metric"] = heatmap_data["Metric"].str.replace("normalized_", "")

# Order metrics like in R
metric_order = ["entropy", "avg_coverage", "supporting_kmers", "#_transcripts_in_cluster"]
heatmap_data["Metric"] = pd.Categorical(heatmap_data["Metric"], categories=metric_order, ordered=True)

# Build heatmap plot
p = (
    ggplot(heatmap_data, aes("SL_ID", "Metric", fill="Normalized_Value"))
    + geom_tile()
    + scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, name="Normalized value")
    + theme_minimal()
    + theme(
        axis_text_x=element_text(angle=90, hjust=1, vjust=0.5, size=8)
    )
)

# Save image
out_file = args.output + "_heatmap.png"
p.save(out_file, dpi=300)
print(f"Saved heatmap to {out_file}")

