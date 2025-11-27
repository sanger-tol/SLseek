#!/usr/bin/env python3
import pandas as pd
from plotnine import *
from mizani.formatters import comma_format, scientific_format, percent_format
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--tsv", required=True, help="Input TSV")
# parser.add_argument("--sep", required=True, help="Separaton '\t' for fastK and ' ' for jellyfish")
parser.add_argument("--output", required=True, help="Output name without extension")
args = parser.parse_args()

# Read data
data = pd.read_csv(
    args.tsv,
    sep=None,
    engine='python',
    header=None,
    names=["abundance", "X.of.kmers"]
)

p1 = (
    ggplot(data, aes(x="X.of.kmers", y="abundance"))
    + geom_col(fill="darkgreen")
    + labs(
        title="Kmers distribution (Frequency vs. Coverage)",
        x="# of Kmers (Frequency)(Log)",
        y="Kmers abundance (Coverage)"
    )
    + theme_minimal()
    + scale_x_continuous(
        trans='log10',
        labels=comma_format()
    )
)



# Save image
out_file = args.output + "_histogram.png"
p1.save(out_file, dpi=300)

