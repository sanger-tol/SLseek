#!/usr/bin/env python3

import argparse
from Bio import SeqIO, Align
import pandas as pd

parser = argparse.ArgumentParser(description="Retrieve putative SLs - highly expressed kmers.")
parser.add_argument("-p", "--putative", required=True, help="putative SLs sequences - .fasta.")
parser.add_argument("-c", "--clstr", required=True, help="Putative SLs clusters.")
parser.add_argument("-a", "--similarity_score", required=False, help="If true, align putative sequences to c elegans SLs to get a pairwise similarity score", default=False, type=bool)
parser.add_argument("-o", "--output", required=True, help="Output file name basename.")
args = parser.parse_args()

# C. elegans SLs:
sl_list = {
    "SL1": ['AGTCGGTTTAATTACCCAAGTTTGAG', 'CTCAAACTTGGGTAATTAAACCGACT'],
    "SL2": [
        'TGGTTTATACCCAGTTAACCAAG', 'CTTGGTTAACTGGGTATAAACCA',
        'TTGGTTTTAACCCAGTTTAACCAAG', 'CTTGGTTAAACTGGGTTAAAACCAA',
        'GTCGGTTTTAACCCAGTTACTCAAG', 'CTTGAGTAACTGGGTTAAAACCGAC',
        'TTGGTTTTAACCCAGTTACCAAG', 'CTTGGTAACTGGGTTAAAACCAA',
        'TTGGTTTTAACCCATATAACCAAG', 'CTTGGTTATATGGGTTAAAACCAA'
    ]
}

putative_dict = SeqIO.index(args.putative, "fasta")

output_file = f"{args.output}_clustered_putativeSLs.tsv"
with open(args.clstr, "r") as clstr_f, open(output_file, "w") as out_f:
    # Header
    if args.similarity_score:
        out_f.write("SL_ID\tputative_sequence\tentropy\tavg_coverage\tsupporting_kmers\t#_transcripts_in_cluster\tsimilarity_to_SL1\tsimilarity_to_SL2\n")
    else:
        out_f.write("SL_ID\tputative_sequence\tentropy\tavg_coverage\tsupporting_kmers\t#_transcripts_in_cluster\n")

    # Temporary storage
    current_cluster = None
    n_transcripts = 0
    cluster_data = {}

    for line in clstr_f:
        line = line.strip()

        if line.startswith(">Cluster"):
            # New cluster: write the previous one
            if current_cluster and "id_name" in cluster_data:
                if args.similarity_score:
                    out_f.write(f"{cluster_data['id_name']}\t{cluster_data['putative_seq']}\t{cluster_data['entropy']}\t"
                                f"{cluster_data['avg_cov']}\t{cluster_data['supporting_kmers']}\t{n_transcripts}\t"
                                f"{cluster_data['sim_sl1']}\t{cluster_data['sim_sl2']}\n")
                else:
                    out_f.write(f"{cluster_data['id_name']}\t{cluster_data['putative_seq']}\t{cluster_data['entropy']}\t"
                                f"{cluster_data['avg_cov']}\t{cluster_data['supporting_kmers']}\t{n_transcripts}\n")

            # Reset for next cluster
            current_cluster = line
            n_transcripts = 0
            cluster_data = {}
            continue

        # Inside a cluster
        if line:
            n_transcripts += 1
            if "*" in line:  # This is the representative sequence
                parts = line.split()
                seq_info = parts[2].replace("...", "").replace(">", "")
                id_name = seq_info.split("_")[0]
                avg_cov = seq_info.split("avgCov;")[1].split("_")[0]
                entropy = seq_info.split("entropy;")[1].split("_")[0]
                supporting_kmers = seq_info.split("supportingKmers;")[1]
                putative_seq = str(putative_dict[seq_info].seq)

                cluster_data = {
                    "id_name": id_name,
                    "avg_cov": avg_cov,
                    "entropy": entropy,
                    "supporting_kmers": supporting_kmers,
                    "putative_seq": putative_seq
                }

                # Compute similarity if requested
                if args.similarity_score:
                    aligner = Align.PairwiseAligner()
                    cluster_data["sim_sl1"] = max(aligner.score(putative_seq, s) for s in sl_list["SL1"])
                    cluster_data["sim_sl2"] = max(aligner.score(putative_seq, s) for s in sl_list["SL2"])

    # After the loop, write the last cluster
    if current_cluster and "id_name" in cluster_data:
        if args.similarity_score:
            out_f.write(f"{cluster_data['id_name']}\t{cluster_data['putative_seq']}\t{cluster_data['entropy']}\t"
                        f"{cluster_data['avg_cov']}\t{cluster_data['supporting_kmers']}\t{n_transcripts}\t"
                        f"{cluster_data['sim_sl1']}\t{cluster_data['sim_sl2']}\n")
        else:
            out_f.write(f"{cluster_data['id_name']}\t{cluster_data['putative_seq']}\t{cluster_data['entropy']}\t"
                        f"{cluster_data['avg_cov']}\t{cluster_data['supporting_kmers']}\t{n_transcripts}\n")

