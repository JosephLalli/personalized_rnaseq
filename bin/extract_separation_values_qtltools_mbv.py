#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd


def parse_mbv_report(file, output_file):
    report = pd.read_csv(file, sep=' ')
    sample_name = file.split('/')[-1].split('.')[0]
    hetero_match_idx = report.perc_het_consistent.argmax()
    homo_match_idx = report.perc_hom_consistent.argmax()
    sorted_hetero_dists = np.sort(report.perc_het_consistent)
    sorted_homo_dists = np.sort(report.perc_hom_consistent)
    top_hetero_separation = (sorted_hetero_dists[-1] - np.mean(sorted_hetero_dists[:-1])) / (sorted_hetero_dists[-2] - np.mean(sorted_hetero_dists[:-2]))
    top_homo_separation = (sorted_homo_dists[-1] - np.mean(sorted_homo_dists[:-1])) / (sorted_homo_dists[-2] - np.mean(sorted_homo_dists[:-2]))

    if hetero_match_idx == homo_match_idx:
        match_idx = hetero_match_idx
        match_name = report.SampleID.values[match_idx]
    else:
        match_name = 'NA'
    output = pd.DataFrame(({'Sample Name': sample_name,
                            'Matched Sample': match_name,
                            "Relative Distance of Heterozygous Matches": top_hetero_separation,
                            "Relative Distance of Homozygous Matches": top_homo_separation},))
    output.to_csv(output_file, sep='\t', index=False)
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Extract the closest match from an mbv report, and the distance to the second-closest match"""
    )
    parser.add_argument("--report", type=str, help="mbv report")
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="top_match.tsv",
        type=str,
    )

    args = parser.parse_args()
    parse_mbv_report(args.report, args.output)