#!/usr/bin/env python3

import argparse
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import seaborn as sns

skiplines = 0

def load_match_data(report_files):
    reports = list()
    samples = None
    for report in tqdm(report_files):
        sample_id = os.path.basename(report).replace('.bamstat.txt', '')
        report = pd.read_csv(report, sep=' ', skiprows=skiplines).rename(columns={'SampleID':'comparison_sample'})
        if samples is None:
            samples = report.iloc[:,0:1]
        report=samples.merge(report, on='comparison_sample', how='left').reset_index()
        report[['perc_het_consistent', 'perc_hom_consistent']] = report[['perc_het_consistent', 'perc_hom_consistent']].fillna(0)
        report = report[['comparison_sample', 'perc_het_consistent', 'perc_hom_consistent']]
        report['bam_sample'] = sample_id
        reports.append(report)
    match_data = pd.concat(reports).sort_values(['bam_sample','comparison_sample'])
    return match_data

def make_heatmaps(match_data):
    het_consistent_array = match_data.pivot(index='bam_sample', columns='comparison_sample', values='perc_het_consistent')
    hom_consistent_array = match_data.pivot(index='bam_sample', columns='comparison_sample', values='perc_hom_consistent')

    #normalize arrays to 0-1 scale
    # hom_consistent_array = (hom_consistent_array-hom_consistent_array.min())/hom_consistent_array.max()
    # het_consistent_array = (het_consistent_array-hom_consistent_array.min())/het_consistent_array.max()

    hom_heatmap_fig, ax1 = plt.subplots(1,1, figsize=(30,30))
    het_heatmap_fig, ax2 = plt.subplots(1,1, figsize=(30,30))
    sns.heatmap(hom_consistent_array, ax=ax1, square=True, xticklabels=1, yticklabels=1)
    sns.heatmap(het_consistent_array, ax=ax2, square=True, xticklabels=1, yticklabels=1)

    ax1.tick_params(axis='x', labelsize=5)
    ax1.tick_params(axis='y', labelsize=5)
    ax2.tick_params(axis='x', labelsize=5)
    ax2.tick_params(axis='y', labelsize=5)

    hom_heatmap_fig.savefig('hom_heatmap.svg')
    het_heatmap_fig.savefig('het_heatmap.svg')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'plot mbv results')

    parser.add_argument("-f", "--input_folder",
                  help="folder with all mbv files to be processed")
    
    parser.add_argument("-o", "--outfile", default='mbv_qtltools_summary_report.tsv',
                  help="outfile for summarized data")

    parser.add_argument("-c", "--conversion_outfile", default='best_dna_match.txt',
                  help="outfile for summarized data")

    # parser.add_argument("-q", "--quiet",
    #               action="store_false", dest="verbose", default=True,
    #               help="don't print status messages to stdout")

    args = parser.parse_args()

    infolder = args.input_folder
    report_files = glob.glob(os.path.join(infolder, '*.bamstat.txt'))

    match_data = load_match_data(report_files).reset_index(drop=True)
    match_data.to_csv(args.outfile, sep='\t', index=False)
    match_data.loc[match_data.groupby('bam_sample').perc_het_consistent.idxmax()][['bam_sample','comparison_sample']].to_csv(args.conversion_outfile, sep='\t', index=False, columns=None)
    make_heatmaps(match_data)
    exit(0)


# loml.exe/run/
# (report)poop.dilithium.grudge.queen/runforestrun/
