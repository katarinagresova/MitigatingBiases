from counting import do_counting
from plotting.plotting import plot_per_base_sequence_comparison, plot_composition_comparison, plot_composition_comparison_boxplot, plot_lenght_comparison

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Bio import SeqIO
import argparse

def parse_fasta(filename):
    """
    Parses a fasta file and returns a generator of SeqIO objects
    """

    with open(filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            yield record.seq

def print_stats(stats):
    # Print stats. Failed in red and passed in green
    for key, value in stats.items():
        if value:
            print(f'{key}: \033[92mPASSED\033[0m')
        else:
            print(f'{key}: \033[91mFAILED\033[0m')

def save_plots(plots, output_file):
    with PdfPages(output_file) as pdf:
        for fig in plots:
            pdf.savefig(fig, bbox_inches='tight')

def parse_args():
    parser = argparse.ArgumentParser(description='Mitigating biases')
    parser.add_argument('--pos', type=str, required=True, help='Path to positive sequences')
    parser.add_argument('--neg', type=str, required=True, help='Path to negative sequences')
    parser.add_argument('--output', type=str, required=True, help='Path to output file')
    return parser.parse_args()

def get_plots_and_stats(
    pos_nucleotides_df,
    neg_nucleotides_df,
    pos_dinucleotides_df,
    neg_dinucleotides_df,
    pos_nucleotides_per_position_df,
    neg_nucleotides_per_position_df,
    pos_nucleotides_per_position_reversed_df,
    neg_nucleotides_per_position_reversed_df,
    p_value_thresh=0.001
):
    plots = []
    stats = {}

    # Plot nucleotides
    fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(10, 8))
    _, passed = plot_composition_comparison(
        pos_nucleotides_df,
        neg_nucleotides_df,
        x_label='Frequency',
        ax=ax,
        stats=True,
        p_value_thresh=p_value_thresh)
    fig.suptitle('Nucleotide composition')
    plots.append(fig)
    stats['Nucleotide composition'] = passed

    # Plot dinucleotides
    fig, ax = plt.subplots(nrows=16, ncols=1, figsize=(10, 32))
    _, passed = plot_composition_comparison(
        pos_dinucleotides_df,
        neg_dinucleotides_df,
        x_label='Frequency',
        ax=ax,
        stats=True,
        p_value_thresh=p_value_thresh)
    fig.suptitle('Dinucleotide composition')
    plots.append(fig)
    stats['Dinucleotide composition'] = passed

    # Plot nucleotides per position
    fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(10, 8))
    _, passed = plot_per_base_sequence_comparison(
        pos_nucleotides_per_position_df, 
        neg_nucleotides_per_position_df, 
        end_position=100,  
        x_label='Position in read (bp)',
        ax=ax,
        stats=True,
        p_value_thresh=p_value_thresh)
    fig.suptitle('Nucleotide per position')
    plots.append(fig)
    stats['Nucleotide per position'] = passed

    # Plot nucleotides per position reversed
    fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(10, 8))
    _, passed = plot_per_base_sequence_comparison(
        pos_nucleotides_per_position_reversed_df, 
        neg_nucleotides_per_position_reversed_df, 
        end_position=100, 
        x_label='Position in read (bp)',
        ax=ax,
        stats=True,
        p_value_thresh=p_value_thresh)
    fig.suptitle('Nucleotide per position reversed')
    plots.append(fig)
    stats['Nucleotides per position reversed'] = passed

    # Plot length distribution
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 3))
    _, passed = plot_lenght_comparison(
        pos_nucleotides_df, 
        neg_nucleotides_df, 
        x_label='Sequence length', 
        ax=ax,
        stats=True,
        p_value_thresh=p_value_thresh)
    fig.suptitle('Length distribution')
    plots.append(fig)
    stats['Length distribution'] = passed

    return plots, stats

def main():

    args = parse_args()

    pos_seq = parse_fasta(args.pos)
    neg_seq = parse_fasta(args.neg)

    pos_nucleotides_df, pos_dinucleotides_df, pos_nucleotides_per_position_df, pos_nucleotides_per_position_reversed_df = do_counting(pos_seq)
    neg_nucleotides_df, neg_dinucleotides_df, neg_nucleotides_per_position_df, neg_nucleotides_per_position_reversed_df = do_counting(neg_seq)

    plots, stats = get_plots_and_stats(
        pos_nucleotides_df,
        neg_nucleotides_df,
        pos_dinucleotides_df,
        neg_dinucleotides_df,
        pos_nucleotides_per_position_df,
        neg_nucleotides_per_position_df,
        pos_nucleotides_per_position_reversed_df,
        neg_nucleotides_per_position_reversed_df
    )

    print_stats(stats)
    save_plots(plots, output_file=args.output)

if __name__ == "__main__":
    main()