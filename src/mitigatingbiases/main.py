from counting import do_counting
from plotting.plotting import plot_per_base_sequence_content, plot_per_base_sequence_comparison, plot_composition_comparison, plot_composition_comparison_boxplot, plot_lenght_comparison

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Bio import SeqIO

def parse_fasta(filename):
    """
    Parses a fasta file and returns a generator of SeqIO objects
    """

    with open(filename, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            yield record.seq


def main():
    pos_seq = parse_fasta("/home/jovyan/MitigatingBiases/data/pos.fasta")
    neg_seq = parse_fasta("/home/jovyan/MitigatingBiases/data/neg.fasta")

    pos_nucleotides_df, pos_dinucleotides_df, pos_nucleotides_per_position_df, pos_nucleotides_per_position_reversed_df = do_counting(pos_seq)
    neg_nucleotides_df, neg_dinucleotides_df, neg_nucleotides_per_position_df, neg_nucleotides_per_position_reversed_df = do_counting(neg_seq)

    plots = []

    # Plot nucleotides
    fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(10, 8))
    plot_composition_comparison(pos_nucleotides_df, neg_nucleotides_df, x_label='Position in read (bp)', ax=ax)
    fig.suptitle('Nucleotide composition')
    plots.append(fig)

    # Plot dinucleotides
    fig, ax = plt.subplots(nrows=16, ncols=1, figsize=(10, 32))
    plot_composition_comparison(pos_dinucleotides_df, neg_dinucleotides_df, x_label='Position in read (bp)', ax=ax)
    fig.suptitle('Dinucleotide composition')
    plots.append(fig)

    # Plot nucleotides per position
    fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(10, 8))
    plot_per_base_sequence_comparison(
        pos_nucleotides_per_position_df, 
        neg_nucleotides_per_position_df, 
        end_position=100,  
        x_label='Position in read (bp)',
        ax=ax,
        stats=True)
    fig.suptitle('Nucleotide per position')
    plots.append(fig)

    # Plot nucleotides per position reversed
    fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(10, 8))
    plot_per_base_sequence_comparison(
        pos_nucleotides_per_position_reversed_df, 
        neg_nucleotides_per_position_reversed_df, 
        end_position=100, 
        x_label='Position in read (bp)',
        ax=ax,
        stats=True)
    fig.suptitle('Nucleotide per position reversed')
    plots.append(fig)

    # Plot length distribution
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 3))
    plot_lenght_comparison(
        pos_nucleotides_df, 
        neg_nucleotides_df, 
        x_label='Sequence length', 
        ax=ax)
    fig.suptitle('Length distribution')
    plots.append(fig)

    with PdfPages('output.pdf') as pdf:
        for fig in plots:
            pdf.savefig(fig, bbox_inches='tight') 

if __name__ == "__main__":
    main()