from collections import Counter
import pandas as pd

def count_nucleotides(sequence) -> dict:
    """
    A function to count the nucleotides in a sequence
    :param sequence: the sequence
    :return: the dictionary of nucleotide counts
    """
    sum_sequence = 0
    counter = Counter()
    for nucleotide in sequence:
        counter[nucleotide] += 1
        sum_sequence += 1  # Update the total count for each nucleotide
    counter = dict(counter)
    counter['sum'] = len(sequence)
    return counter

def count_dinucleotides(sequence) -> dict:
    """
    A function to count the dinucleotides in a sequence
    :param sequence: the sequence
    :return: the dictionary of dinucleotide counts
    """
    sum_sequence = 0
    counter = Counter()
    for i in range(len(sequence) - 1):
        dinucleotide = str(sequence[i:i + 2])
        counter[dinucleotide] += 1
        sum_sequence += 1  # Update the total count for each nucleotide
    counter = dict(counter)
    counter['sum'] = len(sequence)
    return counter

def get_nucleotides_per_position(sequence, nucleotides_per_position, reverse = False):
    """
    A function to count the nucleotides in a sequence per each position
    :param sequences: sequences
    :param reverse: reverse the sequences
    :return: the dictionary of nucleotide counts per position
    """

    if reverse:
        sequence = sequence[::-1]
    for i, nucleotide in enumerate(sequence):
        if i not in nucleotides_per_position:
            nucleotides_per_position[i] = Counter()
        nucleotides_per_position[i][nucleotide] += 1

    return nucleotides_per_position

def do_counting(sequences):

    nucleotides_per_position = {}
    nucleotides_per_position_reversed = {}
    nucleotides = []
    dinucleotides = []

    for sequence in sequences:
        nucleotides.append(count_nucleotides(sequence))
        dinucleotides.append(count_dinucleotides(sequence))
        nucleotides_per_position = get_nucleotides_per_position(sequence, nucleotides_per_position)
        nucleotides_per_position_reversed = get_nucleotides_per_position(sequence, nucleotides_per_position_reversed, reverse = True)

    nucleotides_df = pd.DataFrame(nucleotides).fillna(0)
    dinucleotides_df = pd.DataFrame(dinucleotides).fillna(0)

    nucleotides_per_position_df = pd.DataFrame(nucleotides_per_position).T.fillna(0)
    nucleotides_per_position_df['sum'] = nucleotides_per_position_df.sum(axis=1)

    nucleotides_per_position_reversed_df = pd.DataFrame(nucleotides_per_position_reversed).T.fillna(0)
    nucleotides_per_position_reversed_df['sum'] = nucleotides_per_position_reversed_df.sum(axis=1)

    return nucleotides_df, dinucleotides_df, nucleotides_per_position_df, nucleotides_per_position_reversed_df

        