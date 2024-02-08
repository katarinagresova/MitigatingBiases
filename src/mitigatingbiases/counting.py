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
    counter['sum'] = sum_sequence
    return counter

def get_nucleotides_per_position(sequences, reverse = False):
    """
    A function to count the nucleotides in a sequence per each position
    :param sequences: sequences
    :param reverse: reverse the sequences
    :return: the dictionary of nucleotide counts per position
    """

    nucleotides_per_position = {}
    for sequence in sequences:
        if reverse:
            sequence = sequence[::-1]
        for i, nucleotide in enumerate(sequence):
            if i not in nucleotides_per_position:
                nucleotides_per_position[i] = Counter()
            nucleotides_per_position[i][nucleotide] += 1
    
    # transform dict of counters to pandas dataframe
    nucleotides_per_position = pd.DataFrame(nucleotides_per_position).T.fillna(0)

    return nucleotides_per_position

        