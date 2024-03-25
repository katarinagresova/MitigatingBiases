import random
import string

def generate_bad_data(file_path, num_sequences):
    """
    Function to generate bad data. It is fundamentally the same as the base one, but it is renamed for clarity purposes.
    """
    with open(file_path, 'w') as fasta_file:
        for i in range(num_sequences):
            # generate random length for each sequence
            length = random.randint(20, 100)
            sequence = generate_bad_sequence(length)
            header = f">Sequence {i+1}"
            fasta_file.write(f"{header}\n{sequence}\n")

def generate_bad_sequence(length=10):
    """
    Function to generate false data/sequences. It should randomly generate between non-existent bases,
    special characters, line breaks, and even empty sequences all at once. 
    In order to generate false data, you can freely run `py sample_data/generate_bad_fasta.py` in your terminal.
    WARNING: There is a chance when running the analysis with this data that Python crashes.
    """
    nucleotides = ['A', 'C', 'G', 'T']
    bad_chars = list(set(string.printable) - set(nucleotides))
    sequence = ''.join(random.choice(bad_chars) for _ in range(length))
    return sequence

# Usage example
num_sequences = 5000

generate_bad_data('sample_data/pos_bad.fasta', num_sequences)
generate_bad_data('sample_data/neg_bad.fasta', num_sequences)

