import random

def generate_fasta_file(file_path, num_sequences):
    with open(file_path, 'w') as fasta_file:
        for i in range(num_sequences):
            # generate random length for each sequence
            length = random.randint(10, 100)
            sequence = generate_random_sequence(length)
            header = f">Sequence {i+1}"
            fasta_file.write(f"{header}\n{sequence}\n")

def generate_random_sequence(length=10):
    nucleotides = ['A', 'C', 'G', 'T']
    sequence = ''.join(random.choice(nucleotides) for _ in range(length))
    return sequence

# Usage example
num_sequences = 50

generate_fasta_file('data/pos.fasta', num_sequences)
generate_fasta_file('data/neg.fasta', num_sequences)
