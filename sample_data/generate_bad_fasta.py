import random
import string

def generate_bad_data_file_single(file_path, num_sequences):
    """Only one nucleotide"""
    with open(file_path, 'w') as fasta_file:
        for i in range(num_sequences):
            length = random.randint(20, 100)
            sequence = generate_single_base_sequences(length)
            header = f">Sequence {i+1}"
            fasta_file.write(f"{header}\n{sequence}\n")

def generate_empty_data_file(file_path, num_sequences):
    """Empty file"""
    with open(file_path, 'w') as fasta_file:
        fasta_file.write("")

def generate_bad_data_file_emptyseq(file_path, num_sequences):
    """Empty sequences"""
    with open(file_path, 'w') as fasta_file:
        for i in range(num_sequences):
            length = random.randint(20, 100)
            sequence = generate_empty_sequences(length)
            header = f">Sequence {i+1}"
            fasta_file.write(f"{header}\n{sequence}\n")

def generate_bad_data_file_mixedseq(file_path, num_sequences):
    """Mixed sequence with special chars"""
    with open(file_path, 'w') as fasta_file:
        for i in range(num_sequences):
            length = random.randint(20, 100)
            sequence = generate_mixed_sequences(length)
            header = f">Sequence {i+1}"
            fasta_file.write(f"{header}\n{sequence}\n")

def generate_bad_data_file_specialseq(file_path, num_sequences):
    """Sequence with only special chars"""
    with open(file_path, 'w') as fasta_file:
        for i in range(num_sequences):
            length = random.randint(20, 100)
            sequence = generate_special_characters(length)
            header = f">Sequence {i+1}"
            fasta_file.write(f"{header}\n{sequence}\n")

def generate_bad_data_file_nonexist(file_path, num_sequences):
    """Sequence with non-existent bases"""
    with open(file_path, 'w') as fasta_file:
        for i in range(num_sequences):
            length = random.randint(20, 100)
            sequence = generate_nonexistent_bases(length)
            header = f">Sequence {i+1}"
            fasta_file.write(f"{header}\n{sequence}\n")

def generate_bad_data_file_casemix(file_path, num_sequences):
    """Sequence with upper- and lowercase bases"""
    with open(file_path, 'w') as fasta_file:
        for i in range(num_sequences):
            length = random.randint(20, 100)
            sequence = generate_upper_lower_mix(length)
            header = f">Sequence {i+1}"
            fasta_file.write(f"{header}\n{sequence}\n")

def generate_bad_data_file_whitespace(file_path, num_sequences):
    """Sequence with whitespace"""
    with open(file_path, 'w') as fasta_file:
        for i in range(num_sequences):
            length = random.randint(20, 100)
            sequence = generate_whitespace_sequences(length)
            header = f">Sequence {i+1}"
            fasta_file.write(f"{header}\n{sequence}\n")

def generate_bad_data_file_linebreaks(file_path, num_sequences):
    """Sequence with line breaks"""
    with open(file_path, 'w') as fasta_file:
        for i in range(num_sequences):
            length = random.randint(20, 100)
            sequence = generate_line_breaks_delimiters(length)
            header = f">Sequence {i+1}"
            fasta_file.write(f"{header}\n{sequence}\n")

def generate_bad_data_file_double(file_path, num_sequences):
    """Sequence with duplicate bases"""
    with open(file_path, 'w') as fasta_file:
        for i in range(num_sequences):
            length = random.randint(20, 100)
            sequence = generate_duplicate_bases(length)
            header = f">Sequence {i+1}"
            fasta_file.write(f"{header}\n{sequence}\n")

#____________________________________________CREATING THE SEQUENCES____________________________________________________
def generate_single_base_sequences(length = 10):
    nucleotide = ['A']
    single_sequence = ''.join(map(str, nucleotide))
    return single_sequence

def generate_empty_sequences(length = 10):
    nucleotides = ['']
    empty_sequence = ''.join(random.choice(nucleotides) for _ in range(length))
    return empty_sequence

def generate_mixed_sequences(length = 10):
    bases = ['A', 'C', 'G', 'T']
    special_chars = ['*', '\\', ';', '~']
    mixed_sequence = ''.join(random.choice(bases + special_chars) for _ in range(length))
    return mixed_sequence

def generate_special_characters(length=10):
    special_bases = ['*','~','-','@','/','?']
    special_sequence = ''.join(random.choice(special_bases) for _ in range(length))
    return special_sequence

def generate_nonexistent_bases(length=10):
    false_bases = ['X','Y','Z','Q','R','H']
    false_sequence = ''.join(random.choice(false_bases) for _ in range(length))
    return false_sequence
    
def generate_upper_lower_mix(length=10):
    mixed_bases = ['A','C','G','T','a','c','g','t']
    mixed_sequence = ''.join(random.choice(mixed_bases) for _ in range(length))
    return mixed_sequence

def generate_whitespace_sequences(length=10):
    whitespace_nucleotides = ['A','C','G','T',' ']
    whitespace_sequence = ''.join(random.choice(whitespace_nucleotides) for _ in range(length))
    return whitespace_sequence

def generate_line_breaks_delimiters(length=10):
    linebreaks_nucleotides = ['A','C','G','T','\r\n','\t']
    linebreaks_sequence = ''.join(random.choice(linebreaks_nucleotides) for _ in range(length))
    return linebreaks_sequence
    
def generate_duplicate_bases(length=10):
    double_bases = ['A','CC','G','TT']
    doublebase_seq = ''.join(random.choice(double_bases) for _ in range(length))
    return doublebase_seq

# usage example
num_sequences = 100

generate_bad_data_file_single('pos_single.fasta', num_sequences)
generate_bad_data_file_single('neg_single.fasta', num_sequences)
generate_empty_data_file('pos_empty.fasta', num_sequences)
generate_empty_data_file('neg_empty.fasta', num_sequences)
generate_bad_data_file_emptyseq('pos_emptyseq.fasta', num_sequences)
generate_bad_data_file_emptyseq('neg_emptyseq.fasta', num_sequences)
generate_bad_data_file_mixedseq('pos_mixedseq.fasta', num_sequences)
generate_bad_data_file_mixedseq('neg_mixedseq.fasta', num_sequences)
generate_bad_data_file_specialseq('pos_specialseq.fasta', num_sequences)
generate_bad_data_file_specialseq('neg_specialseq.fasta', num_sequences)
generate_bad_data_file_nonexist('pos_nonexist.fasta', num_sequences)
generate_bad_data_file_nonexist('neg_nonexist.fasta', num_sequences)
generate_bad_data_file_casemix('pos_casemix.fasta', num_sequences)
generate_bad_data_file_casemix('neg_casemix.fasta', num_sequences)
generate_bad_data_file_whitespace('pos_ws.fasta', num_sequences)
generate_bad_data_file_whitespace('neg_ws.fasta', num_sequences)
generate_bad_data_file_double('pos_double.fasta', num_sequences)
generate_bad_data_file_double('neg_double.fasta', num_sequences)
generate_bad_data_file_linebreaks('pos_lb.fasta', num_sequences)
generate_bad_data_file_linebreaks('neg_lb.fasta', num_sequences)