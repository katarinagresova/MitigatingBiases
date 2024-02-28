# Mitigating Biases

## Installation

Clone this repository

```bash
git clone git@github.com:katarinagresova/MitigatingBiases.git
```

Switch to the directory

```bash
cd MitigatingBiases
```

Switch to `dev` branch

```bash
git checkout dev
```

Prepare conda environment

```bash
. prepare_conda.sh
```


## Example

Prepare sample data

```bash
python sample_data/generate_fasta.py
```

Run QC tool on sample data

```bash
python src/mitigatingbiases/main.py --pos sample_data/pos.fasta --neg sample_data/neg.fasta --output out.pdf
```

Sample data are randomly generated from the same distribution, they should pass all test:

```bash
Nucleotide composition: PASSED
Dinucleotide composition: PASSED
Nucleotide per position: PASSED
Nucleotides per position reversed: PASSED
Length distribution: PASSED
```