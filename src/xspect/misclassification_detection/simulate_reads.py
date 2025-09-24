"""
Read simulation for the alignment-based misclassification detection (Used for testing purposes).

Notes:
Developed by Oemer Cetin as part of a Bsc thesis at Goethe University Frankfurt am Main (2025).
(An Integration of Alignment-Free and Alignment-Based Approaches for Bacterial Taxon Assignment)
"""

import random
from Bio import SeqIO

__author__ = "Cetin, Oemer"


def extract_random_reads(
    fasta_file, output_fasta, read_length=150, num_reads=1000, seed=42
) -> None:
    """
    Uniformly extracts reads from a genome and writes them to a FASTA-file.

    Args:
        fasta_file (str): Path to input FASTA file.
        output_fasta (str): Output FASTA file to write simulated reads.
        read_length (int): Length of each read to extract.
        num_reads (int): Total number of reads to extract.
        seed (int): A seed for reproducibility.

    Raises:
        ValueError: If the sequences are shorter than the chosen read length.
    """
    random.seed(seed)
    sequences = [
        record
        for record in SeqIO.parse(fasta_file, "fasta")
        if len(record.seq) >= read_length
    ]
    if not sequences:
        raise ValueError("No sequences long enough for the desired read length.")

    # Probability to extract reads from large contigs is higher
    seq_lengths = [len(rec.seq) for rec in sequences]
    total_length = sum(seq_lengths)
    weights = [single_length / total_length for single_length in seq_lengths]

    with open(output_fasta, "w") as o:
        for i in range(num_reads):
            # random.choices() provides a list!
            selected = random.choices(sequences, weights=weights, k=1)[0]
            seq_length = len(selected.seq)
            start = random.randint(0, seq_length - read_length)
            read_seq = selected.seq[start : start + read_length]
            o.write(
                f">read_{i}_{selected.id}_{start}-{start + read_length}\n{read_seq}\n"
            )
    print("The reads have been simulated successfully.")
