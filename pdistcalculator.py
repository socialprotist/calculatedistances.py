import sys
from Bio import SeqIO
import numpy as np

def calculate_uncorrected_distance(seq1, seq2):
    """Calculate uncorrected pairwise distance between two sequences, ignoring gaps."""
    total_sites = 0
    differing_sites = 0
    
    for a, b in zip(seq1, seq2):
        # Ignore sites with gaps
        if a != '-' and b != '-':
            total_sites += 1
            if a != b:
                differing_sites += 1
                
    if total_sites == 0:
        return 0  # Avoid division by zero if no comparable sites exist
    
    return differing_sites / total_sites

def generate_pairwise_distance_matrix(sequences, sequence_names):
    """Generate an uncorrected pairwise distance matrix."""
    num_sequences = len(sequences)
    distance_matrix = np.zeros((num_sequences, num_sequences), dtype=float)

    for i in range(num_sequences):
        for j in range(i, num_sequences):
            distance = calculate_uncorrected_distance(sequences[i], sequences[j])
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance

    return distance_matrix

def read_fasta_alignment(file_path):
    """Read sequences and names from a FASTA alignment."""
    sequences = []
    sequence_names = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
        sequence_names.append(record.id)
    return sequences, sequence_names

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <fasta_file>")
        sys.exit(1)

    fasta_file_path = sys.argv[1]
    
    # Read sequences and names from the FASTA alignment
    sequences, sequence_names = read_fasta_alignment(fasta_file_path)
    
    # Generate the uncorrected pairwise distance matrix
    distance_matrix = generate_pairwise_distance_matrix(sequences, sequence_names)
    
    # Print the lower half of the distance matrix with sequence names and differences
    print("Uncorrected Pairwise Distance Matrix (lower half, ignoring gaps):")
    print("\t" + "\t".join(sequence_names))
    
    for i in range(len(sequences)):
        print(sequence_names[i], end="\t")
        for j in range(len(sequences)):
            if j <= i:  # Only print the lower half, including the diagonal
                differing_sites = sum(1 for a, b in zip(sequences[i], sequences[j]) if a != b and a != '-' and b != '-')
                print(f"{distance_matrix[i][j]:.4f} ({differing_sites})", end="\t")
            else:
                print("", end="\t")  # Leave empty space for the upper half
        print()

