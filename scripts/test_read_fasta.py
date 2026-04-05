from Bio import SeqIO
import os

# Path to your FASTA file
fasta_path = "data/raw_fastas/test.fasta"

# Check if file exists (prevents confusing errors)
if not os.path.exists(fasta_path):
    print(f"File not found: {fasta_path}")
    print("Make sure you placed a FASTA file in data/raw_fastas/")
    exit()

count = 0

print("Reading FASTA file...\n")

for record in SeqIO.parse(fasta_path, "fasta"):
    count += 1
    print(f"Sequence {count}:")
    print(f"ID: {record.id}")
    print(f"Length: {len(record.seq)}\n")

print(f"Total sequences: {count}")