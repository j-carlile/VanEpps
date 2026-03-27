from Bio import SeqIO

fasta_path = "data/raw_fastas/test.fasta"

count = 0
for record in SeqIO.parse(fasta_path, "fasta"):
    count += 1
    print(record.id)

print(f"Total sequences: {count}")