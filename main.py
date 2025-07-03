import random
from Bio import SeqIO

# generate random sequences with a specified motif (known data for testing)
def generate_dna(length=200):
    return ''.join(random.choice('ATCG') for _ in range(length))

# generate a list of sequences to test with
motif = 'AGGTCCCCA'
seqs = []
for _ in range(10):
    random_seq = generate_dna()
    insert_motif = random.randint(0, len(random_seq) - len(motif))
    seq = random_seq[:insert_motif] + motif + random_seq[insert_motif:]
    seqs.append(seq)

# write to a fasta file
output_filename = 'test_input.fasta'
with open(output_filename, 'w') as file:
    for i in range(len(seqs)):
        file.write('>N' + str(i) + '\n')
        file.write(seqs[i] + '\n')