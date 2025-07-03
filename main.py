from Bio import SeqIO
from collections import Counter
import pandas as pd
import logomaker
import matplotlib.pyplot as plt

file_name = 'test_input.fasta'
k = 9        # length of motif to search for

# read sequences from file
seqs = [record.seq for record in SeqIO.parse(file_name, 'fasta')]

# count k-mers
kmer_counts = Counter()
for seq in seqs:
    # loop through sequence w/ sliding window of size k
    for i in range(len(seq) - k + 1):
        kmer = str(seq[i : i + k])
        kmer_counts[kmer] += 1

if not kmer_counts:
    print("No k-mers were counted. Check the input file and sequences.")
    exit()
# find most common kmers in kmer_counts
# .most_common(N) returns a list of the N most common (element, count) tuples
top_motifs = kmer_counts.most_common(5)
print(f"\nTop 5 most frequent {k}-mers:")
for motif, count in top_motifs:
    print(str(motif)  + ' | ' + str(count))
    
# most frequent k-mer is used as the motif
top_motif = top_motifs[0][0]

# find instanaces of the top motif in the sequences
instances = []
for seq in seqs:
    for i in range(len(seq) - k + 1):
        kmer = str(seq[i : i + k])
        # finds exact matches
        # use hamming dist. to find fuzzy matches -> the seq logo is only rly useful if you have fuzzy matches (else just 100% for each base)
        if kmer == top_motif:
            instances.append(kmer)

# make count matrix
counts_df = pd.DataFrame(0, index=['A', 'C', 'G', 'T'], columns=range(k))
for ins in instances:
    for position, base in enumerate(ins):
        counts_df.loc[base, position] += 1

print("\nPosition Count Matrix")
print(counts_df)

# Convert counts to frequencies (the PWM)
# Add 1 to each count (pseudocount) to avoid zero probabilities
pwm_df = (counts_df + 1) / (len(instances) + 4) 

print("\nPosition Weight Matrix (PWM)")
print(pwm_df)

# 4. Generate a Sequence Logo
print("\nGenerating sequence logo...")
logo = logomaker.Logo(pwm_df.transpose(), # logomaker wants positions as rows
                      font_name='Arial Rounded MT Bold',
                      color_scheme='classic')
logo.ax.set_ylabel('Bits')
logo.ax.set_title(f'Sequence Logo for {top_motif} motif')
plt.show()