from Bio import SeqIO
from collections import Counter

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

# find most common kmers in kmer_counts
if kmer_counts:
    # .most_common(N) returns a list of the N most common (element, count) tuples
    top_motifs = kmer_counts.most_common(5)
    
    print(f"\nTop 5 most frequent {k}-mers:")
    for motif, count in top_motifs:
        print(str(motif)  + ' | ' + str(count))
else:
    print("No k-mers were counted. Check the input file and sequences.")

top_motif = top_motifs[0][0] # the most frequent k-mer

