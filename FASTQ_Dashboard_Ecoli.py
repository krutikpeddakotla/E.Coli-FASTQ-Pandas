## Pandas is not optimal for FASTQ parsing due to high memory usage. Tools like FastQC are optimal.
# This is just a learning project to test Pandas knowledge.

import pandas as pd
import re

def read_fastq(fastq_file, max_reads = 1000):
    reads = []
    with open(fastq_file) as f:
        count = 0
        while True:
            header = f.readline().rstrip()
            if not header:
                break
            seq = f.readline().rstrip()
            plus = f.readline().rstrip()
            qual = f.readline().rstrip()
            reads.append((header, seq, qual))
            count += 1
            if count >= max_reads:
                break
    return pd.DataFrame(reads, columns=['header', 'seq', 'qual'])

df = read_fastq('SRR36764531.fastq', max_reads = 1000)

# per read metrics
df["length"] = df["seq"].str.len()
df["gc_content"] = (df["seq"].str.count("G") + df["seq"].str.count("C")) / df["length"]
df["n_count"] = df["seq"].str.count("N") / df["length"]
df["mean_q"] = df["qual"].apply(lambda a: sum(ord(c) - 33 for c in a)/len(a))

## Count homopolymers
def count_homopolymers(seq, min_length = 3):
    """Count the number of homopolymers in a sequence using regex matching"""
    if pd.isna(seq) or not seq:
        return 0
    # Pattern matching
    pattern = r'(.)\1{' + str(min_length - 1) + r',}'
    matches = re.findall(pattern, seq)
    return len(matches)
df["homopolymer_count"] = df["seq"].apply(lambda x: count_homopolymers(x, min_length = 3))

## k-mer complexity calculation
def kmer_complexity(seq, k):
    """Calculate k-mer complexity efficiently"""
    if len(seq) < k:
        return 0
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    return len(set(kmers)) / len(kmers) if kmers else 0

df["kmer_complexity"] = df["seq"].apply(lambda x: kmer_complexity(x, k=5))

all_kmers = []
for seq in df["seq"]:
    k=5
    all_kmers.extend([seq[i:i+k] for i in range(len(seq) - k + 1)])

kmer_series = pd.Series(all_kmers)
top_kmers = kmer_series.value_counts().head(10)

# Overall metrics
total_reads = len(df)
mean_length = df["length"].mean()
median_length = df["length"].median()
overall_gc = df["gc_content"].mean()
mean_quality = df["mean_q"].mean()

# Tabled results
print(df.to_string(max_rows=10))
#print(df.to_string())
print("total read count:\n", total_reads,
      "\nmean length:\n", mean_length,
      "\nmedian length:\n", median_length,
      "\nOverall GC%:\n", overall_gc * 100,
      "\nMean Q-Score:\n", mean_quality)
print("\nTop 10 Overrepresented k-mers:\n", top_kmers)

df.to_csv("parsed_fastq_ecoli.csv", index=False, escapechar='\\')