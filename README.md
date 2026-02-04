# FASTQ-Pandas: Genomic Quality Control Pipeline

##  Overview
This project is a Python-based exploratory data analysis (EDA) tool designed to parse and analyze **FASTQ** sequencing data using **Pandas**. While specialized tools like FastQC are optimal for production-scale environments, this project serves as a technical demonstration of using vectorized string operations and regex to extract biological insights from raw genomic data.

The pipeline processes **E. coli and Shigella** surveillance data from **NCBI BioProject PRJNA315192**, handling real-world artifacts such as binned quality scores and variable read lengths.

---

##  Technical Features
* **Custom FASTQ Parser:** Built a robust iterator-based parser to transform 4-line FASTQ records into a structured DataFrame.
* **Vectorized Genomics Metrics:** Calculated GC-content, N-content, and length distributions using optimized Pandas `.str` accessors.
* **Homopolymer Detection:** Implemented **Regular Expressions (Regex)** to identify and count low-complexity repeats (e.g., AAAA), a primary source of error in long-read sequencing.
* **K-mer Complexity Analysis:** Developed a sliding-window algorithm to calculate sequence diversity and identify overrepresented 5-mers.
* **Quality Decoding:** Decoded **Phred+33 ASCII** encoding to assess the statistical probability of sequencing errors (e.g., `?` = Q29).

---

##  Calculated Metrics

| Metric | logic | Purpose |
| :--- | :--- | :--- |
| **GC Content** | (G+C) / Length | Taxonomic purity & contamination check |
| **N-Count** | N / Length | Basecalling reliability assessment |
| **Mean Q** | (ASCII - 33) | Statistical accuracy of the sequencing run |
| **K-mer Complexity** | Unique / Total | Library diversity & adapter detection |

---

##  Data Insights (PRJNA315192)
* **Binned Quality Scores:** Identified that the dataset uses binned scores where all bases are reported as `?` (Phred 29 / 99.87% accuracy).
* **Variable Read Lengths:** Detected lengths ranging from 35bp to 99bp, indicating the data was likely pre-processed with quality/adapter trimming.
* **Sequence Bias:** Overrepresented k-mers were flagged to identify potential adapter leakage or common motifs in the *E. coli* genome.

---

##  Installation & Usage
1. **Requirements:**
   ```bash
   pip install pandas
2. **Run the parser:**
   ```bash
   python fastq_parser.py
