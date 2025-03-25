# Assignment 1

## Note on Paired-End Read Processing

As clarified in the communication with the instructor,
paired-end reads can be processed either as independent reads or as paired reads. 
From a biological standpoint, the latter approach is more appropriate, 
as our ultimate goal is to quantify the number of distinct DNA molecules being sequenced. 
Treating paired-end reads as independent units risks inflating counts, 
since both the forward and reverse reads originate from the same DNA fragment. 
If independent reads were sufficient, single-end sequencing or long-read sequencing would have been adequate, 
depending on the target DNA fragment length.

Therefore, in this assignment, we adopt the paired-end perspective.
That is, for a given DNA molecule, even if both the forward and reverse reads are classified as matches, 
we count the pair only once. 
This approach better reflects the experimental intention of metagenomic sequencing. 
That said, our current implementation does not achieve precise quantification of unique molecules, 
as we lack prior information such as whether the forward and reverse reads overlap, 
and if so, the extent of their overlap. 
These factors may influence the actual matched counts. 
A more rigorous strategy would involve merging the forward and reverse reads into a single contiguous sequence 
representing the original DNA molecule before classification. 
Without such merging, our classification is performed on individual reads, 
which may lead to overestimation of the matched counts due to potential double-counting of paired reads.

The implementation of this function is as follows:

```python
from collections import defaultdict

def merge_paired_matches(r1_matches: dict,
                         r2_matches: dict) \
        -> dict:
    merged = defaultdict(set)
    for read_id in set(r1_matches) | set(r2_matches):
        merged[read_id] = r1_matches.get(read_id, set()) | r2_matches.get(read_id, set())
    return merged
```

And the marker of reads direction is removed when reading the FASTQ file:

```python
from Bio import SeqIO

def load_reads_with_id(fastq_file: str):
    reads = []
    for record in SeqIO.parse(fastq_file, "fastq"):
        cleaned_id = record.id.split('/')[0]  # remove here
        reads.append((cleaned_id, str(record.seq)))
    return reads
```

## Task 1

### Task 1.1

In metagenomic classification, a single read may match multiple reference genomes 
due to conserved genomic regions or horizontal gene transfer.
The classification system must efficiently handle these multi-matching reads across an entire collection.

Between the two classes of string matching algorithms discussed in class:

- Aho-Corasick builds an automaton from all query patterns (reads) and scans the reference text(s) once.
- Suffix-based structures (suffix trees / arrays) build an index on the reference genome and support multiple queries.

We recommend using **Aho-Corasick** for this task. The key reason:  
> We are classifying an entire collection of reads, many of which may share identical matches across different genomes.

Why Aho-Corasick fits this goal better:

- Efficiently handles batch matching of thousands of reads at once.
- Allows us to scan the entire reference genome set in a single linear pass per genome.
- Supports multi-genome matching, as each match is annotated with the genome of origin.
- Avoids building and storing large suffix structures for all genomes (memory-efficient).

Identical sequences in multiple genomes are biologically expected (**biological justification**) due to:

- Shared evolutionary origins
- Conserved functional regions (e.g., rRNAs, core enzymes)
- Horizontal gene transfer across species

Hence, we must design classification rules that reflect this biological ambiguity.

When a read matches `N` different reference genomes, we assign `1/N` count to each matched organism.

- **Pros:** Reflects biological uncertainty fairly
- **Cons:** Slight underestimation of low-abundance species

This approach maintains a balanced and interpretable abundance estimation.

For the **computational efficiency**:

- Build an Aho-Corasick automaton once for all reads.
- Iterate over each reference genome, applying the automaton to collect matches.
- Store matches as `(read ID → genome ID)` pairs.
- Aggregate counts using fractional weights for multi-matched reads.

### Task 1.2

In this task, we implemented an exact matching-based read classification pipeline using a custom Aho-Corasick automaton.
All details are shown in the [Python script](https://github.com/HaolingZHANG/CS-249/tree/main/assign1/task_1.2.py) 
(**consistent with the task requirements**).

Specifically, we processed Illumina paired-end sequencing data by separately analyzing R1 and R2 reads. 
To ensure accurate read-level classification, 
we normalized read IDs by removing /1 and /2 suffixes and merged the matches from both ends based on the shared read ID.
Each read was counted once per genome if either of its ends matched. 
This approach reflects realistic metagenomic analysis scenarios where partial information (from one end) 
is sufficient to suggest the presence of an organism. Our Aho-Corasick implementation was done from scratch 
to demonstrate understanding of the trie construction, failure links, and pattern matching.

We evaluated our classifier on five bacterial reference genomes and obtained the following match counts:

| Species         | matched counts |
|-----------------|----------------|
| E. coli         | 3000           |
| B. subtilis     | 500            |
| P. aeruginosa   | 500            |
| S. aureus       | 509            |
| M. tuberculosis | 500            |

### Task 1.3

We implemented a multi-pattern approximate matcher based on a Trie data structure with bounded Hamming distance. 
Each reference genome is scanned with a sliding window, 
and a backtracking DFS is used to traverse the Trie while allowing at most 1 mismatch. 
This method retains the core advantages of Aho-Corasick (multi-pattern indexing) 
while supporting controlled error tolerance.
All details are shown in the [Python script](https://github.com/HaolingZHANG/CS-249/tree/main/assign1/task_1.3.py) 
(**consistent with the task requirements**).

We evaluated our classifier on five bacterial reference genomes and obtained the following match counts:

| Species         | matched counts |
|-----------------|----------------|
| E. coli         | 771            |
| B. subtilis     | 132            |
| P. aeruginosa   | 130            |
| S. aureus       | 143            |
| M. tuberculosis | 111            |

### Task 1.4

I implemented a BLAST-based read mapping method that allows up to one mismatch. 
The workflow consists of three main steps
(all details are shown in the [Python script](https://github.com/HaolingZHANG/CS-249/tree/main/assign1/task_1.4.py)):

1. **Database Construction**: For each reference genome, I used `makeblastdb` to create a nucleotide BLAST database.
2. **Read Alignment**: I performed sequence alignment using `blastn`, outputting the results in tabular (TSV) format using `-outfmt 6`.
3. **Post-Processing**: I parsed the resulting TSV file and filtered the alignments, retaining only those with a mismatch count of one or fewer.

This approach enables a fair comparison with other methods by applying a consistent mismatch threshold.

In addition, execution time and memory usage are estimated by `time` and `memory_profiler`, respectively. 

Based on the obtained results (saved in `./data/task_1`), the summary of BLAST-based method is as follows:

| Species         | matched counts | runtime | memory usage |
|-----------------|----------------|---------|--------------|
| E. coli         | 1004           | 5.89    | 59.9 MiB     |
| B. subtilis     | 162            | 5.09    | 60.0 MiB     |
| P. aeruginosa   | 169            | 5.07    | 60.0 MiB     |
| S. aureus       | 190            | 5.33    | 60.0 MiB     |
| M. tuberculosis | 148            | 5.21    | 60.0 MiB     |

The results of the Trie-based methods are presented in Task 1.3. 
Additionally, both methods exhibit similar performance, 
with runtimes of approximately 10 -- 20 minutes and memory usage around 120 -- 240 MiB.

Although both methods apply a mismatch ≤ 1 criterion, 
BLAST detects local alignments, often matching only a portion of each read. 
In contrast, our Trie-based matcher requires the entire read to align (i.e., global alignment), 
resulting in fewer but stricter matches. 
This fundamental difference explains why BLAST returns more hits even under the same mismatch constraint.

## Task 2

### Task 2.1

To store the k-mer index, we used a nested Python dictionary. 
The outer dictionary maps each unique k-mer string (of length 31) to 
an inner dictionary that records how many times this k-mer appears in each of the five reference genomes. 
This structure supports efficient lookup, insertion, and genome-specific aggregation, 
which is essential for tracking metagenomic signals across multiple genomes. 
Additionally, using `defaultdict` enables safe and concise increment operations during indexing.
All details are shown in the [Python script](https://github.com/HaolingZHANG/CS-249/tree/main/assign1/task_2.1.py).

The total number of unique k-mers in the index is **22104695**.

The number of theoretically possible k-mers of length 31 is **4^31 = 4.61e+18**,
so totoal collected number is only a tiny fraction of the theoretical maximum.

The massive gap between the theoretical and observed number of k-mers is expected. 
While 4.6 quintillion combinations are mathematically possible,
real biological genomes only contain a tiny,
biologically meaningful subset of these. 
The actual number of k-mers is constrained by:
- The total size of the input genomes (roughly 20–25 million base pairs combined);
- Biological redundancy and sequence composition biases (e.g., GC content);
- The fact that many potential k-mers (especially highly random or low-complexity ones) simply do not occur in nature;
- The exclusion of k-mers containing ambiguous nucleotides like `N`.
Thus, the observed index size (~22 million) represents the true diversity of k-mers 
present across these five bacterial genomes, 
and the extremely low coverage ratio (~4.79e-12) is consistent with expectations in genomics.

### Task 2.2

We implemented a classification method based on exact k-mer matches (with k=31). 
For each read, all overlapping 31-mers were extracted and checked against the k-mer index constructed 
from the five reference genomes. 
If at least one k-mer from the read matched a genome, 
the read was considered to be associated with that genome. 
This procedure was applied separately to both R1 and R2 reads, 
and the results were merged based on the shared read ID.
All details are shown in the [Python script](https://github.com/HaolingZHANG/CS-249/tree/main/assign1/task_2.2.py).

The final read classification results were as follows:

| Species         | k-mer Match | Trie Match |
|-----------------|-------------|------------|
| E. coli         | 3003        | 3000       |
| B. subtilis     | 514         | 500        |
| P. aeruginosa   | 533         | 500        |
| S. aureus       | 526         | 509        |
| M. tuberculosis | 506         | 500        |

The results are largely consistent with the previous Trie-based classification method that allowed up to one mismatch. 
However, the k-mer approach yielded slightly higher match counts across all genomes. 
This discrepancy arises because the k-mer method considers **any exact k-mer match**, 
even if only part of a read aligns to a genome, 
while the Trie method attempts to match entire reads with limited mismatches. 
In essence, the k-mer method is more sensitive to partial matches, while the Trie method is stricter.

When a k-mer is found in multiple genomes, we associate the corresponding read with **all matching genomes**. 
Similarly, if multiple k-mers from a single read match different genomes, 
the read is considered to be matched to each of those genomes. 
This approach ensures that ambiguous reads are not prematurely excluded and reflects the biological possibility 
of conserved sequences across different organisms.

### Task 2.3

To reduce memory consumption while maintaining classification accuracy, 
we implemented a minimizer-based index using parameters **k = 31** and **window size w = 10**. 
Instead of storing every overlapping k-mer,
we only retained the lexicographically smallest k-mer (the "minimizer") within each window of `w + k - 1` bases. 
This significantly reduced the number of stored entries.
All details are shown in the [Python script](https://github.com/HaolingZHANG/CS-249/tree/main/assign1/task_2.3.py).

The memory efficiency is reduced as follows:

| Index Type       | Unique Entries | Memory Usage |
|------------------|----------------|--------------|
| Full k-mer Index | 22104664       | ~671 MiB     |
| Minimizer Index  | 4515766        | ~168 MiB     |

The minimizer index reduced the number of stored k-mers by nearly **80%**, and memory usage by approximately **75%**.

The classification accuracy comparison is as follows:

| Species         | Full k-mer Match | Minimizer Match |
|-----------------|------------------|-----------------|
| E. coli         | 3003             | 3003            |
| B. subtilis     | 514              | 514             |
| P. aeruginosa   | 533              | 529             |
| S. aureus       | 526              | 524             |
| M. tuberculosis | 506              | 506             |

Despite the reduction in index size, the classification accuracy remained **virtually identical**, 
with only negligible differences (<= 4 reads per genome). 
This demonstrates that minimizers can preserve nearly all useful signal while greatly reducing storage overhead.

Minimizer-based indexing is a space-efficient alternative to full k-mer indexing in metagenomic classification. 
It achieves comparable accuracy with significantly lower memory requirements,
making it ideal for scaling to large reference datasets.

## Task 3

### Task 3.1

Leveraging Python's capability as a glue language, 
we utilized it to automate the construction of a Kraken2 database based on the five provided FASTA sequences
(detailed in the [Python script](https://github.com/HaolingZHANG/CS-249/tree/main/assign1/task_3.1.py)).

Notably, to construct the database successfully, the head of each FASTA sequence should be modifed, like:

```text
>seq1|kraken:taxid|562
```

Therefore, we successfully built a Kraken2 custom database using the correct NCBI taxonomy IDs 
for the five reference genomes:
- *Escherichia coli* (562)
- *Pseudomonas aeruginosa* (287)
- *Staphylococcus aureus* (1280)
- *Bacillus subtilis* (1423)
- *Mycobacterium tuberculosis* (1773)

The database construction log is:
```text
Masking low-complexity regions of new file... done.
Added "/home/horus/Desktop/Projects/task_3/GCF_000005845.2.fasta" to library (/home/horus/Desktop/Projects/database)
Masking low-complexity regions of new file... done.
Added "/home/horus/Desktop/Projects/task_3/GCF_000009045.1.fasta" to library (/home/horus/Desktop/Projects/database)
Masking low-complexity regions of new file... done.
Added "/home/horus/Desktop/Projects/task_3/GCF_000006765.1.fasta" to library (/home/horus/Desktop/Projects/database)
Masking low-complexity regions of new file... done.
Added "/home/horus/Desktop/Projects/task_3/GCF_000013425.1.fasta" to library (/home/horus/Desktop/Projects/database)
Masking low-complexity regions of new file... done.
Added "/home/horus/Desktop/Projects/task_3/GCF_000195955.2.fasta" to library (/home/horus/Desktop/Projects/database)
Creating sequence ID to taxonomy ID map (step 1)...
Sequence ID to taxonomy ID map complete. [0.014s]
Estimating required capacity (step 2)...
Estimated hash table requirement: 40354376 bytes
Capacity estimation complete. [0.485s]
Building database files (step 3)...
Taxonomy parsed and converted.
CHT created with 6 bits reserved for taxid.
Completed processing of 5 sequences, 22354555 bp
Writing data to disk...  complete.
Database files completed. [9.346s]
Database construction complete. [Total: 9.866s]
```

Then, as shown in Line 12 -- 18 in the [Python script](https://github.com/HaolingZHANG/CS-249/tree/main/assign1/task_3.1.py),
the classification task is executed with the following log:

```text
Loading database information... done.
5000 sequences (1.25 Mbp) processed in 0.087s (3450.9 Kseq/m, 862.72 Mbp/m).
  4997 sequences classified (99.94%)
  3 sequences unclassified (0.06%)
```

Kraken2 produced species-level classifications with high accuracy. 
Only 3 out of 5000 reads remained unclassified.

The differences between species abundance can be compared here:

| Species         | Kraken2 | Minimizer | Full k-mer |
|-----------------|---------|-----------|------------|
| E. coli         | 2993    | 3003      | 3003       |
| B. subtilis     | 500     | 514       | 514        |
| P. aeruginosa   | 500     | 529       | 533        |
| S. aureus       | 500     | 524       | 526        |
| M. tuberculosis | 497     | 506       | 506        |
| Unclassified    | 3       | —         | —          |


### Task 3.2

As shown in the [Python script](https://github.com/HaolingZHANG/CS-249/tree/main/assign1/task_3.2.py),
we utilized it to construct Kraken2 database with Standard-8 DB, 
and then, downloaded metagenomic samples with the following log.

```text
2025-03-24T10:21:22 prefetch.3.2.1: 1) Resolving 'SRR11412973'...
2025-03-24T10:21:26 prefetch.3.2.1: Current preference is set to retrieve SRA Normalized Format files with full base quality scores
2025-03-24T10:21:27 prefetch.3.2.1: 1) Downloading 'SRR11412973'...
2025-03-24T10:21:27 prefetch.3.2.1:  SRA Normalized Format file is being retrieved
2025-03-24T10:21:27 prefetch.3.2.1:  Downloading via HTTPS...
2025-03-24T10:32:23 prefetch.3.2.1:  HTTPS download succeed
2025-03-24T10:32:28 prefetch.3.2.1:  'SRR11412973' is valid: 2359575904 bytes were streamed from 2359565879
2025-03-24T10:32:28 prefetch.3.2.1: 1) 'SRR11412973' was downloaded successfully
2025-03-24T10:32:28 prefetch.3.2.1: 1) Resolving 'SRR11412973's dependencies...
2025-03-24T10:32:28 prefetch.3.2.1: 'SRR11412973' has 0 unresolved dependencies
spots read      : 19,515,475
reads read      : 39,030,950
reads written   : 39,030,950
2025-03-24T10:32:28 prefetch.3.2.1: 1) Resolving 'SRR11412976'...
2025-03-24T10:32:32 prefetch.3.2.1: Current preference is set to retrieve SRA Normalized Format files with full base quality scores
2025-03-24T10:32:33 prefetch.3.2.1: 1) Downloading 'SRR11412976'...
2025-03-24T10:32:33 prefetch.3.2.1:  SRA Normalized Format file is being retrieved
2025-03-24T10:32:33 prefetch.3.2.1:  Downloading via HTTPS...
2025-03-24T10:38:22 prefetch.3.2.1:  HTTPS download succeed
2025-03-24T10:38:26 prefetch.3.2.1:  'SRR11412976' is valid: 2038133951 bytes were streamed from 2038132259
2025-03-24T10:38:26 prefetch.3.2.1: 1) 'SRR11412976' was downloaded successfully
2025-03-24T10:38:26 prefetch.3.2.1: 1) Resolving 'SRR11412976's dependencies...
2025-03-24T10:38:26 prefetch.3.2.1: 'SRR11412976' has 0 unresolved dependencies
spots read      : 16,558,897
reads read      : 33,117,794
reads written   : 33,117,794
2025-03-24T10:38:26 prefetch.3.2.1: 1) Resolving 'SRR11412979'...
2025-03-24T10:38:30 prefetch.3.2.1: Current preference is set to retrieve SRA Normalized Format files with full base quality scores
2025-03-24T10:38:31 prefetch.3.2.1: 1) Downloading 'SRR11412979'...
2025-03-24T10:38:31 prefetch.3.2.1:  SRA Normalized Format file is being retrieved
2025-03-24T10:38:31 prefetch.3.2.1:  Downloading via HTTPS...
2025-03-24T10:45:43 prefetch.3.2.1:  HTTPS download succeed
2025-03-24T10:45:49 prefetch.3.2.1:  'SRR11412979' is valid: 2435341599 bytes were streamed from 2435334699
2025-03-24T10:45:49 prefetch.3.2.1: 1) 'SRR11412979' was downloaded successfully
2025-03-24T10:45:49 prefetch.3.2.1: 1) Resolving 'SRR11412979's dependencies...
2025-03-24T10:45:49 prefetch.3.2.1: 'SRR11412979' has 0 unresolved dependencies
spots read      : 19,795,127
reads read      : 39,590,254
reads written   : 39,590,254
2025-03-24T10:45:49 prefetch.3.2.1: 1) Resolving 'SRR11412980'...
2025-03-24T10:45:51 prefetch.3.2.1: Current preference is set to retrieve SRA Normalized Format files with full base quality scores
2025-03-24T10:45:53 prefetch.3.2.1: 1) Downloading 'SRR11412980'...
2025-03-24T10:45:53 prefetch.3.2.1:  SRA Normalized Format file is being retrieved
2025-03-24T10:45:53 prefetch.3.2.1:  Downloading via HTTPS...
2025-03-24T10:52:08 prefetch.3.2.1:  HTTPS download succeed
2025-03-24T10:52:13 prefetch.3.2.1:  'SRR11412980' is valid: 2502947164 bytes were streamed from 2502936087
2025-03-24T10:52:13 prefetch.3.2.1: 1) 'SRR11412980' was downloaded successfully
2025-03-24T10:52:13 prefetch.3.2.1: 1) Resolving 'SRR11412980's dependencies...
2025-03-24T10:52:13 prefetch.3.2.1: 'SRR11412980' has 0 unresolved dependencies
spots read      : 21,457,804
reads read      : 42,915,608
reads written   : 42,915,608
2025-03-24T10:52:13 prefetch.3.2.1: 1) Resolving 'SRR11412984'...
2025-03-24T10:52:16 prefetch.3.2.1: Current preference is set to retrieve SRA Normalized Format files with full base quality scores
2025-03-24T10:52:17 prefetch.3.2.1: 1) Downloading 'SRR11412984'...
2025-03-24T10:52:17 prefetch.3.2.1:  SRA Normalized Format file is being retrieved
2025-03-24T10:52:17 prefetch.3.2.1:  Downloading via HTTPS...
2025-03-24T10:58:43 prefetch.3.2.1:  HTTPS download succeed
2025-03-24T10:58:48 prefetch.3.2.1:  'SRR11412984' is valid: 2424254138 bytes were streamed from 2424238647
2025-03-24T10:58:48 prefetch.3.2.1: 1) 'SRR11412984' was downloaded successfully
2025-03-24T10:58:48 prefetch.3.2.1: 1) Resolving 'SRR11412984's dependencies...
2025-03-24T10:58:48 prefetch.3.2.1: 'SRR11412984' has 0 unresolved dependencies
spots read      : 19,866,883
reads read      : 39,733,766
reads written   : 39,733,766
2025-03-24T10:58:48 prefetch.3.2.1: 1) Resolving 'SRR21907296'...
2025-03-24T10:58:52 prefetch.3.2.1: Current preference is set to retrieve SRA Normalized Format files with full base quality scores
2025-03-24T10:58:53 prefetch.3.2.1: 1) Downloading 'SRR21907296'...
2025-03-24T10:58:53 prefetch.3.2.1:  SRA Normalized Format file is being retrieved
2025-03-24T10:58:53 prefetch.3.2.1:  Downloading via HTTPS...
2025-03-24T10:59:20 prefetch.3.2.1:  HTTPS download succeed
2025-03-24T10:59:20 prefetch.3.2.1:  'SRR21907296' is valid: 147100412 bytes were streamed from 147091991
2025-03-24T10:59:20 prefetch.3.2.1: 1) 'SRR21907296' was downloaded successfully
2025-03-24T10:59:20 prefetch.3.2.1: 1) Resolving 'SRR21907296's dependencies...
2025-03-24T10:59:20 prefetch.3.2.1: 'SRR21907296' has 0 unresolved dependencies
spots read      : 877,211
reads read      : 1,754,422
reads written   : 1,754,420
reads 0-length  : 2
2025-03-24T10:59:20 prefetch.3.2.1: 1) Resolving 'SRR21907303'...
2025-03-24T10:59:23 prefetch.3.2.1: Current preference is set to retrieve SRA Normalized Format files with full base quality scores
2025-03-24T10:59:24 prefetch.3.2.1: 1) Downloading 'SRR21907303'...
2025-03-24T10:59:24 prefetch.3.2.1:  SRA Normalized Format file is being retrieved
2025-03-24T10:59:24 prefetch.3.2.1:  Downloading via HTTPS...
2025-03-24T11:00:49 prefetch.3.2.1:  HTTPS download succeed
2025-03-24T11:00:50 prefetch.3.2.1:  'SRR21907303' is valid: 207682426 bytes were streamed from 207670807
2025-03-24T11:00:50 prefetch.3.2.1: 1) 'SRR21907303' was downloaded successfully
2025-03-24T11:00:50 prefetch.3.2.1: 1) Resolving 'SRR21907303's dependencies...
2025-03-24T11:00:50 prefetch.3.2.1: 'SRR21907303' has 0 unresolved dependencies
spots read      : 1,113,099
reads read      : 2,226,198
reads written   : 2,226,194
reads 0-length  : 4
2025-03-24T11:00:50 prefetch.3.2.1: 1) Resolving 'SRR21907307'...
2025-03-24T11:00:52 prefetch.3.2.1: Current preference is set to retrieve SRA Normalized Format files with full base quality scores
2025-03-24T11:00:53 prefetch.3.2.1: 1) Downloading 'SRR21907307'...
2025-03-24T11:00:53 prefetch.3.2.1:  SRA Normalized Format file is being retrieved
2025-03-24T11:00:53 prefetch.3.2.1:  Downloading via HTTPS...
2025-03-24T11:01:06 prefetch.3.2.1:  HTTPS download succeed
2025-03-24T11:01:07 prefetch.3.2.1:  'SRR21907307' is valid: 104846027 bytes were streamed from 104833539
2025-03-24T11:01:07 prefetch.3.2.1: 1) 'SRR21907307' was downloaded successfully
2025-03-24T11:01:07 prefetch.3.2.1: 1) Resolving 'SRR21907307's dependencies...
2025-03-24T11:01:07 prefetch.3.2.1: 'SRR21907307' has 0 unresolved dependencies
spots read      : 523,931
reads read      : 1,047,862
reads written   : 1,047,862
2025-03-24T11:01:07 prefetch.3.2.1: 1) Resolving 'SRR21907332'...
2025-03-24T11:01:09 prefetch.3.2.1: Current preference is set to retrieve SRA Normalized Format files with full base quality scores
2025-03-24T11:01:10 prefetch.3.2.1: 1) Downloading 'SRR21907332'...
2025-03-24T11:01:10 prefetch.3.2.1:  SRA Normalized Format file is being retrieved
2025-03-24T11:01:10 prefetch.3.2.1:  Downloading via HTTPS...
2025-03-24T11:01:31 prefetch.3.2.1:  HTTPS download succeed
2025-03-24T11:01:31 prefetch.3.2.1:  'SRR21907332' is valid: 75446983 bytes were streamed from 75438617
2025-03-24T11:01:31 prefetch.3.2.1: 1) 'SRR21907332' was downloaded successfully
2025-03-24T11:01:31 prefetch.3.2.1: 1) Resolving 'SRR21907332's dependencies...
2025-03-24T11:01:31 prefetch.3.2.1: 'SRR21907332' has 0 unresolved dependencies
spots read      : 546,485
reads read      : 1,092,970
reads written   : 1,092,970
2025-03-24T11:01:31 prefetch.3.2.1: 1) Resolving 'SRR21907330'...
2025-03-24T11:01:33 prefetch.3.2.1: Current preference is set to retrieve SRA Normalized Format files with full base quality scores
2025-03-24T11:01:34 prefetch.3.2.1: 1) Downloading 'SRR21907330'...
2025-03-24T11:01:34 prefetch.3.2.1:  SRA Normalized Format file is being retrieved
2025-03-24T11:01:34 prefetch.3.2.1:  Downloading via HTTPS...
2025-03-24T11:01:51 prefetch.3.2.1:  HTTPS download succeed
2025-03-24T11:01:51 prefetch.3.2.1:  'SRR21907330' is valid: 29053515 bytes were streamed from 29040191
2025-03-24T11:01:51 prefetch.3.2.1: 1) 'SRR21907330' was downloaded successfully
2025-03-24T11:01:51 prefetch.3.2.1: 1) Resolving 'SRR21907330's dependencies...
2025-03-24T11:01:51 prefetch.3.2.1: 'SRR21907330' has 0 unresolved dependencies
spots read      : 181,195
reads read      : 362,390
reads written   : 362,390
```

Then, as shown in Line 14 -- 20 in the [Python script](https://github.com/HaolingZHANG/CS-249/tree/main/assign1/task_3.2.py),
multiple classification tasks are executed with the following log:

```text
Loading database information... done.
19515475 sequences (5850.48 Mbp) processed in 285.257s (4104.8 Kseq/m, 1230.57 Mbp/m).
  12372111 sequences classified (63.40%)
  7143364 sequences unclassified (36.60%)
Loading database information... done.
16558897 sequences (4963.68 Mbp) processed in 253.729s (3915.7 Kseq/m, 1173.77 Mbp/m).
  14951645 sequences classified (90.29%)
  1607252 sequences unclassified (9.71%)
Loading database information... done.
19795127 sequences (5933.88 Mbp) processed in 294.441s (4033.8 Kseq/m, 1209.18 Mbp/m).
  14705161 sequences classified (74.29%)
  5089966 sequences unclassified (25.71%)
Loading database information... done.
21457804 sequences (6418.08 Mbp) processed in 287.900s (4471.9 Kseq/m, 1337.56 Mbp/m).
  9840445 sequences classified (45.86%)
  11617359 sequences unclassified (54.14%)
Loading database information... done.
19866883 sequences (5946.67 Mbp) processed in 293.804s (4057.2 Kseq/m, 1214.41 Mbp/m).
  13746688 sequences classified (69.19%)
  6120195 sequences unclassified (30.81%)
Loading database information... done.
877209 sequences (264.92 Mbp) processed in 12.015s (4380.5 Kseq/m, 1322.92 Mbp/m).
  743460 sequences classified (84.75%)
  133749 sequences unclassified (15.25%)
Loading database information... done.
1113096 sequences (336.15 Mbp) processed in 15.482s (4313.8 Kseq/m, 1302.78 Mbp/m).
  1059734 sequences classified (95.21%)
  53362 sequences unclassified (4.79%)
Loading database information... done.
523931 sequences (158.23 Mbp) processed in 6.736s (4667.0 Kseq/m, 1409.44 Mbp/m).
  319731 sequences classified (61.03%)
  204200 sequences unclassified (38.97%)
Loading database information... done.
546485 sequences (159.06 Mbp) processed in 7.369s (4449.5 Kseq/m, 1295.08 Mbp/m).
  517185 sequences classified (94.64%)
  29300 sequences unclassified (5.36%)
Loading database information... done.
181195 sequences (54.72 Mbp) processed in 2.665s (4079.4 Kseq/m, 1231.99 Mbp/m).
  165784 sequences classified (91.49%)
  15411 sequences unclassified (8.51%)
```

Based on the obtained results (saved in `./data/task_3`), the summary table is:

| Sample      | Top Species                                           | Matched Reads |
|-------------|-------------------------------------------------------|---------------|
| SRR11412973 | Phascolarctobacterium faecium                         | 1,232,441     |
| SRR11412976 | Phocaeicola vulgatus                                  | 2,696,037     |
| SRR11412979 | Segatella copri                                       | 8,042,501     |
| SRR11412980 | Bacteroides uniformis                                 | 863,249       |
| SRR11412984 | Segatella copri                                       | 3,534,192     |
| SRR21907296 | Severe acute respiratory syndrome-related coronavirus | 661,508       |
| SRR21907303 | Severe acute respiratory syndrome-related coronavirus | 530,048       |
| SRR21907307 | Severe acute respiratory syndrome-related coronavirus | 240,869       |
| SRR21907330 | Severe acute respiratory syndrome-related coronavirus | 160,535       |
| SRR21907332 | Severe acute respiratory syndrome-related coronavirus | 515,714       |

Based on the results, we can clearly distinguish the sample groups using their top species profiles.

The SRR1141xxxx samples can be identified as gut samples, 
as their dominant species are typical members of the human gut microbiota. 
For instance, Segatella copri and Phocaeicola vulgatus are common commensal bacteria in the human colon,
and Bacteroides uniformis is a well-known anaerobe frequently found in healthy intestinal environments.

In contrast, the SRR2190xxxx samples are clearly wastewater samples, 
as all are dominated by Severe acute respiratory syndrome-related coronavirus (SARS-related coronavirus). 
This suggests that the microbial composition is overwhelmingly influenced by viral sequences, 
with much lower microbial diversity compared to gut samples—characteristic of a typical wastewater profile.

Therefore, we can achieve a highly accurate classification of the samples into two distinct groups
based on the top species reported by Kraken2. 
These two classes reliably correspond to their original environmental sources (gut versus wastewater),
offering a biologically meaningful and data-driven separation strategy.

