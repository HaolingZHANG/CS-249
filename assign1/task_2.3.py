from Bio import SeqIO
from collections import defaultdict
from sys import getsizeof


def extract_kmers(sequence, k):
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]


def extract_minimizers(sequence, k, w):
    minimizers = set()
    if len(sequence) < w + k - 1:
        return minimizers
    for i in range(len(sequence) - (w + k - 1) + 1):
        window = sequence[i:i + w + k - 1]
        kmers = [window[j:j + k] for j in range(w)]
        min_kmer = min(kmers)
        if "N" not in min_kmer:
            minimizers.add(min_kmer)
    return minimizers


def build_minimizer_index(fasta_files, k, w):
    index = defaultdict(lambda: defaultdict(int))  # minimizer → {genome: count}
    for genome_name, filepath in fasta_files.items():
        print("Indexing (minimizer) genome: {}".format(genome_name))
        for record in SeqIO.parse(filepath, "fasta"):
            seq = str(record.seq).upper()
            minimizers = extract_minimizers(seq, k, w)
            for minimizer in minimizers:
                index[minimizer][genome_name] += 1
    return index


def load_reads_with_id(fastq_file):
    reads = []
    for record in SeqIO.parse(fastq_file, "fastq"):
        read_id = record.id.split('/')[0]
        reads.append((read_id, str(record.seq)))
    return reads


def classify_reads_by_kmer(reads, index, k):
    read_to_genomes = defaultdict(set)

    for read_id, seq in reads:
        matched = set()
        for kmer in extract_kmers(seq, k):
            if kmer in index:
                matched.update(index[kmer].keys())
        if matched:
            read_to_genomes[read_id] = matched

    return read_to_genomes


def merge_paired_matches(r1_matches, r2_matches):
    merged = defaultdict(set)
    for read_id in set(r1_matches) | set(r2_matches):
        merged[read_id] = r1_matches.get(read_id, set()) | r2_matches.get(read_id, set())
    return merged


def count_reads_per_genome(merged_matches):
    genome_counts = defaultdict(int)
    for matched_genomes in merged_matches.values():
        for genome in matched_genomes:
            genome_counts[genome] += 1
    return genome_counts


def task():
    k, w = 31, 10
    fasta_files = {
        "E_coli": "./data/task_2/GCF_000005845.2.fasta",
        "B_subtilis": "./data/task_2/GCF_000009045.1.fasta",
        "P_aeruginosa": "./data/task_2/GCF_000006765.1.fasta",
        "S_aureus": "./data/task_2/GCF_000013425.1.fasta",
        "M_tuberculosis": "./data/task_2/GCF_000195955.2.fasta"
    }

    print("Building minimizer index with k = " + str(k) + ", w = " + str(w))
    minimizer_index = build_minimizer_index(fasta_files, k, w)
    print("Total unique minimizers: " + str(len(minimizer_index)))
    print("Estimated index memory (Python dict overhead not included): ~" + str(getsizeof(minimizer_index)) + " bytes")

    reads_r1 = load_reads_with_id("./data/task_2/simulated_reads_no_errors_10k_R1.fastq")
    reads_r2 = load_reads_with_id("./data/task_2/simulated_reads_no_errors_10k_R2.fastq")

    matches_r1 = classify_reads_by_kmer(reads_r1, minimizer_index, k)
    matches_r2 = classify_reads_by_kmer(reads_r2, minimizer_index, k)

    merged_matches = merge_paired_matches(matches_r1, matches_r2)
    counts = count_reads_per_genome(merged_matches)

    print("Classification results using minimizers:")
    for genome in fasta_files:
        print(genome + ": " + str(counts[genome]) + " reads matched.")


if __name__ == "__main__":
    task()
