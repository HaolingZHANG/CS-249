from Bio import SeqIO
from collections import defaultdict
from sys import getsizeof


def extract_kmers(sequence, k):
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]


def build_kmer_index(fasta_files, k):
    index = defaultdict(lambda: defaultdict(int))

    for genome_name, filepath in fasta_files.items():
        print("Indexing genome: {}".format(genome_name))
        for record in SeqIO.parse(filepath, "fasta"):
            seq = str(record.seq).upper()
            for kmer in extract_kmers(seq, k):
                if "N" not in kmer:
                    index[kmer][genome_name] += 1

    return index


def load_reads_with_id(fastq_file):
    reads = []
    for record in SeqIO.parse(fastq_file, "fastq"):
        read_id = record.id.split('/')[0]
        reads.append((read_id, str(record.seq)))
    return reads


def classify_reads_by_kmer(reads, kmer_index, k):
    read_to_genomes = defaultdict(set)

    for read_id, seq in reads:
        matched = set()
        for kmer in extract_kmers(seq, k):
            if kmer in kmer_index:
                matched.update(kmer_index[kmer].keys())
        if matched:
            read_to_genomes[read_id] = matched

    return read_to_genomes


def merge_paired_matches(r1_matches: dict,
                         r2_matches: dict) \
        -> dict:
    merged = defaultdict(set)
    for read_id in set(r1_matches) | set(r2_matches):
        merged[read_id] = r1_matches.get(read_id, set()) | r2_matches.get(read_id, set())
    return merged


def count_reads_per_genome(merged_matches: dict) \
        -> dict:
    genome_counts = defaultdict(int)
    for matched_genomes in merged_matches.values():
        for genome in matched_genomes:
            genome_counts[genome] += 1
    return genome_counts


def task():
    k = 31

    fasta_files = {
        "E_coli": "./data/task_2/GCF_000005845.2.fasta",
        "B_subtilis": "./data/task_2/GCF_000009045.1.fasta",
        "P_aeruginosa": "./data/task_2/GCF_000006765.1.fasta",
        "S_aureus": "./data/task_2/GCF_000013425.1.fasta",
        "M_tuberculosis": "./data/task_2/GCF_000195955.2.fasta"
    }
    kmer_index = build_kmer_index(fasta_files, k)
    print("Total unique indices: {}".format(len(kmer_index)))
    print("Estimated index memory (Python dict overhead not included): ~{} bytes".format(getsizeof(kmer_index)))

    reads_r1 = load_reads_with_id("./data/task_2/simulated_reads_no_errors_10k_R1.fastq")
    reads_r2 = load_reads_with_id("./data/task_2/simulated_reads_no_errors_10k_R2.fastq")

    matches_r1 = classify_reads_by_kmer(reads_r1, kmer_index, k)
    matches_r2 = classify_reads_by_kmer(reads_r2, kmer_index, k)

    merged_matches = merge_paired_matches(matches_r1, matches_r2)
    counts = count_reads_per_genome(merged_matches)

    for genome in fasta_files:
        print("{}: {} reads matched".format(genome, counts[genome]))


if __name__ == "__main__":
    task()
