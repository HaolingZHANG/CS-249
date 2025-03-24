from Bio import SeqIO
from collections import defaultdict
from memory_profiler import profile
from time import perf_counter


class TrieNode:
    def __init__(self):
        self.children = {}
        self.word_ends = set()


def build_trie(reads: dict) \
        -> TrieNode:
    root = TrieNode()
    for read_id, read_seq in reads:
        node = root
        for ch in read_seq:
            node = node.children.setdefault(ch, TrieNode())
        node.word_ends.add(read_id)
    return root


def search_with_mismatch(node, genome, pos, mismatches, max_mismatch, read_len, path_len, matched_ids):
    if path_len == read_len:
        matched_ids.update(node.word_ends)
        return

    if pos >= len(genome):
        return

    genome_char = genome[pos]

    for child_char, child_node in node.children.items():
        mismatch = 0 if child_char == genome_char else 1
        if mismatches + mismatch <= max_mismatch:
            search_with_mismatch(child_node, genome, pos + 1,
                                 mismatches + mismatch, max_mismatch,
                                 read_len, path_len + 1, matched_ids)


@profile
def trie_approximate_match(reads, genomes, max_mismatch):
    read_len = len(reads[0][1])  # assume uniform length
    trie = build_trie(reads)
    read_to_genomes = defaultdict(set)

    for genome_name, sequence in genomes.items():
        start = perf_counter()
        for i in range(len(sequence) - read_len + 1):
            matched = set()
            # start from zero mismatches when entering window match
            search_with_mismatch(trie, sequence, i, 0, max_mismatch, read_len, 0, matched)
            for rid in matched:
                read_to_genomes[rid].add(genome_name)
        end = perf_counter()

        print(genome_name, end - start)
    return read_to_genomes


def load_reads_with_id(fastq_file: str):
    reads = []
    for record in SeqIO.parse(fastq_file, "fastq"):
        cleaned_id = record.id.split('/')[0]
        reads.append((cleaned_id, str(record.seq)))
    return reads


def load_reference_genomes(fasta_files: dict):
    genomes = {}
    for name, path in fasta_files.items():
        sequence = ""
        for record in SeqIO.parse(path, "fasta"):
            sequence += str(record.seq)
        genomes[name] = sequence
    return genomes


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
    fasta_files = {
        "E_coli": "./data/task_1/GCF_000005845.2.fasta",
        "B_subtilis": "./data/task_1/GCF_000009045.1.fasta",
        "P_aeruginosa": "./data/task_1/GCF_000006765.1.fasta",
        "S_aureus": "./data/task_1/GCF_000013425.1.fasta",
        "M_tuberculosis": "./data/task_1/GCF_000195955.2.fasta"
    }

    genomes = load_reference_genomes(fasta_files)

    reads_r1 = load_reads_with_id("./data/task_1/simulated_reads_miseq_10k_R1.fastq")
    reads_r2 = load_reads_with_id("./data/task_1/simulated_reads_miseq_10k_R2.fastq")

    matches_r1 = trie_approximate_match(reads_r1, genomes, max_mismatch=1)
    matches_r2 = trie_approximate_match(reads_r2, genomes, max_mismatch=1)

    merged_matches = merge_paired_matches(matches_r1, matches_r2)
    counts = count_reads_per_genome(merged_matches)

    for genome in fasta_files:
        print(genome + ": " + str(counts[genome]) + " reads matched")


if __name__ == "__main__":
    task()
