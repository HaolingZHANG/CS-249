from Bio import SeqIO
from collections import defaultdict, deque


# ===== Aho-Corasick Automaton =====
class Node:
    def __init__(self):
        self.children = dict()
        self.fail = None
        self.output = set()


class AhoCorasickAutomaton:
    def __init__(self):
        self.root = Node()

    def add_word(self,
                 word: str,
                 read_id: str):
        node = self.root
        for ch in word:
            node = node.children.setdefault(ch, Node())
        node.output.add(read_id)

    def build(self):
        queue = deque()
        for child in self.root.children.values():
            child.fail = self.root
            queue.append(child)

        while queue:
            current = queue.popleft()
            for ch, child in current.children.items():
                queue.append(child)
                fail = current.fail
                while fail and ch not in fail.children:
                    fail = fail.fail
                child.fail = fail.children[ch] if fail and ch in fail.children else self.root
                child.output |= child.fail.output

    def search(self,
               text: str):
        node = self.root
        results = defaultdict(set)
        for index, ch in enumerate(text):
            while node and ch not in node.children:
                node = node.fail
            node = node.children[ch] if node and ch in node.children else self.root
            for read_id in node.output:
                results[read_id].add(index)
        return results


def load_reads_with_id(fastq_file: str):
    reads = []
    for record in SeqIO.parse(fastq_file, "fastq"):
        cleaned_id = record.id.split('/')[0]
        reads.append((cleaned_id, str(record.seq)))
    return reads


def load_reference_genomes(fasta_files: dict) \
        -> dict:
    genomes = {}
    for name, path in fasta_files.items():
        sequence = ""
        for record in SeqIO.parse(path, "fasta"):
            sequence += str(record.seq)
        genomes[name] = sequence
    return genomes


def build_automaton_with_ids(reads: list) \
        -> AhoCorasickAutomaton:
    ac = AhoCorasickAutomaton()
    for read_id, read_seq in reads:
        ac.add_word(read_seq, read_id)
    ac.build()

    return ac


def match_reads_by_id(automaton: AhoCorasickAutomaton,
                      genomes: dict):
    read_to_genomes = defaultdict(set)
    for genome_name, sequence in genomes.items():
        matches = automaton.search(sequence)
        for read_id in matches:
            read_to_genomes[read_id].add(genome_name)
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
    fasta_files = {
        "E_coli": "./data/task_1/GCF_000005845.2.fasta",
        "B_subtilis": "./data/task_1/GCF_000009045.1.fasta",
        "P_aeruginosa": "./data/task_1/GCF_000006765.1.fasta",
        "S_aureus": "./data/task_1/GCF_000013425.1.fasta",
        "M_tuberculosis": "./data/task_1/GCF_000195955.2.fasta"
    }

    genomes = load_reference_genomes(fasta_files)

    reads_r1 = load_reads_with_id("./data/task_1/simulated_reads_no_errors_10k_R1.fastq")
    reads_r2 = load_reads_with_id("./data/task_1/simulated_reads_no_errors_10k_R2.fastq")

    automaton_r1 = build_automaton_with_ids(reads_r1)
    matches_r1 = match_reads_by_id(automaton_r1, genomes)

    automaton_r2 = build_automaton_with_ids(reads_r2)
    matches_r2 = match_reads_by_id(automaton_r2, genomes)

    merged_matches = merge_paired_matches(matches_r1, matches_r2)
    counts = count_reads_per_genome(merged_matches)

    for genome in fasta_files:
        print(genome + ": " + str(counts[genome]) + " reads matched")


if __name__ == "__main__":
    task()
