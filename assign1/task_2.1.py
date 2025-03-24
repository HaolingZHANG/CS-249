from Bio import SeqIO
from collections import defaultdict


def extract_kmers(sequence, k):
    for index in range(len(sequence) - k + 1):
        yield sequence[index: index + k]


def build_kmer_index(fasta_files, k=31):
    kmer_index = defaultdict(lambda: defaultdict(int))
    for genome_name, path in fasta_files.items():
        for record in SeqIO.parse(path, "fasta"):
            seq = str(record.seq).upper()
            for kmer in extract_kmers(seq, k):
                kmer_index[kmer][genome_name] += 1
    return kmer_index


def task():
    fasta_files = {
        "E_coli": "./data/task_2/GCF_000005845.2.fasta",
        "B_subtilis": "./data/task_2/GCF_000009045.1.fasta",
        "P_aeruginosa": "./data/task_2/GCF_000006765.1.fasta",
        "S_aureus": "./data/task_2/GCF_000013425.1.fasta",
        "M_tuberculosis": "./data/task_2/GCF_000195955.2.fasta"
    }
    k = 31

    kmer_index = build_kmer_index(fasta_files, k)

    total_kmers = len(kmer_index)
    theoretical_kmers = 4 ** k

    print("Used k = " + str(k))
    print("Total unique " + str(k) + "-mers in index: " + str(total_kmers))
    print("Theoretical number of " + str(k) + "-mers: " + ("%.2e" % theoretical_kmers))
    print("Coverage ratio: " + ("%.2e" % (total_kmers / theoretical_kmers)))


if __name__ == "__main__":
    task()
