from argparse import ArgumentParser
from Bio import SeqIO
from numpy import ndarray, array, sort, cumsum, searchsorted, sum, max


def compute_n(lengths: ndarray,
              threshold: float):
    sorted_values = sort(lengths)[::-1]
    cumsum_values = cumsum(sorted_values)
    threshold_value = cumsum_values[-1] * threshold
    indices = searchsorted(cumsum_values, threshold_value)

    return sorted_values[indices] if indices < len(sorted_values) else 0


def compute_l(lengths: ndarray,
              threshold: float):
    sorted_values = sort(lengths)[::-1]
    cumsum_values = cumsum(sorted_values)
    threshold_value = cumsum_values[-1] * threshold
    location = searchsorted(cumsum_values, threshold_value)
    return location + 1


def compute_gc_content(sequence: str):
    return (sequence.count("G") + sequence.count("C")) / float(len(sequence)) * 100 if sequence else 0


def evaluate_assembly(fasta_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))

    if not records:
        return {}

    lengths = array([len(record.seq) for record in records])

    total_length = sum(lengths)
    num_contigs = len(records)
    longest = max(lengths)
    gc_content = compute_gc_content("".join(str(rec.seq) for rec in records))

    n50 = compute_n(lengths, 0.5)
    n90 = compute_n(lengths, 0.9)
    l50 = compute_l(lengths, 0.5)

    print("Total Length: %d" % total_length)
    print("Number of Contigs: %d" % num_contigs)
    print("GC Content: %.2f%%" % gc_content)
    print("Longest Contig: %d" % longest)
    print("N50: %d" % n50)
    print("N90: %d" % n90)
    print("L50: %d" % l50)


if __name__ == "__main__":
    parser = ArgumentParser(description="QUAST metrics")
    parser.add_argument("-i", "--input", required=True, help="input FASTA file of assembled contigs")
    args = parser.parse_args()
    evaluate_assembly(args.input)
