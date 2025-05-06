from argparse import ArgumentParser
from collections import defaultdict
from copy import deepcopy
from Bio import Seq, SeqRecord, SeqIO
from typing import List


class OLCAssembler:
    def __init__(self,
                 reads: list,
                 minimum_overlap: int):
        self.reads = reads
        self.reads_dict = {read.id: read for read in reads}
        self.minimum_overlap = minimum_overlap
        self.graph = defaultdict(list)
        self.layout_order = []

    def compute_overlap(self,
                        sequence_1: str,
                        sequence_2: str):
        start = 0
        while True:
            start = sequence_1.find(sequence_2[:self.minimum_overlap], start)
            if start == -1:
                return 0
            if sequence_2.startswith(sequence_1[start:]):
                return len(sequence_1) - start
            start += 1

    def build_overlap_graph(self):
        for index_1, sequence_1 in enumerate(self.reads):
            for index_2, sequence_2 in enumerate(self.reads):
                if index_1 != index_2:
                    overlap_length = self.compute_overlap(str(sequence_1.seq), str(sequence_2.seq))
                    if overlap_length > 0:
                        self.graph[sequence_1.id].append((sequence_2.id, overlap_length))

    def build_layouts(self):
        visited, layouts, graph_backup = set(), [], deepcopy(self.graph)
        for start_node in list(self.graph):
            if start_node in visited:
                continue

            layout_order = []
            stack = [start_node]
            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    layout_order.append(node)
                    neighbors = sorted(self.graph[node], key=lambda x: -x[1])
                    for neighbor, _ in neighbors:
                        if neighbor not in visited:
                            stack.append(neighbor)

            if layout_order:
                layouts.append(layout_order)

        self.layout_order, self.graph = layouts, graph_backup

    def generate_consensus(self):
        contigs = []
        for layout in self.layout_order:

            if not layout:
                continue

            consensus = str(self.reads_dict[layout[0]].seq)
            for i in range(1, len(layout)):
                prev, curr = layout[i - 1], layout[i]
                for neighbor, olen in self.graph[prev]:
                    if neighbor == curr:
                        consensus += str(self.reads_dict[curr].seq)[olen:]
                        break

            contigs.append(consensus)

        return contigs

    def __call__(self, *args, **kwargs):
        self.build_overlap_graph()
        self.build_layouts()
        return self.generate_consensus()


def read_fastq(fastq_path: str) \
        -> List[SeqRecord.SeqRecord]:
    return list(SeqIO.parse(fastq_path, "fastq"))


def write_fasta(contigs: list,
                output_path: str):
    records = []
    for index, sequence in enumerate(contigs):
        record = SeqRecord.SeqRecord(Seq.Seq(sequence),
                                     id="contig_%d" % (index + 1),
                                     description="length=%d" % len(sequence))
        records.append(record)
    SeqIO.write(records, output_path, "fasta")


def main():
    parser = ArgumentParser(description="overlap-layout-consensus genome assembler")
    parser.add_argument("-i", "--input", required=True, help="input FASTQ file")
    parser.add_argument("-o", "--output", required=True, help="output FASTA file for contigs")
    parser.add_argument("-n", "--minimum_overlap", type=int, required=True, help="minimum overlap length")
    args = parser.parse_args()

    reads = read_fastq(args.input)
    assembler = OLCAssembler(reads, minimum_overlap=args.minimum_overlap)
    contigs = assembler()

    write_fasta(contigs, args.output)
    print("[INFO] Assembly complete. %d contig(s) written to %s." % (len(contigs), args.output))


if __name__ == "__main__":
    main()
