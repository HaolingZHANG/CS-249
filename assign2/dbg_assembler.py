from argparse import ArgumentParser
from Bio import Seq, SeqRecord, SeqIO
from collections import defaultdict
from copy import deepcopy
from typing import List


class DeBruijnGraph:

    def __init__(self,
                 k: int,
                 reads: list):
        self.k, self.graph, self.in_degree, self.out_degree = k, defaultdict(list), defaultdict(int), defaultdict(int)
        self.reads = reads

    def add_kmer(self,
                 kmer: str):
        prefix, suffix = kmer[:-1], kmer[1:]

        self.graph[prefix].append(suffix)
        self.out_degree[prefix] += 1
        self.in_degree[suffix] += 1

    def build_graph(self):
        for read in self.reads:
            seq = str(read.seq)
            for location in range(len(seq) - self.k + 1):
                kmer = seq[location:location + self.k]
                self.add_kmer(kmer)

    def assemble_contigs(self):
        contigs, visited_edges, graph_backup = [], set(), deepcopy(self.graph)

        for node in list(self.graph.keys()):
            while self.graph[node]:
                path = self.find_eulerian_path(node, visited_edges)
                if path and len(path) >= self.k:
                    contig = path[0]
                    for next_node in path[1:]:
                        contig += next_node[-1]
                    contigs.append(contig)

        self.graph = graph_backup

        return contigs

    def find_eulerian_path(self, start_node, visited_edges):
        path = []
        stack = [start_node]

        while stack:
            current = stack[-1]
            if self.graph[current]:
                neighbor = self.graph[current][-1]
                edge_id = (current, neighbor)
                if edge_id not in visited_edges:
                    visited_edges.add(edge_id)
                    self.graph[current].pop()
                    stack.append(neighbor)
                else:
                    self.graph[current].pop()
            else:
                path.append(stack.pop())

        path.reverse()
        return path

    def export_gfa(self,
                   output_path: str):
        with open(output_path, "w") as file:
            for node in self.graph:
                file.write("S\t%s\t*\n" % node)
            for node in self.graph:
                for neighbor in self.graph[node]:
                    file.write("L\t%s\t+\t%s\t+\t0M\n" % (node, neighbor))

    def __call__(self, *args, **kwargs) \
            -> list:
        self.build_graph()
        return self.assemble_contigs()


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
    parser = ArgumentParser(description="de Bruijn Graph genome assembler")
    parser.add_argument("-i", "--input", required=True, help="input FASTQ file")
    parser.add_argument("-o", "--output", required=True, help="output FASTA file for contigs")
    parser.add_argument("-k", "--kmer", type=int, required=True, help="k-mer size (recommended: 21–55)")
    parser.add_argument("-g", "--gfa", required=False, help="output GFA file (optional)")
    args = parser.parse_args()

    reads = read_fastq(args.input)

    dbg = DeBruijnGraph(args.kmer, reads)
    contigs = dbg()

    write_fasta(contigs, args.output)

    if args.gfa:
        dbg.export_gfa(args.gfa)

    print("[INFO] Assembly complete. %d contig(s) written to %s." % (len(contigs), args.output))


if __name__ == "__main__":
    main()
