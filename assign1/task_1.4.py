from Bio import SeqIO
from os import makedirs, system, remove
from time import perf_counter
from memory_profiler import profile


def run_blast(query_fastq, ref_fasta, genome_name, output_dir="./temp/"):
    makedirs(output_dir, exist_ok=True)

    db_prefix = output_dir + genome_name + "db"
    output_file = "./data/task_1/" + genome_name + ".tsv"

    system("makeblastdb -in " + ref_fasta + " -dbtype nucl -out " + db_prefix + " -logfile NUL")

    blast_command = (
        "blastn -query " + query_fastq + " -db " + db_prefix +
        " -outfmt \"6 qseqid sseqid mismatch length evalue bitscor\" "
        " -out " + output_file
    )
    system(blast_command)

    remove("./temp/" + genome_name + "db.ndb")
    remove("./temp/" + genome_name + "db.nhr")
    remove("./temp/" + genome_name + "db.nin")
    remove("./temp/" + genome_name + "db.njs")
    remove("./temp/" + genome_name + "db.not")
    remove("./temp/" + genome_name + "db.nsq")
    remove("./temp/" + genome_name + "db.ntf")
    remove("./temp/" + genome_name + "db.nto")


def merge_pe_read_ids(file1: str,
                      file2: str,
                      mismatch_threshold: int = 1):
    matched = set()

    def process_file(tsv: str):
        with open(tsv) as file:
            for line in file:
                cols = line.strip().split('\t')
                raw_id = cols[0]
                mismatches = int(cols[2])
                if mismatches <= mismatch_threshold:
                    clean_id = raw_id.split('/')[0]  # remove /1 or /2
                    matched.add(clean_id)

    process_file(file1)
    process_file(file2)

    return matched


@profile
def unit_task(species: str,
              fasta_file: str,
              fastq_file_pair: list):
    start = perf_counter()
    run_blast(fastq_file_pair[0], fasta_file, species + "_R1")
    run_blast(fastq_file_pair[1], fasta_file, species + "_R2")

    matched = merge_pe_read_ids("./data/task_1/" + species + "_R1.tsv", "./data/task_1/" + species + "_R2.tsv")

    end = perf_counter()

    return len(matched), end - start


def task():
    fasta_files = {
        "E_coli": "./data/task_1/GCF_000005845.2.fasta",
        "B_subtilis": "./data/task_1/GCF_000009045.1.fasta",
        "P_aeruginosa": "./data/task_1/GCF_000006765.1.fasta",
        "S_aureus": "./data/task_1/GCF_000013425.1.fasta",
        "M_tuberculosis": "./data/task_1/GCF_000195955.2.fasta"
    }

    pair = [
        "./data/task_1/simulated_reads_miseq_10k_R1.fastq",
        "./data/task_1/simulated_reads_miseq_10k_R2.fastq"
    ]

    SeqIO.convert(pair[0], "fastq", "./temp/R1.fasta", "fasta")
    SeqIO.convert(pair[1], "fastq", "./temp/R2.fasta", "fasta")

    for species, fasta_file in fasta_files.items():
        print(species, unit_task(species, fasta_file, ["./temp/R1.fasta", "./temp/R1.fasta"]))

    remove("./temp/R1.fasta")
    remove("./temp/R2.fasta")


if __name__ == "__main__":
    task()

