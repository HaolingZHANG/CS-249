from os import system, makedirs, remove


def build_database(database_name: str,
                   fasta_files: dict):
    makedirs(database_name + "/taxonomy", exist_ok=True)
    system("curl -o " + database_name + "/taxdump.tar.gz https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz")
    system("tar -xzvf " + database_name + "/taxdump.tar.gz -C " + database_name + "/taxonomy")
    for name, fasta_file in fasta_files.items():
        system("kraken2-build --add-to-library " + fasta_file + " --db " + database_name)
    system("kraken2-build --build --db " + database_name)


def run_classification(database_name: str,
                       paired_ends: list,
                       output_path: str,
                       report_path: str,):
    system("kraken2 --db " + database_name + " --threads 1 " +
           "--report " + report_path + " --output " + output_path +
           " --paired " + " ".join(paired_ends))


def task():
    fasta_files = {
        "E_coli": "./data/task_3/GCF_000005845.2.fasta",
        "B_subtilis": "./data/task_3/GCF_000009045.1.fasta",
        "P_aeruginosa": "./data/task_3/GCF_000006765.1.fasta",
        "S_aureus": "./data/task_3/GCF_000013425.1.fasta",
        "M_tuberculosis": "./data/task_3/GCF_000195955.2.fasta"
    }
    build_database(database_name="./data/task_3/database_1", fasta_files=fasta_files)
    run_classification(database_name="./data/task_3/database_1",
                       paired_ends=["./data/task_3/simulated_reads_no_errors_10k_R1.fastq",
                                    "./data/task_3/simulated_reads_no_errors_10k_R2.fastq"],
                       output_path="./data/task_3/output.txt", report_path="./data/task_3/report.txt")
    print("Please check the result in \"./data/task_3/report_1.txt\" file.")
    remove("./data/task_3/output.txt")


if __name__ == "__main__":
    task()
