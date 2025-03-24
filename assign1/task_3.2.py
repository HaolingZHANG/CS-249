from os import system, makedirs, path, remove


def build_database(database_name: str):
    if not path.exists(database_name):
        makedirs(database_name)
        system("curl -o ./temp/ https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20241228.tar.gz")
        system("tar -xvzf ./temp/k2_standard_08gb_20241228.tar.gz -C " + database_name)
        remove("./temp/k2_standard_08gb_20241228.tar.gz")


def download_samples(accessions: list):
    for accession in accessions:
        if not path.exists("./temp/" + accession + "_1.fastq"):
            system("prefetch " + accession + " --output-directory ./temp/")
            system("fasterq-dump ./temp/" + accession + "/" + accession + ".sra" + " --split-files --outdir ./temp/")
            remove("./temp/" + accession + "/" + accession + ".sra")


def run_classification(database_name: str,
                       paired_ends: list,
                       output_path: str,
                       report_path: str):
    system("kraken2 --db " + database_name + " --threads 1 " +
           "--report " + report_path + " --output " + output_path +
           " --paired " + " ".join(paired_ends))


def task():
    accessions = ["SRR11412973", "SRR11412976", "SRR11412979", "SRR11412980", "SRR11412984",
                  "SRR21907296", "SRR21907303", "SRR21907307", "SRR21907332", "SRR21907330"]
    build_database(database_name="./data/task_3/database_2")
    download_samples(accessions=accessions)
    for accession in accessions:
        if not path.exists("./data/task_3/o_" + accession + ".txt"):
            print("run " + accession + ":")
            run_classification(database_name="./data/task_3/database_2",
                               paired_ends=["./temp/" + accession + "_1.fastq", "./temp/" + accession + "_2.fastq"],
                               output_path="./data/task_3/o_" + accession + ".txt",
                               report_path="./data/task_3/r_" + accession + ".txt")
            print("Please check the result in \"./data/task_3/r_" + accession + ".\" file.")
            remove("./data/task_3/o_" + accession + ".txt")


if __name__ == "__main__":
    task()
