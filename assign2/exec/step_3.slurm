#!/bin/bash
#SBATCH --job-name=194913_lizard_clean
#SBATCH --account=cs249
#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --time=48:00:00
#SBATCH --partition=batch

module load fastp
module load kraken2
module load seqtk
module load hifiasm

reads_file="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver.fastq.gz"
source_path="/ibex/user/zhanh0m/proj/cs249a/"
raw_folder="${source_path}raw/"
clean_1_folder="${source_path}clean_1/"
clean_2_folder="${source_path}clean_2/"
out_file="${source_path}raw/lizard.asm"

zcat "$reads_file" > "${raw_folder}reads.fastq"

fastp -i "${raw_folder}reads.fastq" -o "${clean_1_folder}reads.fastq" \
      --qualified_quality_phred 20 --length_required 1000 --thread 32 \
      --html "${clean_1_folder}fastp_report.html" \
      --json "${clean_1_folder}qc_report/fastp_report.json "\
      --report_title "HiFi Read QC (Liver)"

cd "$clean_2_folder"
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz
tar -xzf k2_standard_08gb_20230605.tar.gz
rm k2_standard_08gb_20230605.tar.gz

cd "$source_path"
kraken2 --db "$clean_2_folder" --threads 32 \
        --output "${clean_2_folder}kraken2.out" \
        --report "${clean_2_folder}kraken2.report" \
        "${clean_1_folder}reads.fastq"

awk '($1 == "U") || ($1 == "C" && $3 != 2 && $3 != 10239 && $3 != 9606 && $3 != 4751)'  "${clean_2_folder}kraken2.out" | cut -f2 > "${clean_2_folder}keep_clean_ids.txt"
seqtk subseq "${clean_1_folder}reads.fastq" "${clean_2_folder}keep_clean_ids.txt" > "${clean_2_folder}reads.fastq"

hifiasm -o "$out_file" -t 32 "${clean_2_folder}reads.fastq"
