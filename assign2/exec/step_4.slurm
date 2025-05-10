#!/bin/bash
#SBATCH --job-name=194913_lizard
#SBATCH --account=cs249
#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --time=48:00:00
#SBATCH --partition=batch

module load merqury
module load inspector

fasta_file="/ibex/user/zhanh0m/proj/cs249a/raw/lizard.asm.bp.p_ctg.fa"
reads_file="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver.fastq.gz"

project_folder="/ibex/user/zhanh0m/proj/cs249a/"
merqury_database="${project_folder}raw/reads.meryl"
inspector_output="${project_folder}inspector/"

if [ ! -d "$merqury_database" ]; then
    meryl k=21 count output "$merqury_database" "$reads_file"
fi
merqury.sh "$merqury_database" "$fasta_file" "merqury/lizard"

python /ibex/sw/rl9c/inspector/1.3/rl9.4_conda3/Miniconda3/envs/python2.7/bin/inspector.py \
  -c "${project_folder}raw/lizard.asm.bp.p_ctg.fa" -r "$reads_file" -d hifi -t 32 -o "$inspector_output"
