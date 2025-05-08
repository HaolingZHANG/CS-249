#!/bin/bash
#SBATCH --job-name=194913_lizard
#SBATCH --account=cs249
#SBATCH --output=run.out
#SBATCH --error=run.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --partition=batch

module load hifiasm

in_file="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver.fastq.gz"
out_file="/ibex/user/zhanh0m/proj/cs249/raw/lizard.asm"

hifiasm -o "$out_file" -t 32 "$in_file"