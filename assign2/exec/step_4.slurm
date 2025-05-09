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

fasta_file="/ibex/user/zhanh0m/proj/cs249a/raw/lizard.asm.bp.p_ctg.fa"

merqury_database="/ibex/user/zhanh0m/proj/cs249/raw/reads.meryl"

merqury.sh "$merqury_database" "$fasta_file" "merqury/lizard"