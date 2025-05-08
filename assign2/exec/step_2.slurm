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

module load quast
module load busco
module load merqury
module load minimap2
module load samtools
module load flagger

fasta_file="/ibex/user/zhanh0m/proj/cs249/raw/lizard.asm.bp.p_ctg.fa"
reads_file="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver.fastq.gz"

project_folder="/ibex/user/zhanh0m/proj/cs249/"
quast_output="${project_folder}quast/"
merqury_output="${project_folder}merqury/"
merqury_database="${project_folder}raw/reads.meryl"
flagger_output="${project_folder}flagger/"
samtool_output="${project_folder}raw/aln.sorted.bam"

# Step 1: QUAST
quast.py "$fasta_file" -o "$quast_output" -t 32

# Step 2: BUSCO
export JAVA_TOOL_OPTIONS="-Xmx64g"
export BBMAP_JAVA_MEM="-Xmx64g"
busco -i "$fasta_file" -o busco -l "sauropsida_odb10" -m genome -c 32 --out_path "$project_folder" -f --no_stats

# Step 3: Merqury
if [ ! -d "$merqury_database" ]; then
    meryl k=21 count output "$merqury_database" "$reads_file"
fi
merqury.sh "$merqury_database" "$fasta_file" "$merqury_output"

# Step 4: Flagger
minimap2 -t 32 -ax map-hifi "$fasta_file" "$reads_file" | samtools sort -@ 16 -o "$samtool_output"
samtools index "$samtool_output"
flagger -r "$fasta_file" -b "$samtool_output" -o "$flagger_output"