# move to project location
cd /ibex/user/zhanh0m/proj/cs249a/

# run step 3
sbatch step_3.slurm

# convert the GFA file to FASTA file
awk '/^S/{print ">"$2"\n"$3}' raw/lizard.asm.bp.p_ctg.gfa > raw/lizard.asm.bp.p_ctg.fa

# run step 4
sbatch step_4.slurm