# move to project location
cd /ibex/user/zhanh0m/proj/cs249/

# run step 1
sbatch step_1.slurm
# convert the GFA file to FASTA file
awk '/^S/{print ">"$2"\n"$3}' lizard.asm.bp.p_ctg.gfa > lizard.asm.bp.p_ctg.fa

# run step 2
sbatch step_2.slurm