# Pre-validation
python quast_metrics.py -h
python dbg_assembler.py -h
python olc_assembler.py -h

# Task 1.3.1
python dbg_assembler.py -i toy_dataset/reads_b.fastq -o results/dbg/b_k40_contigs.fasta -k 40 -g results/gfa/b_k40_graph.gfa
Bandage image results/gfa/b_k40_graph.gfa results/images/b_k40_graph.png
python quast_metrics.py -i results/dbg/b_k40_contigs.fasta
python quast_metrics.py -i toy_dataset/reference_b.fasta

# Task 1.3.2
python dbg_assembler.py -i toy_dataset/reads_r.fastq -o results/dbg/r_k35_contigs.fasta -k 35 -g results/gfa/r_k35_graph.gfa
python quast_metrics.py -i results/dbg/r_k35_contigs.fasta
Bandage image results/gfa/r_k35_graph.gfa results/images/r_k35_graph.png
python dbg_assembler.py -i toy_dataset/reads_r.fastq -o results/dbg/r_k45_contigs.fasta -k 45 -g results/gfa/r_k45_graph.gfa
python quast_metrics.py -i results/dbg/r_k45_contigs.fasta
Bandage image results/gfa/r_k45_graph.gfa results/images/r_k45_graph.png

# Task 1.3.3
python dbg_assembler.py -i synthetic_dataset/reads/no_error_reads_hiseq_5k.fastq -o results/dbg/hiseq_noerr.fasta -k 40
python dbg_assembler.py -i synthetic_dataset/reads/reads_hiseq_5k.fastq          -o results/dbg/hiseq_err.fasta   -k 40
python dbg_assembler.py -i synthetic_dataset/reads/no_error_ont_hq_50x.fastq     -o results/dbg/ont_noerr.fasta   -k 40
python dbg_assembler.py -i synthetic_dataset/reads/ont_hq_50x.fastq              -o results/dbg/ont_err.fasta     -k 40
python olc_assembler.py -i synthetic_dataset/reads/no_error_reads_hiseq_5k.fastq -o results/olc/hiseq_noerr.fasta -n 30
python olc_assembler.py -i synthetic_dataset/reads/reads_hiseq_5k.fastq          -o results/olc/hiseq_err.fasta   -n 30
python olc_assembler.py -i synthetic_dataset/reads/no_error_ont_hq_50x.fastq     -o results/olc/ont_noerr.fasta   -n 50
python olc_assembler.py -i synthetic_dataset/reads/ont_hq_50x.fastq              -o results/olc/ont_err.fasta     -n 50
python quast_metrics.py -i results/dbg/hiseq_noerr.fasta
python quast_metrics.py -i results/dbg/hiseq_err.fasta
python quast_metrics.py -i results/dbg/ont_noerr.fasta
python quast_metrics.py -i results/dbg/ont_err.fasta
python quast_metrics.py -i results/olc/hiseq_noerr.fasta
python quast_metrics.py -i results/olc/hiseq_err.fasta
python quast_metrics.py -i results/olc/ont_noerr.fasta
python quast_metrics.py -i results/olc/ont_err.fasta
python quast_metrics.py -i synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna

# Task 1.3.4
canu -p genome -d output genomeSize=31k useGrid=false -nanopore-raw synthetic_dataset/reads/no_error_ont_hq_50x.fastq
mv output/genome.contigs.fasta results/canu/ont_noerr.fasta
canu -p genome -d output genomeSize=31k useGrid=false -nanopore-raw synthetic_dataset/reads/ont_hq_50x.fastq
mv output/genome.contigs.fasta results/canu/ont_err.fasta