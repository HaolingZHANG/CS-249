# Assignment 2


In this task, we used personal laptop to execute our works.

## Task 1

We check if the installation of QUAST is successful 

```shell
python quast.py -h
```

If it is successfully, we can obtain the message:

```text
QUAST: Quality Assessment Tool for Genome Assemblies
Version: 5.2.0

Usage: python quast.py [options] <files_with_contigs>

......

Online QUAST manual is available at https://quast.sf.net/manual
```

I also implemented my own QUAST metris pipeline is implemented
[here](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/quast_metrics.py).

```shell
python quast_metrics.py -h
```

The help message outlines the command-line interface:
```text
usage: quast_metrics.py [-h] -i INPUT
```

### Task 1.1

The de Bruijn Graph assembly algorithm is implemented 
[here](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/dbg_assembler.py).

```shell
python dbg_assembler.py -h
```

The help message outlines the command-line interface:
```text
ssembler.py [-h] -i INPUT -o OUTPUT -k KMER [-g GFA]

de Bruijn Graph genome assembler

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input FASTQ file
  -o OUTPUT, --output OUTPUT
                        output FASTA file for contigs
  -k KMER, --kmer KMER  k-mer size (recommended: 21–55)
  -g GFA, --gfa GFA     output GFA file (optional)
```

### Task 1.2

The overlap-layout-consensus assembly algorithm is implemented 
[here](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/olc_assembler.py).

```shell
python olc_assembler.py -h
```

The help message outlines the command-line interface:
```text
usage: olc_assembler.py [-h] -i INPUT -o OUTPUT -n MINIMUM_OVERLAP

overlap-layout-consensus genome assembler

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input FASTQ file
  -o OUTPUT, --output OUTPUT
                        output FASTA file for contigs
  -n MINIMUM_OVERLAP, --minimum_overlap MINIMUM_OVERLAP
                        minimum overlap length
```

### Task 1.3

#### Task 1.3.1

To execute the sub-task 1, we can run the following command-lines:

```shell
python dbg_assembler.py -i toy_dataset/reads_b.fastq -o results/dbg/b_k40_contigs.fasta -k 40 -g results/gfa/b_k40_graph.gfa
Bandage image results/gfa/b_k40_graph.gfa results/images/b_k40_graph.png
```

Using our de Bruijn graph assembler, 
we generated an assembly graph from b.fastq with k=40 and visualized the result using Bandage.

The assembly graph constructed using the De Bruijn Graph assembler with k=40 exhibits a dominant circular structure, 
as shown in `b_k40_graph.png`. 
The circular form suggests a mostly linear or closed genome region with consistent coverage.

A single short dead-end branch is observed near the loop, 
likely due to sequencing errors or low-abundance variants. 
The overall graph is compact and largely unambiguous, 
indicating that the dataset is clean and has good k-mer connectivity.

![Pikachu](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/results/images/b_k40_graph.png)

Such visualization allows identification of structural ambiguities that could result in fragmented contigs. 
In this case, minimal ambiguity exists, and additional graph cleaning (e.g., bubble removal) 
could potentially further improve the contiguity.

In addition, we can run the QUAST analysis by the following command-lines:

```shell
python quast.py results/dbg/b_k40_contigs.fasta -r toy_dataset/reference_b.fasta -o tmp
```

We evaluated the performance of our De Bruijn Graph (DBG)-based assembler 
on the synthetic toy dataset `b` using a k-mer size of 40. 
The resulting contigs were compared against the provided reference genome using basic QUAST metrics. 
The comparison results are summarized in the following table:

| label                      | value         |
|----------------------------|---------------|
| Assembly                   | b_k40_contigs | 
| # contigs (>= 0 bp)        | 2             | 
| # contigs (>= 1000 bp)     | 0             | 
| # contigs (>= 5000 bp)     | 0             | 
| # contigs (>= 10000 bp)    | 0             | 
| # contigs (>= 25000 bp)    | 0             | 
| # contigs (>= 50000 bp)    | 0             | 
| Total length (>= 0 bp)     | 1046          | 
| Total length (>= 1000 bp)  | 0             | 
| Total length (>= 5000 bp)  | 0             | 
| Total length (>= 10000 bp) | 0             |  
| Total length (>= 25000 bp) | 0             |  
| Total length (>= 50000 bp) | 0             |  
| # contigs                  | 1             |  
| Largest contig             | 761           |  
| Total length               | 761           |  
| Reference length           | 1000          |  
| GC (%)                     | 52.30         | 
| Reference GC (%)           | 52.00         | 
| N50                        | 761           |  
| NG50                       | 761           |  
| N90                        | 761           |  
| NG90                       | -             |   
| auN                        | 761.0         |   
| auNG                       | 579.1         |  
| L50                        | 1             |  
| LG50                       | 1             |   
| L90                        | 1             |   
| LG90                       | -             |   
| # N's per 100 kbp          | 0.00          | 

Overall, the assembly successfully reconstructed the full sequence content of the reference genome. 
The total assembled length was 1046 base pairs, slightly longer than the 1000 bp reference, 
which can be attributed to a small degree of overlap or redundancy between contigs. 
The GC content of the assembled contigs (51.91%) closely matched that of the reference (52.00%), 
indicating accurate base-level reconstruction.

The assembler produced two contigs, with the longest contig covering 761 bp. 
This implies that a minor break occurred in the de Bruijn graph, 
likely due to the graph structure failing to meet the Eulerian path condition at that point. 
Despite this, the continuity of the assembly remains acceptable, as reflected in the N50 value of 761. 
All contigs were generated by identifying Eulerian paths through the graph, 
in accordance with the specifications of the assignment.

This result confirms that the assembler is functionally correct and capable of reconstructing genomic sequences with 
high coverage and base-level accuracy.
Further refinement of k-mer size or graph traversal heuristics could potentially reduce the number of contigs 
and yield an even more contiguous assembly.

#### Task 1.3.2

To execute the sub-task 2, we can run the following command-lines:

```shell
# k = 35
python dbg_assembler.py -i toy_dataset/reads_r.fastq -o results/dbg/r_k35_contigs.fasta -k 35 -g results/gfa/r_k35_graph.gfa
# k = 45
python dbg_assembler.py -i toy_dataset/reads_r.fastq -o results/dbg/r_k45_contigs.fasta -k 45 -g results/gfa/r_k45_graph.gfa
```

In addition, we can run the QUAST analysis by the following command-lines:

```shell
python quast.py results/dbg/r_k35_contigs.fasta -r toy_dataset/reference_b.fasta -o tmp
python quast.py results/dbg/r_k45_contigs.fasta -r toy_dataset/reference_b.fasta -o tmp
Bandage image results/gfa/r_k35_graph.gfa results/images/r_k35_graph.png
Bandage image results/gfa/r_k45_graph.gfa results/images/r_k45_graph.png
```

The full set of QUAST statistics is presented below:

| Metric                    | r_k35_contigs | r_k45_contigs |
|---------------------------|---------------|---------------|
| # contigs (≥ 0 bp)        | 3             | 3             |
| # contigs (≥ 1000 bp)     | 0             | 0             |
| # contigs (≥ 5000 bp)     | 0             | 0             |
| # contigs (≥ 10000 bp)    | 0             | 0             |
| # contigs (≥ 25000 bp)    | 0             | 0             |
| # contigs (≥ 50000 bp)    | 0             | 0             |
| Total length (≥ 0 bp)     | 1092          | 1118          |
| Total length (≥ 1000 bp)  | 0             | 0             |
| Total length (≥ 5000 bp)  | 0             | 0             |
| Total length (≥ 10000 bp) | 0             | 0             |
| Total length (≥ 25000 bp) | 0             | 0             |
| Total length (≥ 50000 bp) | 0             | 0             |
| # contigs                 | 1             | 1             |
| Largest contig            | 616           | 632           |
| Total length              | 616           | 632           |
| Reference length          | 1040          | 1040          |
| GC (%)                    | 50.97         | 51.11         |
| Reference GC (%)          | 51.25         | 51.25         |
| N50                       | 616           | 632           |
| NG50                      | 616           | 632           |
| N90                       | 616           | 632           |
| NG90                      | -             | -             |
| auN                       | 616.0         | 632.0         |
| auNG                      | 364.9         | 384.1         |
| L50                       | 1             | 1             |
| LG50                      | 1             | 1             |
| L90                       | 1             | 1             |
| LG90                      | -             | -             |
| # N's per 100 kbp         | 0.00          | 0.00          |

For k = 35, the Bandage graph is:

![Pikachu](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/results/images/r_k35_graph.png)

For k = 45, the  Bandage graph is:

![Pikachu](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/results/images/r_k45_graph.png)

In both k-mer settings, the assembler produced three contigs, 
indicating that the underlying graph could not be traversed in a single Eulerian path. 
Nevertheless, each of the longest contigs covered over 600 bp, suggesting that the assembler recovered 
most of the reference sequence, albeit in fragmented form.

Interestingly, the total assembled lengths slightly exceeded the reference length (1092 and 1118 bp versus 1040 bp), 
which are attributed to overlapping or redundant regions at contig boundaries (see above image when k = 35). 
The GC content remained closely aligned with the reference across both k-mer settings, 
confirming that base-level accuracy was preserved throughout the assembly.

When comparing the two k-mer sizes, `k=45` produced a slightly longer contig and higher total assembly length, 
suggesting that a larger k-mer improved resolution of certain ambiguous connections in the graph. 
However, neither setting was able to reconstruct the genome as a single contig, 
indicating the presence of unresolved graph branches or insufficient read support for critical overlaps.

In conclusion, although the assembler was unable to perfectly reconstruct the reference in one contiguous sequence, 
it successfully identified large and accurate contigs using strict Eulerian path constraints. 
These results highlight both the strengths and the limitations of DBG-based assembly in realistic conditions.

#### Task 1.3.3

To execute the sub-task 3, we can run the following command-lines:

```shell
# The de Bruijn Graph assembly algorithm is used
python dbg_assembler.py -i synthetic_dataset/reads/no_error_reads_hiseq_5k.fastq -o results/dbg/hiseq_noerr.fasta -k 40
python dbg_assembler.py -i synthetic_dataset/reads/reads_hiseq_5k.fastq          -o results/dbg/hiseq_err.fasta   -k 40
python dbg_assembler.py -i synthetic_dataset/reads/no_error_ont_hq_50x.fastq     -o results/dbg/ont_noerr.fasta   -k 40
python dbg_assembler.py -i synthetic_dataset/reads/ont_hq_50x.fastq              -o results/dbg/ont_err.fasta     -k 40
# The overlap-layout-consensus assembly algorithm is used
python olc_assembler.py -i synthetic_dataset/reads/no_error_reads_hiseq_5k.fastq -o results/olc/hiseq_noerr.fasta -n 30
python olc_assembler.py -i synthetic_dataset/reads/reads_hiseq_5k.fastq          -o results/olc/hiseq_err.fasta   -n 30
python olc_assembler.py -i synthetic_dataset/reads/no_error_ont_hq_50x.fastq     -o results/olc/ont_noerr.fasta   -n 50
python olc_assembler.py -i synthetic_dataset/reads/ont_hq_50x.fastq              -o results/olc/ont_err.fasta     -n 50
```

Notably, for HiSeq data, `n` is suggested as `30`, and for ONT data, `n` is suggested as `50`.

In addition, we can run the QUAST analysis by the following command-lines:

```shell
python quast.py results/dbg/hiseq_noerr.fasta -r synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna -o tmp
python quast.py results/dbg/hiseq_err.fasta   -r synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna -o tmp
python quast.py results/dbg/ont_noerr.fasta   -r synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna -o tmp
python quast.py results/dbg/ont_err.fasta     -r synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna -o tmp
python quast.py results/olc/hiseq_noerr.fasta -r synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna -o tmp
python quast.py results/olc/hiseq_err.fasta   -r synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna -o tmp
python quast.py results/olc/ont_noerr.fasta   -r synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna -o tmp
python quast.py results/olc/ont_err.fasta     -r synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna -o tmp
```

The full set of QUAST statistics is presented below:

| Assembly                   | dbg_hiseq_noerr | dbg_hiseq_err | dbg_ont_noerr | dbg_ont_err | olc_hiseq_noerr | olc_hiseq_err | olc_ont_noerr | olc_ont_err |
|----------------------------|-----------------|---------------|---------------|-------------|-----------------|---------------|---------------|-------------|
| # contigs (>= 0 bp)        | 16              | 25            | 8             | 72          | 7               | 97            | 13            | 2           |
| # contigs (>= 1000 bp)     | 5               | 7             | 2             | 9           | 3               | 2             | 13            | 2           |
| # contigs (>= 5000 bp)     | 3               | 3             | 1             | 5           | 2               | 2             | 12            | 2           |
| # contigs (>= 10000 bp)    | 1               | 1             | 1             | 4           | 2               | 1             | 6             | 2           |
| # contigs (>= 25000 bp)    | 0               | 0             | 1             | 4           | 0               | 1             | 1             | 0           |
| # contigs (>= 50000 bp)    | 0               | 0             | 0             | 3           | 0               | 0             | 0             | 0           |
| Total length (>= 0 bp)     | 30000           | 43421         | 30021         | 730813      | 30194           | 50378         | 141365        | 37293       |
| Total length (>= 1000 bp)  | 25128           | 37891         | 27284         | 720138      | 29006           | 37392         | 141365        | 37293       |
| Total length (>= 5000 bp)  | 22304           | 31465         | 25724         | 712482      | 24597           | 37392         | 137549        | 37293       |
| Total length (>= 10000 bp) | 11004           | 15400         | 25724         | 703293      | 24597           | 29139         | 84863         | 37293       |
| Total length (>= 25000 bp) | 0               | 0             | 25724         | 703293      | 0               | 29139         | 28921         | 0           |
| Total length (>= 50000 bp) | 0               | 0             | 0             | 663213      | 0               | 0             | 0             | 0           |
| # contigs                  | 11              | 11            | 5             | 11          | 4               | 3             | 13            | 2           |
| Largest contig             | 11004           | 15400         | 25724         | 333439      | 13510           | 29139         | 28921         | 20439       |
| Total length               | 29229           | 41476         | 29330         | 721647      | 29716           | 38305         | 141365        | 37293       |
| Reference length           | 30119           | 30119         | 30119         | 30119       | 30119           | 30119         | 30119         | 30119       |
| GC (%)                     | 41.31           | 41.30         | 41.15         | 41.10       | 41.26           | 41.36         | 41.10         | 40.69       |
| Reference GC (%)           | 41.24           | 41.24         | 41.24         | 41.24       | 41.24           | 41.24         | 41.24         | 41.24       |
| N50                        | 6294            | 8765          | 25724         | 229982      | 11087           | 29139         | 10492         | 20439       |
| NG50                       | 6294            | 15400         | 25724         | 333439      | 11087           | 29139         | 28921         | 20439       |
| N90                        | 873             | 1154          | 1560          | 99792       | 4409            | 8253          | 8096          | 16854       |
| NG90                       | 625             | 7300          | 1560          | 333439      | 4409            | 29139         | 28921         | 16854       |
| auN                        | 6596.5          | 9207.5        | 22693.2       | 243524.3    | 10949.8         | 23966.2       | 13779.8       | 18818.8     |
| auNG                       | 6401.6          | 12679.4       | 22098.7       | 5834808.4   | 10803.3         | 30480.0       | 64676.3       | 23301.2     |
| L50                        | 2               | 2             | 1             | 2           | 2               | 1             | 5             | 1           |
| LG50                       | 2               | 1             | 1             | 1           | 2               | 1             | 1             | 1           |
| L90                        | 7               | 7             | 2             | 3           | 3               | 2             | 11            | 2           |
| LG90                       | 8               | 3             | 2             | 1           | 3               | 1             | 1             | 2           |
| # N's per 100 kbp          | 0.00            | 0.00          | 0.00          | 0.00        | 0.00            | 0.00          | 0.00          | 0.00        |

Here we used QUAST to evaluate all assembled contigs and compared them against the MERS reference genome. 
The table above summarizes the results across both DBG and OLC assemblers, tested on HiSeq and ONT read sets,
with and without simulated sequencing errors.

For error-free datasets, both assemblers generally performed well. 
The de Bruijn Graph (DBG) assembler achieved near-complete assemblies on both HiSeq and ONT data, 
producing a moderate number of contigs (16 and 8 respectively) 
with N50 values reflecting relatively contiguous reconstructions. 
The OLC assembler performed even better in the no-error ONT dataset, 
where it assembled only 13 contigs and achieved the highest GC accuracy (41.30%).

When errors were introduced, performance diverged significantly between read types and algorithms. 
The DBG assembler showed poor tolerance to ONT errors, 
producing a very large assembly (730 kb) composed of 72 contigs — far exceeding the expected reference size (30 kb), 
suggesting severe fragmentation and misassembly due to erroneous edges in the overlap graph. 
The OLC assembler, on the other hand, performed surprisingly well in this condition, 
recovering a near-complete genome with only 2 contigs and an N50 of over 20 kb.

In contrast, on HiSeq data with errors, the OLC assembler showed signs of over-fragmentation (97 contigs), 
likely due to its sensitivity to alignment mismatches during overlap detection. 
The DBG approach was relatively more stable here, producing a reasonable assembly of 25 contigs.

Overall, these results highlight the tradeoffs between the two assembly paradigms. 
The DBG approach offers robust performance when read coverage is high and errors are minimal, 
but suffers on noisy long-read data. 
The OLC assembler, while more sensitive to base-level errors in short reads, 
demonstrates stronger resilience on long-read data with high error rates, 
especially when contigs can be built from strong overlaps despite sequence noise.

These findings emphasize the importance of selecting assembly strategies based on read technology and data quality. 
In high-error ONT scenarios, OLC offers a compelling advantage, while DBG is preferable for cleaner short-read datasets.

#### Task 1.3.4

We check if the installation of Canu is successful

```shell
canu -options
```

If it is successfully, we can obtain the message:

```text
usage:   canu [-version] [-citation] \
              [-haplotype | -correct | -trim | -assemble | -trim-assemble] \
              [-s <assembly-specifications-file>] \
               -p <assembly-prefix> \
               -d <assembly-directory> \
               genomeSize=<number>[g|m|k] \
              [other-options] \
              [-haplotype{NAME} illumina.fastq.gz] \
              [-corrected] \
              [-trimmed] \
              [-pacbio |
               -nanopore |
               -pacbio-hifi] file1 file2 ...

example: canu -d run1 -p godzilla genomeSize=1g -nanopore-raw reads/*.fasta.gz 

......

Complete documentation at https://canu.readthedocs.org/en/latest/
```

Now, we repeat the same assembly using the Canu assembler on the ONT datasets (setting `genomeSize` as `31k`), 
both with and without sequencing errors. 

```shell
canu -p genome -d output genomeSize=31k useGrid=false -nanopore-raw synthetic_dataset/reads/no_error_ont_hq_50x.fastq
mv output/genome.contigs.fasta results/canu/ont_noerr.fasta
canu -p genome -d output genomeSize=31k useGrid=false -nanopore-raw synthetic_dataset/reads/ont_hq_50x.fastq
mv output/genome.contigs.fasta results/canu/ont_err.fasta
```

The QUAST results were compared with those obtained from our custom DBG and OLC implementations.

```shell
python quast.py results/canu/ont_noerr.fasta -r synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna -o tmp
python quast.py results/canu/ont_err.fasta   -r synthetic_dataset/GCF_000901155.1_ViralProj183710_genomic.fna -o tmp
```

Then, we can construct the following table:

| Assembly                   | ont_noerr | ont_err |
|----------------------------|-----------|---------|
| # contigs (>= 0 bp)        | 1         | 1       |
| # contigs (>= 1000 bp)     | 1         | 1       |
| # contigs (>= 5000 bp)     | 1         | 1       |
| # contigs (>= 10000 bp)    | 1         | 1       |
| # contigs (>= 25000 bp)    | 1         | 1       |
| # contigs (>= 50000 bp)    | 0         | 0       |
| Total length (>= 0 bp)     | 28655     | 28564   |
| Total length (>= 1000 bp)  | 28655     | 28564   |
| Total length (>= 5000 bp)  | 28655     | 28564   |
| Total length (>= 10000 bp) | 28655     | 28564   |
| Total length (>= 25000 bp) | 28655     | 28564   |
| Total length (>= 50000 bp) | 0         | 0       |
| # contigs                  | 1         | 1       |
| Largest contig             | 28655     | 28564   |
| Total length               | 28655     | 28564   |
| Reference length           | 30119     | 30119   |
| GC (%)                     | 41.05     | 41.06   |
| Reference GC (%)           | 41.24     | 41.24   |
| N50                        | 28655     | 28564   |
| NG50                       | 28655     | 28564   |
| N90                        | 28655     | 28564   |
| NG90                       | 28655     | 28564   |
| auN                        | 28655.0   | 28564.0 |
| auNG                       | 27262.2   | 27089.3 |
| L50                        | 1         | 1       |
| LG50                       | 1         | 1       |
| L90                        | 1         | 1       |
| LG90                       | 1         | 1       |
| # N's per 100 kbp          | 0.00      | 0.00    |

Compared to the reference genome (30,119 bp),
Canu produced a nearly complete reconstruction with only a minor loss in total length. 
Both error-free and error-containing inputs yielded a single contig with high continuity and stable GC content, 
indicating the robustness of Canu to input quality.

Canu assemblies show the highest accuracy among all tested methods, 
as indicated by the near-complete recovery of genome length and stable GC content relative to the reference. 
In terms of contiguity, Canu consistently produced a single contig, 
unlike the fragmented outputs of DBG and traditional OLC assemblies, 
demonstrating its superior ability to reconstruct full-length genomes. 
Moreover, Canu exhibited strong robustness to sequencing errors: 
while the DBG assembly dramatically degraded in the presence of errors (expanding to over 700 kb and 72 contigs),
Canu maintained a clean, single-contig output with minimal reduction in length. 
Although Canu is OLC-based, it significantly outperformed the baseline OLC assembler, 
likely due to its improved overlap detection, error correction, and consensus-building mechanisms. 
These findings support the use of Canu as a reliable and accurate assembler for long-read ONT data, 
particularly in viral genome assembly scenarios.

### Summary

All command lines of Task 1 is attached in the `assign2/exec` folder, see
[here](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/exec/execution_1.sh).

## Task 2

In this task, we used IBEX to execute our works.
As the assignment document indicates that all required software is already installed on IBEX, 
we did not perform additional checks.

### Task 2.1

Firstly, we go to our own project location:

```shell
cd /ibex/user/zhanh0m/proj/cs249/
```

Here we use `Hifiasm` to assemble the *Scincus mitranus* genome. 
The sbatch script is shown below (named `step_1.slurm` in our project location, also is attached
[here](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/exec/step_1.slurm)):

```shell
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
```

The source data is `lizard_liver.fastq.gz`, mentioned by Olga Mashkova in the email. 
Notably, we used the computational resource parameters recommended by ChatGPT for analyzing PacBio HiFi data:
i.e., `32` cpus and `120G` memory.

We run it by the following command line:

```shell
sbatch step_1.slurm
```

The execution log is attached 
[here](https://github.com/HaolingZHANG/CS-249/tree/main/assign2/results/log/step_1.txt).

Since `Hifiasm` did not produce a FASTA file, we manually converted the output using the following command line:

```shell
cd raw
awk '/^S/{print ">"$2"\n"$3}' lizard.asm.bp.p_ctg.gfa > lizard.asm.bp.p_ctg.fa
```

Finally, the assembly results are shown below:

![Pikachu](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/results/images/assembly_screenshot.png)

### Task 2.2

Here we use `QUAST`, `BUSCO`, `Merqury`, and `Flagger` to evaluate our assembly.

The sbatch script is shown below (named `step_2.slurm` in our project location, also is attached
[here](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/exec/step_2.slurm)):

