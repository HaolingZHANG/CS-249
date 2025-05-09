# Assignment 2

Author: Haoling Zhang (KAUST ID = 194913)

## Repository structure

```html
└── assign2                           // Source codes of Assignment 2.
      ├── exec                        // All executable scripts and batch executable scripts.
      ├── results                     // All result files, including organized logs during execution.
      ├── toy_dataset                 // Dataset provided by the instructor.
      ├── synthetic_dataset           // Dataset provided by the instructor.
      ├── tmp                         // Temp folder, used to store some data that needs to be organized.
      ├── dbg_assembler.py            // de Bruijn Graph (DBG) Assembly tool (for Task 1.1).
      ├── olc_assembler.py            // Overlap-Layout-Consensus (OLC) Assembly tool (for Task 1.2).
      ├── quast_metrics.py            // QUAST tool implemented by the student.
      ├── merge_quast.py              // The process of merging multiple QUAST reports to Markdown table.
      ├── distribution_plot.py        // The process of plotting the k-mer distribution results.
      └── REPORT.md                   // Markdown-type report.
```

So, you can also check the report 
[online](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/REPORT.md)

## Execution usages

If you wish to reproduce all the results of Task 1, 
please navigate to the assign2 directory by
```shell
cd assign2
```
and execute the command lines in `exec/executioner_1.sh` one by one in the console.

If you wish to reproduce all the results of Task 2.1 and 2.2, 
please navigate to my cs249 directory in IBEX by
```shell
cd /ibex/user/zhanh0m/proj/cs249/
```
copy `step_1.slurm` (see 
[here](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/exec/step_1.slurm)) and `step_2.slurm` (see 
[here](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/exec/step_2.slurm)) in this directory,
and execute the command lines in `exec/executioner_2.sh` (see 
[here](https://github.com/HaolingZHANG/CS-249/blob/main/assign2/exec/execution_2.sh)) 
one by one in the console.

## Algorithm usages

If you want to use the assembly algorithm directly, please download the repository and 
```shell
cd assign2
```

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

## Files

All assembly and report files of Task 1 are attached in `assign2/results/` (see 
[here](https://github.com/HaolingZHANG/CS-249/tree/main/assign2/results))

The report files of Task 2.1 and 2.2 are attached in `assign2/results/hifiasm` (see
[here](https://github.com/HaolingZHANG/CS-249/tree/main/assign2/results/hifiasm)).

The assembly and detailed report files of Task 2.1 and 2.2 are attached in `/ibex/user/zhanh0m/proj/cs249/`.

The assembly and detailed report files of Task 2.3 are attached in `/ibex/user/zhanh0m/proj/cs249a/`.

## LLM usage declaration

I used ChatGPT (GPT-4o) in following specific situations:

- ChatGPT recommended Computational resource parameters for all the computational tools (for Task 2 only). 
Given the execution time of each task, I'm unable to incrementally debug from low to high resource allocations.
- Debugging and usage of `Inspector`. This tool cannot be activated correctly through `module load inspector` 
in IBEX (cannot be directly provided as an executable command).