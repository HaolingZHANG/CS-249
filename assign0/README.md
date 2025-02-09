# Assignment 0

## Task Description

First, your program MUST output the total number of exact matches of AluY in the genome. 
Your program MAY output the location of each exact AluY match (chromosome, start, end).
You MUST implement exact pattern matching using only character-to-character comparisons.

Second, your program MUST output the total number of matches of AluY with up to one mismatch in the genome.
Consider these type of mismatches: replacing a character, inserting a character, deleting a character. 
Your program MAY output the location of each match (chromosome, start, end). 
You MUST implement pattern matching with mismatch using only character-to-character comparisons.

You MUST apply your program to both the GRCh38 and CHM13 human assembly, and report the output.

For both tasks, you MUST output the total runtime, and you SHOULD output peak memory usage. 
Both runtime and memory usage SHOULD be determined outside your program, 
e.g., using the time command for time, and a memory profiler like heaptrack or valgrind for memory.

Your program SHOULD read compressed (bgzipped) FASTA files due to file sizes,
and be able to handle chromosomes/sequences larger than available RAM.

# Brief Implementation Description

The entire program design consists of four main components: 
(1) downloading the corresponding genome data from NCBI; 
(2) employing a lazy-loading approach to sequentially retrieve sequence fragments as needed; 
(3) implementing tailored alignment methods that meet specific task requirements 
to analyze the matching between retrieved sequence fragments and a given reference fragment 
(i.e., AluY in this work);
and (4) monitoring the duration and memory of the entire program.

For the [function (1)](https://github.com/HaolingZHANG/CS-249/blob/main/assign0/code.py#L8), 
I utilized Python's built-in requests package to construct NCBI API commands through string concatenation, 
enabling the download of the corresponding compressed file. 
To minimize memory pressure, I set a fixed data chunk size of 1 megabyte per request.
Since the downloaded file is in ZIP format rather than GZIP, 
I directly extracted the FASTA file from the archive for subsequent analysis.
Such a FASTA file is no longer further compressed for compression purposes.

For the [function (2)](https://github.com/HaolingZHANG/CS-249/blob/main/assign0/code.py#L56), 
I used "yield" instead of "return" to sequentially produce sequence fragments 
from chromosomes/sequences larger than available RAM. 
Specifically, yield transforms the function into a generator, 
returning an iterator rather than outputting the entire result at once. 
This design ensures that the gene-loading function does not store all data in memory simultaneously, 
effectively reducing memory usage.

According to the task requirements, I implemented two types of matching functions: 
one for exact matching with no errors and another for fuzzy matching that allows a single edit error. 
To ensure result comparability, I employed three different approaches.
For exact matching, I implemented both the Naïve method (brute-force comparison) 
and the Knuth–Morris–Pratt (KMP) algorithm. For fuzzy matching, 
I extended the Naïve method by incorporating a local alignment approach based on Levenshtein distance 
and also implemented the Q-gram lemma method.
The implementation details of the two Naïve method variants can be found 
[here](https://github.com/HaolingZHANG/CS-249/blob/main/assign0/code.py#L194).
The implementation details of the KMP method are available 
[here](https://github.com/HaolingZHANG/CS-249/blob/main/assign0/code.py#L281).
The implementation details of the Q-gram lemma method can be accessed 
[here](https://github.com/HaolingZHANG/CS-249/blob/main/assign0/code.py#L337).

Since the task description specifies that the monitoring functionality should be determined externally from the program, 
we utilized the recommended "heaptrack" tool and integrated it via shell commands for execution.
Before use, it is necessary to verify whether this tool is installed.
```shell
sudo apt update
sudo apt install heaptrack
heaptrack --version
```

Potentially, HeapTrack can encounter an incompatibility issue with DWARF debugging information during execution, e.g.:

```vbnet
ERROR:Failed to create backtrace state for module /usr/local/miniconda3/bin/../lib/libstdc++.so.6: unrecognized DWARF version ...
ERROR:Failed to create backtrace state for module /usr/local/miniconda3/bin/../lib/libgcc_s.so.1: unrecognized DWARF version ...
```

We can force HeapTrack to use the system's libstdc++.so.6 and libgcc_s.so.1 to resolve this issue:

```shell
export LD_PRELOAD=/lib/x86_64-linux-gnu/libstdc++.so.6:/lib/x86_64-linux-gnu/libgcc_s.so.1
```

## Installation & Execution

You can clone this code repository (using the following command line) 
or just download this code repository with the ZIP format.

```shell
git clone https://github.com/HaolingZHANG/CS-249.git
```
No dependencies other than NumPy. 

Next, you can access the folder of the assignment 0 by

```shell
cd CS-249/assign0/
```

Here is the table of each Python script:


| class                       | script for GRCh38 | script for CHM13v2 |
|-----------------------------|-------------------|--------------------|
| exact matching using Naïve  | run_1.py          | run_2.py           |
| exact matching using KMP    | run_3.py          | run_4.py           |
| fuzzy matching using Naïve  | run_5.py          | run_6.py           |
| fuzzy matching using Q-gram | run_7.py          | run_8.py           |


To execute them, you can use the following command line:

```shell
nohup bash -c "cd /YOUR_TEMP_FOLDER/ && heaptrack python /YOUR_CODE_FOLDER/run_RUN_ID.py" > /hy-tmp/results/output-(RUN_ID).log 2>&1 &
heaptrack_print "/YOUR_TEMP_FOLDER/heaptrack.python.THE_GIVEN_PID.gz" | tail -n 20
```

In my setup, "/YOUR_TEMP_FOLDER/" and "/YOUR_CODE_FOLDER/" correspond to "./assign0/temp/" and "./assign0/", respectively.

The "RUN_ID" ranges from "1" to "8" within the "./assign0/" directory.

The "THE_GIVEN_PID", which is system-generated, corresponds to the following values:

| script   | pid  |
|----------|------|
| run_1.py | 998  |
| run_2.py | 1129 |
| run_3.py | 1258 |
| run_4.py | 1377 |
| run_5.py | 1496 |
| run_6.py | 1613 |
| run_7.py | 1738 |
| run_8.py | 1845 |


## Results
