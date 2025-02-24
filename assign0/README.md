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
nohup bash -c "cd /YOUR_TEMP_FOLDER/ && heaptrack python /YOUR_CODE_FOLDER/run_RUN_ID.py" > /YOUR_CODE_FOLDER/logs/output-(RUN_ID).log 2>&1 &
heaptrack --analyze "/YOUR_TEMP_FOLDER/heaptrack.python.THE_GIVEN_PID.gz"
```

In my setup, "/YOUR_TEMP_FOLDER/" and "/YOUR_CODE_FOLDER/" correspond to "./assign0/temp/" and "./assign0/", respectively.

The "RUN_ID" ranges from "1" to "8" within the "./assign0/" directory.

The "THE_GIVEN_PID", which is system-generated, corresponds to the following values:

| script   | pid    |
|----------|--------|
| run_1.py | 153460 |
| run_2.py | 153461 |
| run_3.py | 153462 |
| run_4.py | 153508 |
| run_5.py | 153509 |
| run_6.py | 153510 |
| run_7.py | 153511 |
| run_8.py | 153512 |


## Results

Using "run_9.py", we can obtain the results.

For the exacting search, the Naïve method and the KMP method 
both find 3 matches in GRCh38 and both find 2 matches in CHM13v2. 

| method | GRCh38 | CHM13v2 |
|--------|--------|---------|
| Naïve  | 3      | 2       |
| KMP    | 3      | 2       |

The runtime and the peak memory are:

Exact matching using Naïve for GRCh38:
- total runtime: 6981.62s.
- bytes allocated in total (ignoring deallocations): 9.91GB (1.42MB/s)
- calls to allocation functions: 821120 (117/s)
- temporary memory allocations: 2805 (0/s)
- peak heap memory consumption: 365.53MB
- peak RSS (including heaptrack overhead): 91.42GB
- total memory leaked: 3.64MB

Exact matching using Naïve for CHM13v2:
- total runtime: 6668.11s.
- bytes allocated in total (ignoring deallocations): 9.23GB (1.38MB/s)
- calls to allocation functions: 782830 (117/s)
- temporary memory allocations: 2802 (0/s)
- peak heap memory consumption: 325.80MB
- peak RSS (including heaptrack overhead): 87.88GB
- total memory leaked: 3.64MB

Exact matching using KMP for GRCh38:
- total runtime: 106033.95s.
- bytes allocated in total (ignoring deallocations): 435.58GB (4.11MB/s)
- calls to allocation functions: 127420488 (1201/s)
- temporary memory allocations: 84408660 (796/s)
- peak heap memory consumption: 365.53MB
- peak RSS (including heaptrack overhead): 91.42GB
- total memory leaked: 3.56MB

Exact matching using KMP for CHM13v2:
- total runtime: 99344.00s.
- bytes allocated in total (ignoring deallocations): 414.95GB (4.18MB/s)
- calls to allocation functions: 121437228 (1222/s)
- temporary memory allocations: 80439455 (809/s)
- peak heap memory consumption: 325.80MB
- peak RSS (including heaptrack overhead): 87.88GB
- total memory leaked: 3.64MB

Ideally, the KMP method should be faster than the Naïve method, but this is not the case. 
A possible explanation is that certain fundamental operations in the Python implementation 
may involve more complex access costs, leading to significant differences in the execution time of individual steps. 
Notably, the number of calls to dynamic memory allocation functions in KMP is higher than in the Naïve approach, 
which may support our hypothesis.

For the fuzzy search, the Naïve method and the Q-gram method 
both find 3 matches in GRCh38 and both find 2 matches in CHM13v2. 

| method | GRCh38 | CHM13v2 |
|--------|--------|---------|
| Naïve  | N.A.   | N.A.    |
| Q-gram | 37324  | 34369   |

It is important to note that these do not represent the final results,
especially considering that for a match without errors, its surrounding shifts are also included.

The runtime and the peak memory are:

Fuzzy matching using Q-gram for GRCh38:
- total runtime: 598349.68s.
- bytes allocated in total (ignoring deallocations): 55.12TB (92.12MB/s)
- calls to allocation functions: 25850921351 (43203/s)
- temporary memory allocations: 957068 (1/s)
- peak heap memory consumption: 368.13MB
- peak RSS (including heaptrack overhead): 91.59GB
- total memory leaked: 3.70MB

Fuzzy matching using Q-gram for CHM13v2:
- total runtime: 571141.50s.
- bytes allocated in total (ignoring deallocations): 53.84TB (94.27MB/s)
- calls to allocation functions: 25250852356 (44211/s)
- temporary memory allocations: 845677 (1/s)
- peak heap memory consumption: 328.40MB
- peak RSS (including heaptrack overhead): 88.05GB
- total memory leaked: 3.64MB

The Naïve method for the fuzzy search is still in progress.
