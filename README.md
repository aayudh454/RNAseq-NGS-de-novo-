# Erharta Transcriptome pipeline 

Login info: ssh aadas@vacc-user2.uvm.edu

#### Where are our Erharta files?

```
cd Erharta/Preston_JP_iLabs_14958_PR75_090820/
```
C= Control, D=Drought, F=Freezing

## Chapter 1: Trimming

#### Copying R1 and R2 files to a folder?

First you should in the folder from where you will copy then direct which folder it is going

```
[aadas@vacc-user2 ~]$ cd ~/Erharta/Preston_JP_iLabs_14958_PR75_090820/ECC1
```

Now to initiate copy-

```
cp ECC1_Illumina13_S13_L002_R1_001.fastq.gz ~/Erharta_data_analysis/Erharta_trimming/
```
 Now to create a blank script

```
[aadas@bluemoon-user2 ~]$ vi trimmomatic.sh
```

Now press **i** to insert 


Copy and paste everything present in the new script from this sample one

1) change workDIR= navigate to the the folder where you have separated the R1 and R2. Now use 'pwd' command.

use workDIR=/users/a/a/aadas/Erharta_data_analysis/Erharta_trimming/Control_1

To reopen **vim** filename.sh for further edit

```
#!/bin/bash

#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --job
#SBATCH --output=Aayudh.out
#SBATCH --mail-user=aadas@uvm.edu
#SBATCH --mail-type=ALL

###LOAD JAVA MODULE AVAILABLE FROM THE CLUSTER, YOU MAY WANT TO CHECK FIRST
module load java-sdk/sun-jdk-1.6.0.12
ulimit -s unlimited
###CHANGE THE DIRECTORY ACCORDINGLY, THE FOLLOWING SETTINGS ARE FOR MY ACCOUNT
SOFTWARE=/users/a/a/aadas/Bin/Trimmomatic-0.36
workDIR=/users/a/a/aadas/Erharta_data_analysis/Erharta_trimming/Control_1
cd $workDIR
#####TRIMMING COMMANDS AND PARAMETERS
java -jar $SOFTWARE/trimmomatic-0.36.jar PE -phred33 $workDIR/ECC1_Illumina13_S13_L002_R1_001.fastq.gz $workDIR/ECC1_Illumina13_S13_L002_R2_001.fastq.gz $workDIR/ECC1_R1_paired.fq.gz $workDIR/ECC1.R1.unpaired.fq.gz $workDIR/ECC1_R2.paired.fq.gz $workDIR/ECC1.R2.unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

```

and **esc**, **:wq** to save and quit. 

**:q!** only quit without saving.

 Make your script executable

you should be the folder where you saved your script

```
[aadas@bluemoon-user2 Ba]$ chmod 700 trimmomatic.sh
```

**700**=file's owner may read, write, and execute the file.


Submit your job and check status of your job

```
[aadas@bluemoon-user2 Ba]$ sbatch filename.sh 
```

Check your status of your job

```
[aadas@bluemoon-user2 Ba]$ squeue -u aadas
```

## What does all this things actually mean?

Based on the sequencing (paired end or signle end)  you need to follow this protocol BUT modify other factors! Don't copy deto!!!!!

**Paired End:**

java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

This will perform the following:

- Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
- Remove leading low quality or N bases (below quality 3) (LEADING:3)
- Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
- Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
- Drop reads below the 36 bases long (MINLEN:36)

**Single End:**

java -jar trimmomatic-0.35.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

This will perform the same steps, using the single-ended adapter file

i. You need to check **Submitting Jobs to the Cluster**-

TO KNOW ABOUT BATCH SYSTEM OF UVM - https://www.uvm.edu/vacc/kb/knowledge-base/write-submit-job-bluemoon/

first part of the script explains that details

ii. Now you need to specify where your software is present i.e. the Trimmomatic-0.36 which is in your main directory **/users/a/a/aadas/Trimmomatic-0.36**

iii. Specify your working directory-**/users/a/a/aadas/Ba** (because your R1 and R2 files are in Ba)

iv. TRIMMING COMMANDS AND PARAMETERS

a. change software version from as **trimmomatic-0.36**

b. Now the **first 2 files are your input file**, so after $workDIR/"name of the file" space [here R1 and R2 is the main change]
c. Last 4 files are **output files**.  

## Chapter 2: Concatenation

Now copy all **R1_paired.fq.gz** and **R2.paired.fq.gz** of all treatment into a concatenation folder 

```
-rw-r--r-- 1 aadas pi-jcpresto   94078350 Nov 24 20:59 ECC1_R1_paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto   97389322 Nov 24 21:00 ECC1_R2.paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto 1522125982 Nov 24 21:02 ECC2_R1_paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto 1567731429 Nov 24 21:02 ECC2_R2.paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto 1361193302 Nov 24 21:03 ECC3_R1_paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto 1411042061 Nov 24 21:03 ECC3_R2.paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto 1440267353 Nov 24 21:05 ECD1_R1_paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto 1484403475 Nov 24 21:06 ECD1_R2.paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto 1698297453 Nov 24 21:12 ECD2_R1_paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto 1752209079 Nov 24 21:12 ECD2_R2.paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto  611186972 Nov 24 21:20 ECD3_R1_paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto  632695412 Nov 24 21:20 ECD3_R2.paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto  130682994 Nov 24 21:35 ECF1_R1_paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto  131250331 Nov 24 21:36 ECF1_R2.paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto  492205587 Nov 24 21:37 ECF2_R1_paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto  513204211 Nov 24 21:37 ECF2_R2.paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto  287248833 Nov 24 21:37 ECF3_R1_paired.fq.gz
-rw-r--r-- 1 aadas pi-jcpresto  299917169 Nov 24 21:38 ECF3_R2.paired.fq.gz
```
**Now concatenate all the R1.** it takes 30in or more, wait until size increase 

```
[aadas@vacc-user2 Concatenation_Erharta]$ zcat *R1_paired.fq.gz > Erharta.R1.trimmed.fq &
[1] 31773
```
**Now concatenate all the R2**

```
[aadas@vacc-user2 Concatenation_Erharta]$ zcat *R2.paired.fq.gz > Erharta.R2.trimmed.fq &
[2] 32386
```

**After finishing the concatenation check the "sequence header for both R1 and R2"-it should be same**

```
[aadas@vacc-user2 Concatenation_Erharta]$ grep -c "@" Erharta.R1.trimmed.fq
150563683
[aadas@vacc-user2 Concatenation_Erharta]$ grep -c "@" Erharta.R2.trimmed.fq
150563683
```
both shows 150563683; That means R1 and R2 has same reads.


## Chapter 3: Assembly

Trinity combines three independent software modules: **Inchworm, Chrysalis, and Butterfly**, applied sequentially to process large volumes of RNA-seq reads. Trinity partitions the sequence data into many individual de Bruijn graphs, each representing the transcriptional complexity at a given gene or locus, and then processes each graph independently to extract full-length splicing isoforms and to tease apart transcripts derived from paralogous genes. Briefly, the process works like so:

**Inchworm** assembles the RNA-seq data into the unique sequences of transcripts, often generating full-length transcripts for a dominant isoform, but then reports just the unique portions of alternatively spliced transcripts.

**Chrysalis** clusters the Inchworm contigs into clusters and constructs complete de Bruijn graphs for each cluster. Each cluster represents the full transcriptonal complexity for a given gene (or sets of genes that share sequences in common). Chrysalis then partitions the full read set among these disjoint graphs.

**Butterfly** then processes the individual graphs in parallel, tracing the paths that reads and pairs of reads take within the graph, ultimately reporting full-length transcripts for alternatively spliced isoforms, and teasing apart transcripts that corresponds to paralogous genes.



Move concatenated R1 and R2 trimmed files to a different folder

```
[aadas@vacc-user2 Concatenation_Erharta]$ cp Erharta.R1.trimmed.fq /users/a/a/aadas/Erharta_data_analysis/assembly_erharta
[aadas@vacc-user2 Concatenation_Erharta]$ cp Erharta.R2.trimmed.fq /users/a/a/aadas/Erharta_data_analysis/assembly_erharta
```

**Sript for assembly for single ended by trinity 2.1.1**

```
#!/bin/bash

#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=30:00:00
#SBATCH --mem=256G
#SBATCH --job
#SBATCH --output=Aayudh_assembly.out
#SBATCH --mail-user=aadas@uvm.edu
#SBATCH --mail-type=ALL

module load samtools-1.10-gcc-7.3.0-pdbkohx
module load bowtie2-2.3.5.1-gcc-7.3.0-ycvyrzc
export PATH="/users/a/a/aadas/Bin/bowtie-1.1.1:$PATH"
#ulimit -s unlimited

SOFTWAREDIR=/users/a/a/aadas/Bin/trinityrnaseq-2.1.1
WORKINGDIR=/users/a/a/aadas/Erharta_data_analysis/assembly_erharta
cd $WORKINGDIR

/users/a/a/aadas/Bin/trinityrnaseq-2.1.1/Trinity --seqType fq --normalize_reads --max_memory 256G --left /users/a/a/aadas/Erharta_data_analysis/assembly_erharta/ Erharta.R1.trimmed.fq –right /users/a/a/aadas/Erharta_data_analysis/assembly_erharta/Erharta.R2.trimmed.fq --CPU 24
```

Make your script executable. You should be the folder where you saved your script
```
[aadas@vacc-user2 assembly_erharta]$ chmod 700 assembly_2.11_erharta.sh
```
Submit your job and check status of your job
```
[aadas@vacc-user2 assembly_erharta]$ sbatch assembly_2.11_erharta.sh 
```
Check your status of your job
```
[aadas@vacc-user2 assembly_erharta]$ squeue -u aadas
```

**What all those mean?**
- Use all reads from an individual (all conditions) to capture most genes
- Read files may be gzipped (as in this example) or not (then they should not have the “.gz” ending)
- Paired-end reads specified with **--left** and -**-right**. If only single-end, use **--single** instead.
- 256G is the maximum memory to be used at any stage which allows memory limitation (jellyfish, sorting, etc.)
- At most 24 CPU cores will be used in any stage.

It might take at least 10hours and if takes more than 30hrs just resubmit it again in the server.


```
[aadas@vacc-user2 assembly_erharta]$ grep ">" Trinity.fasta | less
[aadas@vacc-user2 assembly_erharta]$ grep ">" Trinity.fasta | sed "s/_i[0-9]\{1,2\} len.*//g" | less
[aadas@vacc-user2 assembly_erharta]$ grep ">" Trinity.fasta | sed "s/_i[0-9]\{1,2\} len.*//g" | sort -u | less
[aadas@vacc-user2 assembly_erharta]$ grep ">" Trinity.fasta | sed "s/_i[0-9]\{1,2\} len.*//g" | sort -u | wc -l
```
**95447**

**check no. of seq**

```
[aadas@vacc-user2 assembly_erharta]$ grep ">" Trinity.fasta -c
```
**158806**

## Chapter 3: Transcript quantification by RSEM

RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome

Script for transcript quantification.

* Make sure your **Trinity.fasta** file in the directory
* Make sure that **.trimmed.fq.gz** is also in the FASTQ_DIR


```
#!/bin/bash

#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=30:00:00
#SBATCH --mem=8G
#SBATCH --job
#SBATCH --output=Aayudh_RSEM.out
#SBATCH --mail-user=aadas@uvm.edu
#SBATCH --mail-type=ALL

module load samtools-1.10-gcc-7.3.0-pdbkohx

TRINITY_HOME=/users/a/a/aadas/Bin/trinityrnaseq-2.1.1

export PATH=/users/a/a/aadas/Bin/bowtie-1.1.1:$PATH
export PATH=/users/a/a/aadas/Bin/RSEM-1.2.19:$PATH

FASTQ_DIR=/users/a/a/aadas/Erharta_data_analysis/rsem_erharta/trinity_out_dir

cd $FASTQ_DIR

$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left $FASTQ_DIR/ECC1_R1_paired.fq.gz --right $FASTQ_DIR/ECC1_R2.paired.fq.gz --est_method RSEM --aln_method bowtie --thread_count 4 --SS_lib_type RF --trinity_mode --prep_reference --output_dir $FASTQ_DIR --output_prefix Erharta_Control01


$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left $FASTQ_DIR/ECC2_R1_paired.fq.gz --right $FASTQ_DIR/ECC2_R2.paired.fq.gz --est_method RSEM --aln_method bowtie --thread_count 4 --SS_lib_type RF --trinity_mode --output_dir $FASTQ_DIR --output_prefix Erharta_Control02


$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left $FASTQ_DIR/ECC3_R1_paired.fq.gz --right $FASTQ_DIR/ECC3_R2.paired.fq.gz --est_method RSEM --aln_method bowtie --thread_count 4 --SS_lib_type RF --trinity_mode --output_dir $FASTQ_DIR --output_prefix Erharta_Control03
```

**Open R in UNIX**

You should be in the directory where your files seats.

```
[aadas@vacc-user2 results_erharta]$ module load r-3.6.3-gcc-7.3.0-qo3xjgm
[aadas@vacc-user2 results_erharta]$ R
```

Type 'q()' to quit R.

```
> q()
Save workspace image? [y/n/c]: n
```
The above `R` commands list the top **10 highest expressed genes** (shown below).

```
> data = read.table("Erharta_Control01.genes.results", header=T, stringsAsFactors=F)
> idx = order(data[,"TPM"], decreasing=T)
> data[idx[1:10], c("gene_id", "expected_count", "TPM")]

> data[idx[1:10], c("gene_id", "expected_count", "TPM")]
                    gene_id expected_count      TPM
63364 TRINITY_DN47580_c1_g2        2287.17 81161.47
63365 TRINITY_DN47580_c1_g3       19087.94 57618.49
65712 TRINITY_DN48068_c6_g1       15660.55 53514.80
65709 TRINITY_DN48068_c4_g1        4400.00 39752.32
50273 TRINITY_DN44444_c4_g1       15965.00 32329.55
41409 TRINITY_DN41740_c2_g1       34639.00 21155.95
63363 TRINITY_DN47580_c1_g1        5559.00 19867.52
39675 TRINITY_DN41080_c0_g1       11087.00 15337.50
43098 TRINITY_DN42341_c0_g1        1177.00 12890.06
39676 TRINITY_DN41080_c1_g1       10341.00 12117.45
```
**Dowload files to your desktop**

Go to your Mac portal then do this-

```
Aayudhs-MacBook-Pro:~ aayudhdas$ scp aadas@vacc-user2.uvm.edu:/users/a/a/aadas/Erharta_data_analysis/rsem_erharta/results_erharta/Erharta_Control01.genes.results ~/Desktop/
````


