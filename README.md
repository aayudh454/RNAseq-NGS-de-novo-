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
Run for all three treatments and then copy them to a different folder for next analysis

```
-rw-r--r-- 1 aadas pi-jcpresto 8851655 Nov 27 21:55 Erharta_Control01.genes.results
-rw-r--r-- 1 aadas pi-jcpresto 8862588 Nov 27 21:55 Erharta_Control02.genes.results
-rw-r--r-- 1 aadas pi-jcpresto 8870068 Nov 27 21:55 Erharta_Control03.genes.results
-rw-r--r-- 1 aadas pi-jcpresto 8866060 Nov 27 21:55 Erharta_Drought01.genes.results
-rw-r--r-- 1 aadas pi-jcpresto 8869023 Nov 27 21:55 Erharta_Drought02.genes.results
-rw-r--r-- 1 aadas pi-jcpresto 8856867 Nov 27 21:55 Erharta_Drought03.genes.results
-rw-r--r-- 1 aadas pi-jcpresto 8850300 Nov 27 21:55 Erharta_Freezing01.genes.results
-rw-r--r-- 1 aadas pi-jcpresto 8857487 Nov 27 21:55 Erharta_Freezing02.genes.results
-rw-r--r-- 1 aadas pi-jcpresto 8851219 Nov 27 21:55 Erharta_Freezing03.genes.results
```

## Chapter 4: Build Transcript and Gene Expression Matrices

Terms:

**<u>FPKM</u>: fragments per kilobase transcript length per million fragments mapped**

**<u>TPM</u>: transcripts per million transcripts**

Using the transcript and gene-level abundance estimates for each of your samples, construct a matrix of counts and a matrix of normalized expression values using the following script:

#### For genes of *Erharta* control, drought and freezing

```
[aadas@bluemoon-user2 rsem_npulBdis]$ ~/Bin/trinityrnaseq-2.1.1/util/abundance_estimates_to_matrix.pl --est_method RSEM Erharta_Control01.genes.results Erharta_Control02.genes.results Erharta_Control03.genes.results Erharta_Drought01.genes.results Erharta_Drought02.genes.results Erharta_Drought03.genes.results Erharta_Freezing01.genes.results Erharta_Freezing02.genes.results Erharta_Freezing03.genes.results --out_prefix Erharta.genes
 
```
Now smilar way you can do this for isoforms- use **.isoforms.results** 

#### Counting Numbers of Expressed Transcripts or Genes
```
~/Bin/trinityrnaseq-2.1.1/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \ Erharta.genes.TPM.not_cross_norm | tee Erharta.genes.TPM.not_cross_norm.counts_by_min_TPM
```
The above table indicates that we have 95,447 'genes' that are expressed by at least 1 TPM in any one of the many samples in this expression matrix.

Plotting the number of 'genes' (or 'transcripts') as a function of minimum TPM threshold, we can see that the vast majority of all expressed features have very little expression support. Using R (or your own favorite data analysis package), we might extrapolate the number of expressed 'genes' based on the trend prior to the massive influx of lowly expressed transcripts:

Copy the **TPM.not_cross_norm.counts_by_min_TPM** file to MAC and open R in MAC.

```
setwd("~/Desktop/Erharta_data_analysis")
list.files()
data = read.table("Erharta.genes.TPM.not_cross_norm.counts_by_min_TPM", header=T)

tiff("Erharta_gene_matrix.tiff", width = 6.65, height = 3.75, units = 'in', res = 300)
plot(data, xlim=c(-100,0), ylim=c(0,100000), t='b')

# extract the data between 10 TPM and 100 TPM
filt_data = data[data[,1] > -100 & data[,1] < -10,] 
# perform a linear regression on this filtered subset of the data
fit = lm(filt_data[,2] ~ filt_data[,1])
print(fit)
dev.off()

# add the linear regression line to the plot 
abline(fit, col='green', lwd=3)
```
The linear regression allows us to extrapolate (based on the Y-intercept) that we have 13965 'genes', which is a far better guess than our count of 95,447 'genes' having at least 1 TPM in any sample, and certainly better than the 1.4 million 'genes' that were assembled. 


## Chapeter 5: Differential Expression Analysis (Voom and DeSeq2)

Copy **Erharta.genes.counts.matrix** from RSEM output and create a text file with all the parameters

```
conditionA   Erharta_Control01
conditionA   Erharta_Control02
conditionA   Erharta_Control03

conditionB   Erharta_Drought01
conditionB   Erharta_Drought02
conditionB   Erharta_Drought03

conditionC   Erharta_Freezing01
conditionC   Erharta_Freezing02
conditionC   Erharta_Freezing03
```

*First lets get 'R' working* Check what R version is available

```
[aadas@vacc-user2 differential_gene_expression]$ module avail
[aadas@vacc-user2 differential_gene_expression]$ module load r-3.6.3-gcc-7.3.0-qo3xjgm
```
Then type 'R' and load the libraries. Now install all the packages by either install.packages ("package name"). Specially EdgeR and Deseq2 you have to follow this below after opening R in the portal
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

BiocManager::install('limma')

BiocManager::install('ctc')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("qvalue")
```
After installing all packages call them.
```
library('limma')
library('edgeR')
library('ctc')
library('DESeq2')
library('Biobase')
library('gplots')
library('ape')
library ('qvalue')
```

**voom method for DE**

Any of the available methods support analyses containing biological replicates. Here, for example, we again choose voom within the limma package.
```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Erharta.genes.counts.matrix --method voom --samples_file gene_expression.txt
```
Now scp the pdfs to your MAC.

**Extracting and clustering differentially expressed transcripts (HEATMAP)**

An initial step in analyzing differential expression is to extract those transcripts that are most differentially expressed (most significant FDR and fold-changes) and to cluster the transcripts according to their patterns of differential expression across the samples. 

This will extract all genes that have P-values at most 1e-3 and are at least 2^2 fold differentially expressed. For each of the earlier pairwise DE comparisons, this step will generate the following files:

You should do this **inside the voom folder**. 
```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix Erharta.genes.counts.matrix -P 1e-3 -C 2 --samples gene_expression.txt
```

**Automatically Partitioning Genes into Expression Clusters**

There are three different methods for partitioning genes into clusters:

- use K-means clustering to define K gene sets. (use the -K parameter). This does not leverage the already hierarchically clustered genes as shown in the heatmap, and instead uses a least-sum-of-squares method to define exactly k gene clusters.
- cut the hierarchically clustered genes (as shown in the heatmap) into exactly K clusters.
- (Recommended) cut the hierarchically clustered gene tree at --Ptree percent height of the tree.

Do it inside voom folder and call module load r-3.6.3-gcc-7.3.0-qo3xjgm

```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl -R diffExpr.P1e-3_C2.matrix.RData --Ptree 60
```

**4. Principal Component Analysis (PCA)**

The --prin_comp 3 indicates that the first three principal components will be plotted, as shown above, with PC1 vs. PC2 and PC2 vs. PC3. In this example, the replicates cluster tightly according to sample type, which is very reassuring.
```
~/Bin/trinityrnaseq-2.1.1/Analysis/DifferentialExpression/PtR --matrix Erharta.genes.counts.matrix \
    -s gene_expression.txt --log2 --prin_comp 3
```
If you have replicates that are clear outliers, you might consider removing them from your study as potential confounders. If it's clear that you have a [batch effect](http://www.nature.com/nrg/journal/v11/n10/full/nrg2825.html), you'll want to eliminate the batch effect during your downstream analysis of differential expression.


## Chapeter 6: Extended version of DESeq2 analysis in R

**scp the gene.matrix to your MAC and then arrange the treatment to the fashion**. For installing softwares if error is coming do all the steps of installing that you did in the server.
**TREATMENT DETAILS**
```
sample	cond	pop	rep
Erharta_Control01	Control	ErhartaC	1
Erharta_Control02	Control	ErhartaC	2
Erharta_Control03	Control	ErhartaC	3
Erharta_Drought01	Drought	ErhartaD	1
Erharta_Drought02	Drought	ErhartaD	2
Erharta_Drought03	Drought	ErhartaD	3
Erharta_Freezing01	Freezing	ErhartaF	1
Erharta_Freezing02	Freezing	ErhartaF	2
Erharta_Freezing03	Freezing	ErhartaF	3
```
**Run DESeq2**
```
setwd("~/Desktop/Erharta_data_analysis/Deseq2")
list.files()
library("DESeq2")
library("ggplot2")

#Erharta Control_vs_Freezing
countsTable <- read.delim('Erharta.genes.counts.matrix.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
#Choose control and freezing replicates column from the table
countData=as.matrix(countsTable[,c(1,2,3,7,8,9)]) #Choose control and freezing replicates from the table
storage.mode(countData) = "integer"
head(countData)

conds <- read.delim("TREATMENT DETAILS.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds[c(1,2,3,7,8,9),])
head(colData)

dim(countData)
dim(colData)
#model 
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ pop)
dds
dim(dds)
# Filtering to remove rows with 0 reads
dds <- dds[ rowSums(counts(dds)) > 1, ]
dim(dds)
dds <- DESeq(dds) 
res <- results(dds)
#sorts according to pvalue
res <- res[order(res$padj),]
head(res)
summary(res)
summary(sub_data)
#Subsetting based on pvalue 0.01
x.sub <- subset(res, padj < 0.01)
sub_data=as.data.frame(x.sub)
y.sub <- subset(res, pvalue < 0.001)
y.sub1 = as.data.frame(y.sub)
z.sub <- subset(res, padj < 0.0001)
z.sub1 =  as.data.frame(z.sub)
#Up and Downregulated genes
table(sign(y.sub1$log2FoldChange))
#save as csv
write.csv(z.sub1, file = "z.sub1.csv")

#MA plot

plotMA(res, main="Erharta_controlvsFreezing", ylim=c(-12,12), xlim=c(1e-01, 1e+04))
abline(h=c(-1:1), col="red")
#Volcano plot
tiff("Erharta_controlvsFreezing.tiff", width = 9.97, height = 6.5, units = 'in', res = 300)
par(family="Times")
volcanoData <- as.data.frame(sub_data[,c(2,5,6)])
with(volcanoData, plot(log2FoldChange, -log10(pvalue), pch=20,lwd=4, main="Erharta Control vs Freezing", xlim=c(-12,13), ylim = c(0,30)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(volcanoData, padj<.001 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(volcanoData, abs(log2FoldChange)>5), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(volcanoData, padj<.001 & abs(log2FoldChange)>5), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(volcanoData, padj<.0001 & abs(log2FoldChange)>5), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
legend( x="topright", 
        legend=c("padj<0.001 and log2(FC)>5","log2(FC)>5","padj<0.001","padj<0.01 and log2(FC)>2"), 
        col=c("green","orange","red","black"), lwd=4, lty=c(NA,NA,NA,NA), 
        pch=c(20,20,20,20), merge=FALSE, cex = 0.75)
dev.off()
```
## Chapeter 7: Transdecoder and Gene Annotation

***Step 1: extract the long open reading frames by Transdecoder***

```
#!/bin/bash

#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=30:00:00
#SBATCH --mem=10G
#SBATCH --job
#SBATCH --output=Aayudh_transDecoder-step1.out
#SBATCH --mail-user=aadas@uvm.edu
#SBATCH --mail-type=ALL


export PATH="/users/a/a/aadas/Bin/TransDecoder-3.0.1/transdecoder_plugins/cdhit:$PATH"
export PATH="/users/a/a/aadas/Bin/TransDecoder-3.0.1:$PATH"
export PATH="/users/a/a/aadas/Bin/hmmer-3.1b2-linux-intel-x86_64/binaries:$PATH"
export PATH="/users/a/a/aadas/Bin/ncbi-blast-2.6.0+/bin:$PATH"

transDecoder_dir=/users/a/a/aadas/Bin/TransDecoder-3.0.1
INPUT_DIR=/users/a/a/aadas/Erharta_data_analysis/transdecoder
cd $INPUT_DIR

$transDecoder_dir/TransDecoder.LongOrfs -t $INPUT_DIR/Trinity.fasta
```


***Step 2 (run two scripts):BlastP and pfam search together***

This single line using the blastp command below will compare your transcript fasta file

(-query) to the already formatted uniref90 database (-db).

You can enter 'blastp --help' for a list of the parameters.

We choose the tab-delimited output format (6) and to only help the top hit (-max_target_seqs) and only if it has a minimum evalue of 0.001.

**batch_blastp_hmmscan.pl SCRIPT**

```
#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

 
my $file;		#BLAST query sequences
my $size= 2000;		#Number of sequences per task
my $type;		#Blast program.
my $database;		#path to database
my $eval= 1e-5;		#BLAST e-value cutoff
my $outputFormat= 6;	#BLAST output format
my $outputDir="blast_out";	#output directory
my $help;
my $input_dir=`pwd`;
my $blast_dir="/users/a/a/aadas/Bin/ncbi-blast-2.6.0+/bin";
my $hmm_dir="/users/a/a/aadas/Bin/hmmer-3.1b2-linux-intel-x86_64/binaries";
my $hmmDB_dir="/users/a/a/aadas/blastp_Brachyleytrum/database";
my $hmm="hmmscan";
my $hmmoutputDir="hmmscan_out";
GetOptions(
	'query=s' 	=> \$file,
	'num_seqs=i'	=> \$size,
	'program=s'	=> \$type,
	'database=s'	=> \$database,
	'eval=f' 	=> \$eval,
	'm=i'		=> \$outputFormat,
	'output_dir=s'	=> \$outputDir,
#	'rcc_queue=s'   => \$queue,
#	'combine'	=> \$autoCombine,
	'help'		=> \$help
);

my $usage = <<__EOUSAGE__;

###################################################################################
#	Batch Blast: Task Array
###################################################################################
#
#  --query <string>		File containing query sequences
# 
#  --program <string>		The BLAST program to use
#
#  --database <string>		The location of the BLAST database
#
# Optional:
# 
#  --num_seqs <integer>		The number of sequences per data sub-set
#				Default: 1000
#
#  --eval <float>           	The BLAST e-value cutoff
#				Default: 1e-10
#
#  --m <integer>		BLAST output format (8 for tablular)
#				Default: 8 - Tabular Output
#
#  --output_dir <string>	Output directory for BLAST results
#				Default: blast_out
#
#  --combine			Automatically combine and remove the split.*
#				directories.
#				Default: False
#
#  --job_name_prefix		A prefix to differentiate this batch_blast run
#				from another.
#				
#  --rcc_queue			SGE queue to run the BLAST job on
#				Default: rcc-30d
#
#  --help  
###################################################################################
#  
#  To pass other arguments to blastall use -- <args> AFTER all required
#  arguments.
#
###################################################################################

__EOUSAGE__
;

if(!$file || !$type || !$database || $help){
	die($usage);
}

my $seqid;
my $seq;
 
my $seq_counter=0;
my $split_count=1;
 
open my $infile, "<", $file;
while(<$infile>){
	chomp;
	#is this line a new sequence header?
	if(/^>/){
		#if there is a storred sequence then write it out to file first
		if($seqid){
			mkdir "split.".$split_count;
			chdir ("split.$split_count");
			open my $OUT, ">>", "split.$split_count.fasta";
			print $OUT "$seqid\n$seq\n";
			close $OUT;
			if($seq_counter == $size){
				$seq_counter = 0;
				$split_count++;
				
			}
			chdir ("../");
		}
		$seq=();
		$seqid = $_;
		$seq_counter++;
	}
	#if not then keep building the sequence
	else{
		$seq .= $_;
	}
}
#necessarily this loop exits before the last sequence is written. Write it now.
if($seqid){
	mkdir "split.".$split_count;
	chdir ("split.$split_count");
	open my $OUTsplit, ">>", "split.$split_count.fasta";
	print $OUTsplit "$seqid\n$seq\n";
	close $OUTsplit;
	if($seq_counter == $size){
		$seq_counter = 0;
		$split_count++;
	}
	chdir ("../");
}
close $infile;
#if the output directory doesn't exist then make it
if(! -e $outputDir){
	`mkdir $outputDir`;
}

#print the task-array script
for(my $i=1;$i<=$split_count;$i++){
     open my $SUB_SCRIPTS, ">", "$type-part-$i.sh" or die();
     print $SUB_SCRIPTS "#!/bin/bash\n",
			"#SBATCH --job-name=out.$type.part-$i\n",
			"#SBATCH --nodes=1\n",
			"#SBATCH --ntasks=24\n",
			"#SBATCH --mem=10G\n",
			"#SBATCH --time=30:00:00\n",
			"#SBATCH --mail-user=aadas\@uvm.edu\n",
			"#SBATCH --mail-type=ALL\n",
			"\n",
			"INPUT_DIR=$input_dir\n",
			"BLAST_DIR=$blast_dir\n",
			"cd \$INPUT_DIR\n",
			"\n",
			"\$BLAST_DIR/$type ",
			"-query \$INPUT_DIR/split.$i/split.$i.fasta ",
			"-db $database ",
			"-outfmt $outputFormat ",
			"-evalue $eval ",
			"-num_threads 1 ",
			"-max_target_seqs 1 ",
			"\> $outputDir/split.$i.$type",
			"\n",
			"exit\n";
}
#print the task-array scripts hmmscan
if(! -e $hmmoutputDir){
        `mkdir $hmmoutputDir`;
}

for(my $i=1;$i<=$split_count;$i++){
     open my $HMM_SCRIPTS, ">", "$hmm-part-$i.sh" or die();
     print $HMM_SCRIPTS "#!/bin/bash\n",
			"#SBATCH --job-name=out.$type.part-$i\n",
                        "#SBATCH --nodes=1\n",
                        "#SBATCH --ntasks=24\n",
                        "#SBATCH --mem=10G\n",
                        "#SBATCH --time=30:00:00\n",
                        "#SBATCH --mail-user=aadas\@uvm.edu\n",
                        "#SBATCH --mail-type=ALL\n",
                        "\n",
                        "INPUT_DIR=$input_dir\n",
                        "HMM_DIR=$hmm_dir\n",
                        "cd \$INPUT_DIR\n",
                        "\n",
                        "\$HMM_DIR/$hmm ",
                        "--cpu 1 ",
                        "--domtblout $hmmoutputDir/split.$i.domtblout ",
                        "$hmmDB_dir/Pfam-A.hmm ",
                        "\$INPUT_DIR/split.$i/split.$i.fasta ",
                        "\n",
                        "exit\n";
}
#and submit it


```

Then (**When you run the 3rd script; generally it's submitting 300 jobs to VACC but it's not going to run; So, comment out blast-then it's only submitting hmm scan; Even if doesn't work then edit blast-part* and manually put part-1* then part-2*  **)

```
#!/bin/bash
cd `pwd`

perl batch_blastp_hmmscan.pl -q /users/a/a/aadas/Erharta_data_analysis/transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep -p blastp -d /users/a/a/aadas/blastp_Brachyleytrum/database/uniprot_sprot.pep

chmod 700 *.sh
for i in blastp-part*; do
  sbatch $i &
done

for i in hmmscan-part*; do
  sbatch $i &
done
```
**Run the job with the below code. NOT by sbatch**

```
./run_blastp_hmmscan.sh
```

