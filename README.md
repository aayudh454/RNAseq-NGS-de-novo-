# Erharta Transcriptome pipeline 

Login info: ssh aadas@vacc-user2.uvm.edu

#### Where are our Erharta files?

```
cd Erharta/Preston_JP_iLabs_14958_PR75_090820/
```
C= Control, D=Drought, F=Freezing

## Table of contents    
* [Page 1: 2020-12-05](#id-section1). Chapter 1: Moving files and trimming by Trimmomatic-0.36

* [Page 2: 2020-12-05](#id-section2). Chapter 2: Concatenation

* [Page 3 2020-12-06](#id-section3). Chapter 3: Assembly by Trinity 2.1.1

* [Page 4 2020-12-07](#id-section4). Chapter 4: Transcript quantification by RSEM

* [Page 5 2020-12-08](#id-section5). Chapter 5: Build Transcript and Gene Expression Matrices

* [Page 6 2020-12-08](#id-section6). Chapeter 6: Differential Expression Analysis (Voom and DeSeq2)

* [Page 7 2020-12-14](#id-section7). Chapeter 7: Extended version of DESeq2 analysis in R

* [Page 8 2020-12-16](#id-section8). Chapeter 8: Transdecoder and Gene Annotation (blastp)

* [Page 9 2020-12-22](#id-section9). Chapter 9: Go annotation (uniprot)

* [Page 10 2020-12-22](#id-section10). Chapter10: Orthofinder

* [Page 11 2021-04-06](#id-section11). Chapter11: Enrichment Analysis

* [Page 12 2021-04-06](#id-section12). Chapter12: KEGG Analysis

* [Page 13 2023-05-17](#id-section13). Chapter13: NGS FAQ




------
<div id='id-section1'/>

## Chapter 1: Moving files and trimming by Trimmomatic-0.36

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

-ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
-SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
-LEADING: Cut bases off the start of a read, if below a threshold quality
-TRAILING: Cut bases off the end of a read, if below a threshold quality
-CROP: Cut the read to a specified length
-HEADCROP: Cut the specified number of bases from the start of the read
-MINLEN: Drop the read if it is below a specified length
-TOPHRED33: Convert quality scores to Phred-33
-TOPHRED64: Convert quality scores to Phred-64

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

-----
<div id='id-section2'/>


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

***Another way of counting fasta files***

```
 echo $(cat Brachy.R1.trimmed.fq|wc -l)/4|bc
```

------

<div id='id-section3'/>


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

/users/a/a/aadas/Bin/trinityrnaseq-2.1.1/Trinity --seqType fq --normalize_reads --max_memory 256G --left /users/a/a/aadas/Erharta_data_analysis/assembly_erharta/Erharta.R1.trimmed.fq --right /users/a/a/aadas/Erharta_data_analysis/assembly_erharta/Erharta.R2.trimmed.fq --CPU 24
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

------

<div id='id-section4'/>

## Chapter 4: Transcript quantification by RSEM

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

If you want to download the **whole folder** then
```
scp -r aadas@vacc-user2.uvm.edu://users/a/a/aadas/Nassella_data_analysis/differential_gene_expression/voom.3044.dir ~/Desktop/Nassella_Brachypodium_data\ analysis/Nassella/
```
**MAC to UNix copy**
```
Aayudhs-MacBook-Pro:Downloads aayudhdas$ scp -rp /Users/aayudhdas/Downloads/fastme-2.1.5.tar.gz aadas@vacc-user2.uvm.edu:/users/a/a/aadas/Bin
```

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

------

<div id='id-section5'/>

## Chapter 5: Build Transcript and Gene Expression Matrices

Terms:

**<u>FPKM</u>: fragments per kilobase transcript length per million fragments mapped**

**<u>TPM</u>: transcripts per million transcripts**

Using the transcript and gene-level abundance estimates for each of your samples, construct a matrix of counts and a matrix of normalized expression values using the following script:

**Now first load R. Without loading R you will have full of erros**
```
[aadas@vacc-user2 gene_expression_matrices]$ module load r-3.6.3-gcc-7.3.0-qo3xjgm
```
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


------

<div id='id-section6'/>


## Chapeter 6: Differential Expression Analysis (Voom and DeSeq2)

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

------

<div id='id-section7'/>

## Chapeter 7: Extended version of DESeq2 analysis in R

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

#Erharta Control_vs_Freezing (Same as example of Nassella)
setwd("~/OneDrive - University of Vermont/Freezing tolerance/Nassella Freezing")
list.files()
library("DESeq2")
library("ggplot2")

#Nassella Control_vs_Freezing
countsTable <- read.delim('Nassella.genes.counts.matrix.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
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
#sorts according to padj
res <- res[order(res$padj),]
head(res)
summary(res)

#Subsetting based on padj 0.05 and LogFC2
x.sub <- subset(res, padj < 0.05)
y.sub <- as.data.frame(x.sub)

y.sub1 <- subset(y.sub, log2FoldChange>2)
z.sub <- subset(y.sub, log2FoldChange < -2)
sub_data <- rbind (y.sub1,z.sub)

summary(sub_data)
dim(sub_data)
head (sub_data)

#Up and Downregulated genes
table(sign(sub_data$log2FoldChange))

#save as csv
write.csv(sub_data, file = "Nassella_freezing_p0.05_LogFC2.csv")


#plotMA(res, main="Nassella_controlvscold", ylim=c(-12,12), xlim=c(1e-01, 1e+04))
#abline(h=c(-1:1), col="red")
#Volcano plot
tiff("Nassella_controlvsfreezing.tiff", width = 5, height = 5, units = 'in', res = 300)
par(family="Times")
par(pty="s")
volcanoData <- as.data.frame(sub_data[,c(2,5,6)])
with(volcanoData, plot(log2FoldChange, -log10(padj), pch=20,lwd=4, main="Nassella Control vs Freezing", 
                       xlim=c(-12,12), ylim = c(0,50)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(volcanoData, padj<.001 ), points(log2FoldChange, -log10(padj), pch=20, col="red"))
with(subset(volcanoData, abs(log2FoldChange)>5), points(log2FoldChange, -log10(padj), pch=20, col="orange"))
with(subset(volcanoData, padj<.001 & abs(log2FoldChange)>5), points(log2FoldChange, -log10(padj), pch=20, col="green"))
with(subset(volcanoData, padj<.0001 & abs(log2FoldChange)>5), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
#legend( x="topright", 
        #legend=c("padj<0.001 and log2(FC)>2","log2(FC)>5","padj<0.001","padj<0.001 and log2(FC)>5"), 
        #col=c("green","orange","red","blue"), lwd=4, lty=c(NA,NA,NA,NA), 
        #pch=c(20,20,20,20), merge=FALSE, cex = 0.75)
dev.off()


#merge with blast output
list.files()
data <- read.csv("Nassella_blastp.csv")
data1 <- read.csv("Nassella_freezing_p0.05_LogFC2.csv")
head(data)
head(data1)
#merge data by ID
merged_data<- merge(data1,data, by="Trinity_ID")
head(merged_data)
dim(merged_data)

#Delete repeated column
library(dplyr)
merged_data_final <- merged_data %>% distinct
head(merged_data_final)
dim(merged_data_final)

merged_data_final <- as.data.frame(merged_data)
head(merged_data_final)

#Save as csv file
write.csv(merged_data_final, file = "Nassella_freezing_uniprot.csv")

##Now follow the steps for GO terms and GO terms csv.

#Now make all data in one file 
list.files()
data <- read.csv("Nassella_freezing_uniprot.csv")
data1 <- read.csv("Nassella_freezing_GO.csv")
#merge data by ID
merged_data<- merge(data,data1, by="Uniprot_ID")
head(merged_data)
dim(merged_data)
merged_data_final <- as.data.frame(merged_data)
head(merged_data_final)
dim(merged_data_final)

library(data.table)
MASTER_data <- unique(setDT(merged_data_final)[order(Uniprot_ID, -Uniprot_ID)], by = "Uniprot_ID")

#Save as csv file
write.csv(MASTER_data, file = "Nassella_FREEZING_MasterDATA.csv")

Nas_cold <- read.csv("Nassella_FREEZING_MasterDATA.csv")
#Up and Downregulated genes
table(sign(Nas_cold$log2FoldChange))

```

![alt text](https://github.com/aayudh454/Erharta/blob/main/Freezing-Drought_deseq.png)

**Making heatmap using R using log fold change values**

```
setwd("~/OneDrive - University of Vermont/Freezing tolerance/Freezing conserved")
list.files()

library("gplots")

data <-  read.table("Ehr_Npul_freezing_heatmap.csv", sep=",", header=TRUE, row.names=1)
mat_data <- data.matrix(data[,1:ncol(data)])
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

setwd("~/OneDrive - University of Vermont/Freezing tolerance/Freezing conserved")
tiff("Ehr_Npul_freezing_conserved.tiff", width = 8, height = 15, units = 'in', res = 300)
par(family="Times")
heatmap.2(mat_data,
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(6,6),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",cexRow=1)  
dev.off()
```

------

<div id='id-section8'/>


## Chapeter 8: Transdecoder and Gene Annotation

***Step 1: extract the long open reading frames by Transdecoder***

Move Trinity.fasta to a separate folder

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

*Before running the 3rd script, copy the trinity fasta and longest_orfs ../ one folder back and run the 3rd script from their NOT inside of Trinity.fasta.transdecoder_dir*

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

**Step 3 - Final step: predict the likely coding regions**

```
#!/bin/bash

#SBATCH --partition=bigmem
#SBATCH --job-name=out.finalstep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH --mail-user=aadas@uvm.edu
#SBATCH --mail-type=ALL


export PATH="/users/a/a/aadas/Bin/TransDecoder-3.0.1/transdecoder_plugins/cdhit:$PATH"
export PATH="/users/a/a/aadas/Bin/TransDecoder-3.0.1:$PATH"
export PATH="/users/a/a/aadas/Bin/hmmer-3.1b2-linux-intel-x86_64/binaries:$PATH"
export PATH="/users/a/a/aadas/Bin/ncbi-blast-2.6.0+/bin:$PATH"

transDecoder_dir=/users/a/a/aadas/Bin/TransDecoder-3.0.1
INPUT_DIR=/users/a/a/aadas/Erharta_data_analysis/transdecoder
cd $INPUT_DIR
######################################################################
###concatenate outputs for blastp and hmmscan searches
#####################################################################
cat blast_out/split.* > blastp.outfmt6
cat hmmscan_out/* > pfam.domtblout
####################################################################################################
####remove files that are no longer needed
###################################################################################################
rm -r split.* blastp-part-* hmmscan-part-*
#####################################################################################################
###submit final step of TransDecoder searching for potential coding regions of the transcripts
#####################################################################################################
$transDecoder_dir/TransDecoder.Predict -t $INPUT_DIR/Trinity.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6
```

***Export the outfmt6 in a csv or text file***

1. Just rename the outfmt6 file with .txt
2. Then open is any text viewer 
3. Copy by command + C. Now open a excel sheet and paste it (make sure you leave one row empty for the headers).
4. Use the header

| query_ID | subject_id | %_identity | alignment length | mismatches | gap opens | q. start | q. end | s. start | s. end | evalue | bit score |
| -------- | ---------- | ---------- | ---------------- | ---------- | --------- | -------- | ------ | -------- | ------ | ------ | --------- |
|          |            |            |                  |            |           |          |        |          |        |        |           |

1. Save the file as csv.

2. Now open it in any text reader (BBEdit), now replace ''::'' with '',''. Now save it and close it. Open the csv now and adjust the header and delete column 2,3,4. Save and replace the csv.

3. Rename both file 1st column with "query_ID"

4. Now use R to merge
```
setwd("~/Desktop/Erharta_data_analysis/BLASTP")
list.files()
data <- read.csv("Erharta_blastp.outfmt6.csv")
data1 <- read.csv("Erharta_drought_freezing_conserved.csv")
head(data)
head(data1)
#merge data by ID
merged_data<- merge(data1,data, by="query_ID")
head(merged_data)
merged_data_final <- as.data.frame(merged_data[,1:14])
head(merged_data_final)
#Delete repeated column
library(dplyr)
data2 <- merged_data_final %>% distinct
#Save as csv file
write.csv(data2, file = "Erharta_Merged_freezing_conserved_annotated.csv")
```
5. Now cut the column 8 and paste it after column 1.

------
<div id='id-section9'/>

## Chapter 9: Go annotation (uniprot)

1. Take the **Merged_freezing_conserved_annotated.xlsx** file and copy the query ID
2. Now go to the uniprot website [http://www.uniprot.org/uploadlists/](http://www.uniprot.org/uploadlists/).
3. Upload your file.
4. Select From: UniProtKB AC/ID, To: UniProtKB (default)
5. Click GO.
6. Click on the “columns” to add all five Gene Ontology columns, KEGG (under 'Genome Annotation'), and PANTHER and Pfam (under 'Family and Domains'). 
7. Click on the “download” button to download your data as a tab separated format.
8. Open in text viewer and copy it in a excel and save as csv.
9. Extract gene names from uniprot iD
```
=LEFT(A2, FIND("_",A2)-1)
```
to count NA in excel

``
=COUNTIF(E2:E150, "NA")
``

------
<div id='id-section10'/>


## Chapter10: Orthofinder

1. Go to https://github.com/davidemms/OrthoFinder/releases and copy link address (OrthoFinder_glibc-2.15.tar.gz) then to **wget** to download to your software directory.

```
[aadas@vacc-user2 Bin]$ tar -zxvf OrthoFinder_glibc-2.15.tar.gz 
[aadas@vacc-user2 OrthoFinder]$ ~/Bin/OrthoFinder/orthofinder -h

OrthoFinder version 2.5.1 Copyright (C) 2014 David Emms

SIMPLE USAGE:
Run full OrthoFinder analysis on FASTA format proteomes in <dir>
  orthofinder [options] -f <dir>

Add new species in <dir1> to previous run in <dir2> and run new analysis
  orthofinder [options] -f <dir1> -b <dir2>

OPTIONS:
 -t <int>        Number of parallel sequence search threads [Default = 4]
 -a <int>        Number of parallel analysis threads
 -d              Input is DNA sequences
 -M <txt>        Method for gene tree inference. Options 'dendroblast' & 'msa'
                 [Default = dendroblast]
 -S <txt>        Sequence search program [Default = diamond]
                 Options: blast, mmseqs, blast_gz, diamond, blast_nucl
 -A <txt>        MSA program, requires '-M msa' [Default = mafft]
                 Options: muscle, mafft
 -T <txt>        Tree inference method, requires '-M msa' [Default = fasttree]
                 Options: iqtree, raxml-ng, fasttree, raxml
 -s <file>       User-specified rooted species tree
 -I <int>        MCL inflation parameter [Default = 1.5]
 -x <file>       Info for outputting results in OrthoXML format
 -p <dir>        Write the temporary pickle files to <dir>
 -1              Only perform one-way sequence search
 -X              Don't add species names to sequence IDs
 -y              Split paralogous clades below root of a HOG into separate HOGs
 -z              Don't trim MSAs (columns>=90% gap, min. alignment length 500)
 -n <txt>        Name to append to the results directory
 -o <txt>        Non-default results directory
 -h              Print this help text

WORKFLOW STOPPING OPTIONS:
 -op             Stop after preparing input files for BLAST
 -og             Stop after inferring orthogroups
 -os             Stop after writing sequence files for orthogroups
                 (requires '-M msa')
 -oa             Stop after inferring alignments for orthogroups
                 (requires '-M msa')
 -ot             Stop after inferring gene trees for orthogroups 

WORKFLOW RESTART COMMANDS:
 -b  <dir>         Start OrthoFinder from pre-computed BLAST results in <dir>
 -fg <dir>         Start OrthoFinder from pre-computed orthogroups in <dir>
 -ft <dir>         Start OrthoFinder from pre-computed gene trees in <dir>

LICENSE:
 Distributed under the GNU General Public License (GPLv3). See License.md

CITATION:
 When publishing work that uses OrthoFinder please cite:
 Emms D.M. & Kelly S. (2019), Genome Biology 20:238

 If you use the species tree in your work then please also cite:
 Emms D.M. & Kelly S. (2017), MBE 34(12): 3267-3278
 Emms D.M. & Kelly S. (2018), bioRxiv https://doi.org/10.1101/267914
```

**Dependencies**
1. **BLAST+**: Go to (https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/). Do  

2. **MCL**: Go https://micans.org/mcl/. **wget** to download to your software directory and do **tar -zxvf** to unzip. 
```
[aadas@vacc-user2 Bin]$ cd mcl-14-137/
[aadas@vacc-user2 mcl-14-137]$ ./configure --prefix=$HOME/local
[aadas@vacc-user2 mcl-14-137]$ make install
```
3. **FastME**: Go http://www.atgc-montpellier.fr/fastme/binaries.php. 
```
[aadas@vacc-user2 Bin]$ cd fastme-2.1.5/
[aadas@vacc-user2 fastme-2.1.5]$ ./configure
[aadas@vacc-user2 fastme-2.1.5]$ make
```
4. **DLCpar**: Go https://www.cs.hmc.edu/~yjw/software/dlcpar/.

Now run orthofinder on the example data

```
[aadas@vacc-user2 OrthoFinder]$ ~/Bin/OrthoFinder/orthofinder -f ExampleData/
```
**FOR YOUR DATA**

1. Make a folder and copy all the **longest_orfs.pep** files of individual species and rename them with **species_name.faa**.
2. Now submit this job either in the server or at your remote option in the folder where all these faa files are stored.

**AT YOUR UNIX**
```
~/Bin/OrthoFinder/orthofinder -f /users/a/a/aadas/Erharta_data_analysis/orthofinder_try
```
**SUBMISSION AT VACC**

```
#!/bin/bash

#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=30:00:00
#SBATCH --mem=256G
#SBATCH --job
#SBATCH --output=Aayudh_orthofinder
#SBATCH --mail-user=aadas@uvm.edu
#SBATCH --mail-type=ALL


export PATH="/users/a/a/aadas/Bin/ncbi-blast-2.6.0+/bin:$PATH"
export PATH="/users/a/a/aadas/Bin/mcl-14-137:$PATH"
export PATH="/users/a/a/aadas/Bin/fastme-2.1.5/bin:$PATH"
export PATH="/users/a/a/aadas/Bin/dlcpar-2.0.1/bin:$PATH"
export PATH="/users/a/a/aadas/Bin/OrthoFinder:$PATH"

OrthoFinder_dir=/users/a/a/aadas/Bin/OrthoFinder
INPUT_DIR=/users/a/a/aadas/Erharta_data_analysis/orthofinder_try
cd $INPUT_DIR

$OrthoFinder_dir/orthofinder -f /users/a/a/aadas/Erharta_data_analysis/orthofinder_try
```

**SRA Toolkit**
Getting data from SRA (https://www.ncbi.nlm.nih.gov/sra). 

```
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR6127933/SRR6127933.1
module load sratoolkit-2.9.6-gcc-7.3.0-65lpczt
fastq-dump --split-files --origfmt --gzip SRR6127933.1
```
or after wget just run a script

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

workDIR=/users/a/a/aadas/2nd_chapter_RNASEQ/nassella_trimming/6wkcold_nassella_rep2

module load sratoolkit-2.9.6-gcc-7.3.0-65lpczt
ulimit -s unlimited

fastq-dump --split-files --origfmt --gzip $workDIR/SRR6127937.1
```
------
<div id='id-section11'/>

## Chapter11: Enrichment analysis

1. Paste uniprot ID column (format: Q9SJW9) in http://geneontology.org/
2. Use Biological process and Arabidopsis as organism.
3. launch and download the table.
4. To capitalize the first letter of the BIOlogical processes use =UPPER(LEFT(A2,1))&RIGHT(A2,LEN(A2)-1)
5. Then use the R code-

```
setwd("~/OneDrive - University of Vermont/Freezing tolerance/Nassella Freezing/Enrichment_nassella freezing")
list.files()
data <- read.csv("Nassella_freezing_enrichment.csv")
head(data)
library("ggplot2")


tiff("Nassella_freezing_enrichment.tiff", width = 9, height = 7.5, units = 'in', res = 300)
ggplot(data, aes(x =fold.Enrichment, y = GO.biological.process.complete, 
                 size = P.value, fill = FDR)) +
  geom_point(shape = 21) +
  labs(x = "Fold Enrichment", y = "Biological Process (Nassella_freezing)",
       size = "pvalue", fill = "FDR")+
  scale_fill_gradient(low="#CC6666", high="#9999CC")+
  theme_classic() +
  theme(legend.position="right",text=element_text(size=18, colour = "black", family="Times"),
        axis.line = element_line(size=0.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 18),
        axis.text.y=element_text(colour="black", size = 10),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()
```
------
<div id='id-section12'/>

## Chapter12: KEGG analysis

1. KEGG ID conversion: GO to https://www.genome.jp/kegg/tool/conv_id.html. Then use 20 inputs of uniport id (eg: D9J0Z7). Select **Outside DB: Uniport, organism ath**, then exec. Input organism code if it's a different one.
2. If you are using uniprot based data sheet then from KEGG IDs remove the ':' in the end by using =LEFT(A2, LEN(A2) - 1).
3. Now go to https://www.genome.jp/kegg/mapper/search.html. Use organism ath for Arabidopsis in **other org** box. paste codes (eg-ath:AT3G28610)
and exec.

------
<div id='id-section13'/>

## Chapter13: NGS FAQ

### cDNA library prep:
To perform RNA sequencing (RNA-seq), a library preparation is required. Library preparation involves converting the RNA molecules from a sample into a library of cDNA (complementary DNA) fragments that can be sequenced. 

#### Steps
To create a cDNA (complementary DNA) library, several steps are involved. The cDNA library is a collection of cloned cDNA fragments that represent the expressed genes of an organism [4]. Here are the general steps to make a cDNA library:

Isolation of mRNA: The first step is to isolate cellular mRNA from the sample. This mRNA extract should represent all of the transcripts in the cells at the time of isolation, known as the transcriptome [1].

Reverse transcription: The isolated mRNA is then reverse transcribed into cDNA using reverse transcriptase enzyme and oligo(dT) primers. Reverse transcription converts the mRNA into complementary DNA (cDNA) molecules, which represent the genetic information of the expressed genes [2].

Second-strand synthesis: After reverse transcription, the single-stranded cDNA is converted into double-stranded cDNA using DNA polymerase and DNA ligase. This step involves the synthesis of the complementary strand to the cDNA template [2].

Cloning: The double-stranded cDNA fragments are then ligated into a suitable cloning vector, such as a plasmid or a viral vector. This step allows for the incorporation of the cDNA fragments into host cells, forming a library [4].

Transformation: The ligated cDNA vectors are introduced into host cells, typically bacteria such as Escherichia coli. This process is known as transformation, where the host cells take up the vectors carrying the cDNA fragments [6].

Screening and validation: The final step involves screening and validating the cDNA library to ensure the presence of the desired cDNA fragments and correct plasmids. Various screening methods, such as transformation and replica plating, are used to select transformed cells carrying the cDNA of interest [6].

#### Adapters used 

The adapter commonly used in Illumina HiSeq sequencing platforms is known as the **"TruSeq Adapter."** The TruSeq Adapter is designed to facilitate the attachment of sequencing primers to the DNA fragments during library preparation. It contains specific sequences required for proper binding and amplification during the sequencing process. The TruSeq Adapter is compatible with various Illumina library preparation kits and is widely used in Illumina HiSeq sequencing workflows.

#### How to calculate the coverage per nucleotide in RNA-seq?

Obtain the total number of reads: Determine the total number of reads generated from your RNA-seq experiment. This information can usually be obtained from the sequencing data files.

Align reads to the reference genome/transcriptome: Use a bioinformatics tool, such as STAR, HISAT2, or Bowtie, to align the sequencing reads to a reference genome or transcriptome. This step assigns each read to a specific genomic location.

Count the number of reads per nucleotide: Once the reads are aligned, you can count the number of reads that align to each nucleotide position in the reference. This can be done using tools like SAMtools or BEDtools. The resulting output will provide read counts per nucleotide.

Normalize the counts: Normalize the read counts by the total number of aligned reads and the length of the reference sequence. This step accounts for differences in sequencing depth and allows for comparison between samples. The normalized value represents the coverage per nucleotide.

The formula to calculate coverage per nucleotide is:

**Coverage per nucleotide = (Number of reads aligned to a specific nucleotide) / (Total number of aligned reads * Reference sequence length)**

How to Estimate and Achieve Your Desired NGS Coverage Level?
Estimate Sequencing Runs:

The Lander/Waterman equation1 is a method for computing genome coverage. The general equation is: **C (coverage) = L(read length) * N (the number of reads) / G (haploid genome length)**

### RNA-seq quality check

Performing quality checks on RNA-seq data is crucial to ensure reliable and accurate results. Here are some commonly used methods for RNA-seq quality check:

FastQC: FastQC is a widely used tool for assessing the quality of sequencing data. It provides a comprehensive analysis of various quality metrics such as per-base sequence quality, sequence length distribution, GC content, overrepresented sequences, adapter content, and more. FastQC generates detailed reports and visualizations to identify potential issues in the data.

Read Distribution: Examining the distribution of reads across different genomic features can provide insights into the quality of RNA-seq data. Tools like RNA-SeQC or RSeQC can be used to assess the distribution of reads across exons, introns, intergenic regions, ribosomal RNA (rRNA), and other genomic features.

Mapping Statistics: Aligning RNA-seq reads to a reference genome or transcriptome can provide information about mapping efficiency and potential biases. Tools such as STAR, HISAT2, or Bowtie can be used to generate mapping statistics, including alignment rates, multi-mapping rates, and strand specificity.

Gene Body Coverage: Analyzing the coverage of reads across gene bodies can assess the uniformity of coverage and potential biases. Tools like RSeQC or QualiMap can generate gene body coverage plots and metrics to evaluate the evenness of coverage.

Transcript Integrity: Assessing transcript integrity helps to identify potential degradation or bias in the RNA-seq data. Tools like RSeQC or RNA-SeQC can estimate transcript integrity based on metrics such as 3' to 5' bias, coverage uniformity, or junction saturation analysis.

Reproducibility Analysis: If you have replicate samples, assessing the reproducibility between them is important. Tools like MultiQC or DESeq2 can be used to generate correlation plots, Principal Component Analysis (PCA), or perform differential gene expression analysis to evaluate the consistency between replicates.

It's important to note that the choice of tools and methods may vary depending on the specific requirements and analysis goals of the RNA-seq experiment.
