# Erharta Transcriptome pipeline 

Login info: ssh aadas@vacc-user2.uvm.edu

#### Where are our Erharta files?

```
cd Erharta/Preston_JP_iLabs_14958_PR75_090820/
```
C= Control, D=Drought, F=Freezing

## Trimming

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

## Concatenation

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


## Assembly

Move concatenated R1 and R2 trimmed files to a different folder

```
[aadas@vacc-user2 Concatenation_Erharta]$ cp Erharta.R1.trimmed.fq /users/a/a/aadas/Erharta_data_analysis/assembly_erharta
[aadas@vacc-user2 Concatenation_Erharta]$ cp Erharta.R2.trimmed.fq /users/a/a/aadas/Erharta_data_analysis/assembly_erharta
```

