# Erharta

Login info: ssh aadas@vacc-user2.uvm.edu

#### Where are our Erharta files?

```
cd Erharta/Preston_JP_iLabs_14958_PR75_090820/
```
C= Control, D=Drought, F=Freezing

#### Copying R1 and R2 files to a folder?

First you should in the folder from where you will copy then direct which folder it is going

```
cp ECC1_Illumina13_S13_L002_R1_001.fastq.gz ~/Erharta_data_analysis/Erharta_trimming/
```
 Now to create a blank script

```
[aadas@bluemoon-user2 ~]$ vi trimmomatic.sh
```
To reopen **vim** filename.sh
Now press **i** to insert 

Copy and paste everything present in the new script from this sample one

1) change workDIR= navigate to the the folder where you have separated the R1 and R2. Now use 'pwd' command.

use workDIR=/users/a/a/aadas/Erharta_data_analysis/Erharta_trimming/Control_1


```
#!/bin/bash

######## This job needs 1 nodes, 2 processors total
#PBS -q poolmemq
#PBS -l nodes=1:ppn=2,mem=16gb,vmem=18gb
# it needs to run for 6 hours
#PBS -l walltime=06:00:00
#PBS -N job
#PBS -j oe
#PBS -M aadas@uvm.edu
#PBS -m bea
###LOAD JAVA MODULE AVAILABLE FROM THE CLUSTER, YOU MAY WANT TO CHECK FIRST
module load java-sdk/sun-jdk-1.6.0.12
ulimit -s unlimited
###CHANGE THE DIRECTORY ACCORDINGLY, THE FOLLOWING SETTINGS ARE FOR MY ACCOUNT
SOFTWARE=/users/a/a/aadas/Trimmomatic-0.36
workDIR=/users/a/a/aadas/Erharta_data_analysis/Erharta_trimming/Control_1
cd $workDIR
#####TRIMMING COMMANDS AND PARAMETERS
java -jar $SOFTWARE/trimmomatic-0.36.jar PE -phred33 $workDIR/ECC1_Illumina13_S13_L002_R1_001.fastq.gz $workDIR/ECC1_Illumina13_S13_L002_R2_001.fastq.gz $workDIR/ECC1_R1_paired.fq.gz $workDIR/ECC1.R1.unpaired.fq.gz $workDIR/ECC1_R2.paired.fq.gz $workDIR/ECC1.R2.unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

and **esc**, **:wq** to save and quit. 

**:q!** only quit without saving.

look how I edited the things
