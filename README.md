# TRGN510---Final-Project

# RNAseq Analysis of CTC Bone Metastasis Samples

## Description:
   I will be processing RNA-seq data from Brx68 Parentals and Brx68 Bone Metastasis samples using TopHat. Brx68 is an ex vivo cultured circulating tumor cell (CTC) cell line. Previously we had processed the same samples using STAR, so I will be comparing results of TopHAT to STAR. In addition I will also be processing samples Brx50 Parentals and Brx50 Bone Metastasis. I will then compare TopHAT results of  Brx50 and Brx68 analysis to find commonly upregulated genes in bone metastasis.     

## 1. Transfer files from local machine to HPC Server TRGN510 Project Folder
```
rsync -a -P ~/TRGN510project sganesan@hpc-transfer.usc.edu:~/project/TRGN510project

```
## 2. Unzip transferred files using gunzip
```
gunzip *.tar.gz &

```
## 3. Converted filename.fastq.tar to filename.fastq
```
tar -xf filename.fastq.tar

```
## 4. The previous step was to test the untar conversion. But, to untar mutiple files at the same time:
```  
vim untar.sh

```

```
#!/bin/bash
#loop to untar all files
for f in *.tar; do
 tar xf $f &
done
wait
```

## 5. I downloaded FASTQC from www.bioinformatics.babraham.ac.uk
###   Then performed these set of commands to download, unzip and make it executable

```
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc_v0.11.7.zip
chmod 755 fastqc

```

## 6. I ran a sample analysis using fastqc
###   For this, I CD'ed into the folder where my fastq files reside and then ran this command:

```
~/project/NGSTools/FastQC/fastqc 68H.fastq

```
#### This generated a 68H_fastqc.zip, which I unzipped it. There were summary fastqc, data fastqc and other files
#### When I cat the summary fastqc file, it displayed this:

	PASS	Basic Statistics	68H.fastq
	PASS	Per base sequence quality	68H.fastq
	PASS	Per tile sequence quality	68H.fastq
	PASS	Per sequence quality scores	68H.fastq
	FAIL	Per base sequence content	68H.fastq
	WARN	Per sequence GC content	68H.fastq
	PASS	Per base N content	68H.fastq
	WARN	Sequence Length Distribution	68H.fastq
	WARN	Sequence Duplication Levels	68H.fastq
	WARN	Overrepresented sequences	68H.fastq
	PASS	Adapter Content	68H.fastq

## 7. Next a created a git hub repo called TRGN510---Final-Project then went to command and performed following commands to create README.md file and push it

```
echo "#TRGN510---Final-Project" >> README.md
git init
git add README.md
git commit -m "test commit"
git remote add origin https://github.com/gsathishk/TRGN510---Final-project.git
git push -u origin master

```
  
## 8. Downloaded GRCh37/hg19 reference genome and annotation from GENCODE
     Website: https://www.gencodegenes.org/releases/28lift37.html

```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.transcripts.fa.gz

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.annotation.gtf.gz

```

## 9. Dowloaded TopHat 2.1.1
     Website : https://ccb.jhu.edu/software/tophat/index.shtml

```
wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz

```	

## 10. Downloaded Bowtie2
       Website : https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1

```
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip

```

## 11. Testing TopHat2 and Bowtie2
       Test data from : https://ccb.jhu.edu/software/tophat/tutorial.shtml
 
```
wget  https://ccb.jhu.edu/software/tophat/downloads/test_data.tar.gz

```
       Ran the following command to test:

```
tophat2 -r 20 test_ref reads_1.fq reads_2.fq

```
      Shown below is a part of the output:

```
[2018-04-11 22:34:49] Beginning TopHat run (v2.1.1)
-----------------------------------------------
[2018-04-11 22:34:49] Checking for Bowtie
		  Bowtie version:	 2.3.4.1
[2018-04-11 22:34:49] Checking for Bowtie index files (genome)..
	Found both Bowtie1 and Bowtie2 indexes.
[2018-04-11 22:34:49] Checking for reference FASTA file
[2018-04-11 22:34:49] Generating SAM header for test_ref
[2018-04-11 22:34:49] Preparing reads
	 left reads: min. length=75, max. length=75, 100 kept reads (0 discarded)
	right reads: min. length=75, max. length=75, 100 kept reads (0 discarded)
[2018-04-11 22:34:49] Mapping left_kept_reads to genome test_ref with Bowtie2
[2018-04-11 22:34:49] Mapping left_kept_reads_seg1 to genome test_ref with Bowtie2 (1/3)
[2018-04-11 22:34:49] Mapping left_kept_reads_seg2 to genome test_ref with Bowtie2 (2/3)
[2018-04-11 22:34:49] Mapping left_kept_reads_seg3 to genome test_ref with Bowtie2 (3/3)
[2018-04-11 22:34:49] Mapping right_kept_reads to genome test_ref with Bowtie2
[2018-04-11 22:34:50] Mapping right_kept_reads_seg1 to genome test_ref with Bowtie2 (1/3)
……..

[2018-04-11 22:34:51] Mapping right_kept_reads_seg1 to genome segment_juncs with Bowtie2 (1/3)
[2018-04-11 22:34:51] Mapping right_kept_reads_seg2 to genome segment_juncs with Bowtie2 (2/3)
[2018-04-11 22:34:51] Mapping right_kept_reads_seg3 to genome segment_juncs with Bowtie2 (3/3)
[2018-04-11 22:34:52] Joining segment hits
[2018-04-11 22:34:52] Reporting output tracks
-----------------------------------------------
[2018-04-11 22:34:52] A summary of the alignment counts can be found in ./tophat_out/align_summary.txt
[2018-04-11 22:34:52] Run complete: 00:00:03 elapsed

```
       It created a tophat_out directory, as expected. An example of the expected results from junctions.bed is shown here:

```
track name=junctions description="TopHat junctions"
test_chromosome 179     400     JUNC00000001    37      +       179     400     255,0,0 2       71,50   0,171
test_chromosome 350     550     JUNC00000002    37      +       350     550     255,0,0 2       50,50   0,150

```
      The test output confirms TopHat2 and Bowtie2 were installed properly and function as expected

 
