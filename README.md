# TRGN510---Final-Project

# RNAseq Analysis of CTC Bone Metastasis Samples

## Description:
   I will be processing RNA-seq data from Brx68 Parentals and Brx68 Bone Metastasis samples using TopHat. Brx68 is an ex vivo cultured circulating tumor cell (CTC) cell line. Previously we had processed the same samples using STAR, so I will be comparing results of TopHAT to STAR. Outside of the class I will do these processes: In addition I will also be processing samples Brx50 Parentals and Brx50 Bone Metastasis. I will then compare TopHAT results of  Brx50 and Brx68 analysis to find commonly upregulated genes in bone metastasis.     

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

## 12. Trimming fastq files
       Before using the fastq files for alignment via tophat, I trimmed using Trim Galore:
       
       First I downloaded Trim Galore!
       Website: https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

```
wget https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.zip

```

     Then unzipped and installed it.
   
     Trim Galore! needs 'cutadapt' to function properly, so I downloaded that next.
     Website: https://cutadapt.readthedocs.io/en/stable/installation.html#quick-installation

```
wget https://pypi.python.org/packages/68/73/2ae48245bbf6d84a24bdf29540ed01669f09c5d21c26258f2ce07e13c767/cutadapt-1.16.tar.gz#md5=222d062207e2a3a2a0760caa83bf14ff

```
    To install cutadapt, I need pip, so i downloaded that next. 
    
```
wget https://bootstrap.pypa.io/get-pip.py

```
    Next used this command as per installation guide :

```
python get-pip.py
```
    But this gave an error and so I tried this:
```
python get-pip.py --user
```
    pip was then successfully installed.
    Next continued with installign cutadapt.
```
pip install --user --upgrade cutadapt
```
    cutadapt was installed successfully!
    I moved cutadapt, trimgalore to bin folder
    Next, ran trimgalore on one file to check, it worked without any errors, so used the following slurm scirpt to batch run trimgalore on all fasta files
```
#!/bin/bash

#SBATCH --ntasks=16

#SBATCH --time=18:00:00

#SBATCH --mail-type=ALL

#SBATCH --mail-user=sganesan@usc.edu

cd ~/project/TRGN510project/TRGN510project
~/project/bin/trim_galore *.fastq

```
   A portion of output summary for one of the processed files is shown below:

```
Writing report to '68H.R3.fastq_trimming_report.txt'

SUMMARISING RUN PARAMETERS
==========================
Input filename: 68H.R3.fastq
Trimming mode: single-end
Trim Galore version: 0.4.4_dev
Cutadapt version: 1.16
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; auto-detected)
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length before a sequence gets removed: 20 bp

Writing final adapter and quality trimmed output to 68H.R3_trimmed.fq


  >>> Now performing quality (cutoff 20) and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file 68H.R3.fastq <<<
10000000 sequences processed
20000000 sequences processed
30000000 sequences processed
This is cutadapt 1.16 with Python 2.7.5
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC 68H.R3.fastq
Running on 1 core
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 915.04 s (27 us/read; 2.26 M reads/minute).

=== Summary ===

Total reads processed:              34,453,893
Reads with adapters:                 3,635,993 (10.6%)
Reads written (passing filters):    34,453,893 (100.0%)

Total basepairs processed: 2,602,788,033 bp
Quality-trimmed:               6,993,547 bp (0.3%)
Total written (filtered):  2,591,058,972 bp (99.5%)

=== Adapter 1 ===

Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 3635993 times.

No. of allowed errors:
0-9 bp: 0; 10-13 bp: 1

Bases preceding removed adapters:
  A: 35.0%
  C: 27.2%
  G: 19.3%
  T: 18.5%
  none/other: 0.0%

Overview of removed sequences
length	count	expect	max.err	error counts
1	2971303	8613473.2	0	2971303
2	377417	2153368.3	0	377417
3	236432	538342.1	0	236432
4	36067	134585.5	0	36067
5	10322	33646.4	0	10322

```
	 
## 13. Alignment using TopHat2 and Bowtie2.
       I didn't use the GENCODE annotation and sequence files that I downloaded earlier.
       Tophat2 provides already annotation and bowtie2 indexed files for GRCh27/hg19 so I downloaded those
       Website: https://ccb.jhu.edu/software/tophat/igenomes.shtml
```
wget  ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz
```
      Next, I ran tophat on one of my files to check if it works:
```
	To do tophat on every trimmed file in the folder:
	
	#!/bin/bash
	
	#SBATCH --ntasks=24
	
	#SBATCH --time=24:00:00
	
	#SBATCH --mail-type=ALL
	
	#SBATCH --mail-user=sganesan@usc.edu
	
	cd ~/project/TRGN510project/TRGN510project/TrimmedFiles
        ~/project/bin/tophat2 --GTF ~/project/TRGN510project/TRGN510project/TrimmedFiles/genes.gtf --no-novel-juncs genome 68bone-21-1_trimmed.fq
```
       When I submitted the job by sbatch, It couldn't complete the alignment even in 5 hours, so I add 'tophat2 -p 8' option to make it faster.
       Now, the alignment took only around 2 hours.
       Slurm output is shown below:
```
[2018-04-13 16:49:19] Beginning TopHat run (v2.1.1)
-----------------------------------------------
[2018-04-13 16:49:19] Checking for Bowtie
		  Bowtie version:	 2.3.4.1
[2018-04-13 16:49:20] Checking for Bowtie index files (genome)..
[2018-04-13 16:49:20] Checking for reference FASTA file
[2018-04-13 16:49:20] Generating SAM header for genome
[2018-04-13 16:49:49] Reading known junctions from GTF file
[2018-04-13 16:50:09] Preparing reads
	 left reads: min. length=20, max. length=76, 30076813 kept reads (4953 discarded)
[2018-04-13 16:56:49] Building transcriptome data files ./tophat_out/tmp/genes
[2018-04-13 16:57:48] Building Bowtie index from genes.fa
[2018-04-13 17:16:02] Mapping left_kept_reads to transcriptome genes with Bowtie2
[2018-04-13 17:52:12] Resuming TopHat pipeline with unmapped reads
[2018-04-13 17:52:13] Mapping left_kept_reads.m2g_um to genome genome with Bowtie2
[2018-04-13 18:19:38] Mapping left_kept_reads.m2g_um_seg1 to genome genome with Bowtie2 (1/3)
[2018-04-13 18:21:45] Mapping left_kept_reads.m2g_um_seg2 to genome genome with Bowtie2 (2/3)
[2018-04-13 18:23:46] Mapping left_kept_reads.m2g_um_seg3 to genome genome with Bowtie2 (3/3)
[2018-04-13 18:26:04] Retrieving sequences for splices
[2018-04-13 18:27:11] Indexing splices
Building a SMALL index
[2018-04-13 18:27:27] Mapping left_kept_reads.m2g_um_seg1 to genome segment_juncs with Bowtie2 (1/3)
[2018-04-13 18:27:44] Mapping left_kept_reads.m2g_um_seg2 to genome segment_juncs with Bowtie2 (2/3)
[2018-04-13 18:28:03] Mapping left_kept_reads.m2g_um_seg3 to genome segment_juncs with Bowtie2 (3/3)
[2018-04-13 18:28:24] Joining segment hits
[2018-04-13 18:30:22] Reporting output tracks
-----------------------------------------------
[2018-04-13 18:55:21] A summary of the alignment counts can be found in ./tophat_out/align_summary.txt
[2018-04-13 18:55:21] Run complete: 02:06:01 elapsed
```
       Then I modified the script to do the following:
       a) Use a for loop to process all files ending in .fq
       b) Include the name of the input file to the tophat output file so I can distinguish output from all these different trimmed .fq files
       c) use move command to rename accepted_hits.bam to filename.bam
```
      	#!/bin/bash
	
	#SBATCH --ntasks=24
	
	#SBATCH --time=24:00:00
	
	#SBATCH --mail-type=ALL
	
	#SBATCH --mail-user=sganesan@usc.edu
	
	cd ~/project/TRGN510project/TRGN510project/TrimmedFiles
	
	for file in *.fq; do
	    samplename=$(basename $file.fq)
	    ~/project/bin/tophat2 -p 8 -o $samplename-tophat --GTF ~/project/TRGN510project/TRGN510project/TrimmedFiles/genes.gtf --no-novel-juncs genome $file; mv $samplename-tophat/accepted_hits.bam $samplename-tophat/$samplename.bam
	    done
```   
       After 24 hours, remaining files were  processed by re-running the script.
       All files were successfully processed.

## 14. Analysis using CuffDiff
       I downloaded CuffDiff to perform differential analysis from mapped reads.
       Website: http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/

``` 
wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz

```

      Then I ran CuffDiff on my unit test data using the following script:
```
#!/bin/bash

#SBATCH --ntasks=24

#SBATCH --time=24:00:00

#SBATCH --mail-type=ALL

#SBATCH --mail-user=sganesan@usc.edu

cd ~/project/TRGN510project/TRGN510project/TrimmedFiles/CompletedTophat/

~/project/bin/cuffdiff -o ~/project/TRGN510project/TRGN510project/TrimmedFiles/CompletedTophat/50BoneVsParental_TEST -L 50Bone,50Parental -p 4 --frag-bias-correct ~/project/TRGN510project/TRGN510project/TrimmedFiles/genome.fa --multi-read-correct --library-norm-method quartile ~/project/TRGN510project/TRGN510project/TrimmedFiles/genes.gtf 50Bone 50Parental         
```

     A sample of output file gene_exp.diff from the test analysis comparing 50Bone and parental is shown here:
```
sganesan@hpc-login3:~/project/TRGN510project/TRGN510project/TrimmedFiles/CompletedTophat/50BoneVsParental_TEST/50ParentalVsBone_TEST cat gene_exp.diff | head
test_id	gene_id	gene	locus	sample_1	sample_2	status	value_1	value_2	log2(fold_change)	test_stat	p_value	q_value	significant
ENSG00000000003	ENSG00000000003	TSPAN6	X:99883666-99894988	50Bone.bam	50Parental.bam	OK	10.6042	2.07391	-2.35421	-2.06356	0.042	0.994919	no
ENSG00000000005	ENSG00000000005	TNMD	X:99839798-99854882	50Bone.bam	50Parental.bam	NOTEST	0	0	0	0	1	1	no
ENSG00000000419	ENSG00000000419	DPM1	20:49505584-49575092	50Bone.bam	50Parental.bam	OK	54.0528	48.7805	-0.148064	-0.0653279	0.9213	0.994919	no
ENSG00000000457	ENSG00000000457	SCYL3	1:169631244-169863408	50Bone.bam	50Parental.bam	OK	5.9338	12.042	1.02104	0.816778	0.41145	0.994919	no
ENSG00000000460	ENSG00000000460	C1orf112	1:169631244-169863408	50Bone.bam	50Parental.bam	OK	8.43574	2.41329	-1.80551	-0.732915	0.1979	0.994919	no
ENSG00000000938	ENSG00000000938	FGR	1:27938574-27961788	50Bone.bam	50Parental.bam	NOTEST	0	0.0431237	inf	0	1	1	no
ENSG00000000971	ENSG00000000971	CFH	1:196621007-196716634	50Bone.bam	50Parental.bam	NOTEST	0	0.0607381	inf	0	1	1	no
ENSG00000001036	ENSG00000001036	FUCA2	6:143771943-143832827	50Bone.bam	50Parental.bam	OK	52.3181	25.7441	-1.02307	-0.844985	0.32195	0.994919	no
ENSG00000001084	ENSG00000001084	GCLC	6:53362138-53481969	50Bone.bam	50Parental.bam	OK	34.3053	27.3201	-0.328469	-0.262155	0.71015	0.994919	no
```
 I used rsync to transfer this file to local machine and then imported into R for analysis.
CuffDiff analysis will be run on Brx68Parental and Bone samples

## 15. Data Manipulation and Visualization using R 

      After transferring data to local machine, I removed those fields with FPKM values of 0 to avoid genes that show 'infinite' as log2fold change value.
      Then I imported the data into R and performed the following analysis:
      1. First I created a newdata frame by subsetting only 'gene', 'value_1', 'value-2', 'log2fold change' and 'p-value'.
         Values 1 and 2 refer to coloumns with FPKM values for bone and parental respectively.
      2. Next, I created a volcano plot to visualize changes in gene experssion of bone vs parental samples. Log2fold change and p-values are used for plotting the data.
         
```
getwd()
setwd('/Users/sathishkumarganesan/TRGN510project/68BoneVsParentals')
#Extract Table from CuffDiff analysis into R.
Brx68CuffDiff <- read.table("gene_exp.diff", header = TRUE)
#Making a new table with 'gene', 'value_1', 'value-2', 'log2fold change' and 'p-value'.
GeneExpression <- Brx68CuffDiff[c(1:63657),c(3,8:10,12)]
#Removing genes with FPKM values less than 0.9 in Brx68 parental
GeneExpression_Value1 <- subset(GeneExpression, value_1 > 0.9)
# From that subset, removing genes with FPKM values less than 0.9 in Brx68 Bone 
Brx68ParentalvsBone <- subset(GeneExpression_Value1, value_2 > 0.9)
## To make Volcano Plot 
with(Brx68ParentalvsBone, plot(log2.fold_change., -log10(p_value), pch=20, col="grey", main="Brx68 Parental vs Bone", xlim=c(-7,7), ylim=c(0,4)))
#To colour genes with pvalue less than 0.05 and log2fold change in expression more than 1.5
with(subset(Brx68ParentalvsBone, p_value<.05 & abs(log2.fold_change.)>1.5), points(log2.fold_change., -log10(p_value), pch=20, col="darkcyan"))
# To add text to genes that have pvalue less than 0.01 and log2fold change greater than 4.0
library(calibrate)
with(subset(Brx68ParentalvsBone, p_value<.01 & abs(log2.fold_change.)>2.5), textxy(log2.fold_change., -log10(p_value), labs=gene, cex=.3))
```
## 16. R Shiny and Deployment

      I created a new app.R document with the following code to generate plot of log2Foldchange vs -log10pvalue. 
      Added a radio button to choose colours for the plot
      Added option to click on a datapoint which would then retreive all the row data (gene name, FPKM values, log2foldchange, p-value and -log10pvalue)
```
library(shiny)
ParentalvsBoneMets <- read.table("Data.txt", header = TRUE)

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Gene Expression Changes in Bone Metastatic CTCs"),
  # Copy the line below to make a set of radio buttons
  radioButtons("radio", label = h5("Plot Color"),
               choices = list("Blue", "Orange", "Black"), 
               selected = "Blue"),
  
 
  plotOutput("plot1", click = "plot_click"),
  verbatimTextOutput("info")
  
)

server <- function(input, output) {
  
  # You can access the values of the widget (as a vector)
  # with input$radio, e.g.
  output$plot1 <- renderPlot({
    plot(ParentalvsBoneMets$log2.fold_change., ParentalvsBoneMets$log10pvalue, col=input$radio)
  })
  
  output$info <- renderPrint({
    # With base graphics, need to tell it what the x and y variables are.
    nearPoints(ParentalvsBoneMets, input$plot_click, xvar = "log2.fold_change.", yvar = "log10pvalue")
    # nearPoints() also works with hover and dblclick events
  })
}
  
# Run the application 
shinyApp(ui = ui, server = server)
```

     Next, to deploy:
     Used rsync to copy the app.R and the dataset into server shiny folder to deploy. Verified by checking on the browser.

    Volcano Plot from 15 has been uploaded as image in Github.
    R shiny deployment has been uploaded as PDF showing the server name and the deployed plot.
