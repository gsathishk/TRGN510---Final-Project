# TRGN510---Final-Project

## 1. Transfer files from local machine to HPC Server TRGN510 Project Folder
```
rsync -a -P ~/TRGN510project sganesan@hpc-transfer.usc.edu:~/project/TRGN510project

```
## 2. Unzip transferred files using gunzip
###   gunzip *.tar.gz &

## 3. Converted filename.fastq.tar to filename.fastq
### tar -xf filename.fastq.tar

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
  
 
