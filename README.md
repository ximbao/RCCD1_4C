# RCCD1_4C

Submission examples and codes

4C Overall Pipeline (Cis Interactions)

- Required to perform data pre-processing - **forward and reverse primer, index seq, cutter names**
- **Required Softwares and tools:** fastqc, cutadapt, trimmomatic, bwa mem, samtools
- **First** - **conda activate 4C_analysis** (Activating conda environment in the terminal where the required softwares are installed)

Before doing the alignment of the fastq sequences, it is necessary to remove the PCR primers and Illumina adapters that is included in some of the sequences of the 4C sequence libraries.

**Software** - Cutadapt

**Note:** Please run your .fastq files using FastQC to check the quality of your sequences. Run FastQC before cutadapt and after doing cutadapt.

```bash
cutadapt -g ^SEQUENCE_OF_THE_FORWARD_PRIMER -o /path/to/the/desired/folder/output.fq /path/to/the/source/input.fq

cutadapt -g ^SEQUENCE_OF_THE_REVERSE_PRIMER -o /path/to/the/desired/folder/output.fq /path/to/the/source/input.fq

cutadapt -a SEQUENCE_OF_ILLUMINA_INDEX -o /path/to/the/desired/folder/output.fq /path/to/the/source/input.fq
```

**Note**: Please note that the first input file is the original .fq file. After trimming the forward primer and obtaining the output file, **that resulting output file will be the input file for the next step, which is trimming of the REVERSE primer**.

After removing the PCR primer sequences and the Illumina adapters (indeces), depending on the FastQC report, a quality filtering may be ideal for removing bad quality reads from the libraries. This can be done using Trimmomatic and using the following command on the following command line.

**IMPORTANT NOTE:** In the command line, please make sure you are in the Trimmomatic directory when using Trimmomatic before using the program.

```bash
java -jar trimmomatic-0.39.jar SE -phred33 /path/to/the/input/file.fq /path/to/the/desired/output/filename.fq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:54
```
Here is a brief description of Trimmomatic and the settings from the above command. For a more detailed version, please refer to the link provided above.

Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data.The selection of trimming steps and their associated parameters are supplied on the command line.

The current trimming steps are:

- ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
- SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
- LEADING: Cut bases off the start of a read, if below a threshold quality
- TRAILING: Cut bases off the end of a read, if below a threshold quality CROP: Cut the read to a specified length HEADCROP: Cut the specified number of bases from the start of the read
- MINLEN: Drop the read if it is below a specified length

After Trimming: Trimmomatic will also tell you how many sequences "survived" and how many were "dropped"

Input Reads: 2677386 Surviving: 2389762 (89.26%) Dropped: 287624 (10.74%)


After quality filtering and checking the quality using the **FastQC**, and if the sequence the sequences are in good condition, the alignment can proceed. For the reference genome, the **hg19** build will be used. Please download the reference genome **hg19** from the UCSC website.

**Note:** Indexing can be done only once and it creates set of indexed hg19 files such as hg19.fa

**Indexing the reference genome**

Before the alignment, hg19 has to be indexed first using the following command:

```bash
./bwa mem index.  Or using  **bwa index -a bwtsw hg19.fa**
```

**Aligning files to the indexed reference genome**

```bash
bwa mem /path/to/the/folder/with/indexed/genome.fa /path/to/the/input/processed/file.fq > /path/to/the/desired/folder/file/aligned_file.sam
```

**Conversion of the .sam output of alignment to .bam**

The .sam output from the BWA alignment is being converted to .bam because: it significantly reduces the amount of the output takes from the memory due to its binary nature. Additionally, and more importantly, the main 4C analysis using r3Cseq, requires the input file to be  in .bam form

```bash
samtools view -S -b input.sam > output.bam
```

**Looking at the statistics of the alignment**

```bash
samtools flagstat output.bam
```

**##The report would appear like this**

3383484 + 0 in total (QC-passed reads + QC-failed reads)

0 + 0 secondary

22068 + 0 supplementary

0 + 0 duplicates

3350647 + 0 mapped (99.03% : N/A)

0 + 0 paired in sequencing

0 + 0 read1

0 + 0 read2

0 + 0 properly paired (N/A : N/A)

0 + 0 with itself and mate mapped

0 + 0 singletons (N/A : N/A)

0 + 0 with mate mapped to a different chr

0 + 0 with mate mapped to a different chr (mapQ>=5)

**Sorting the .bam file**

```bash
samtools flagstat input.bam > output.sorted.bam
```

Fourth Step: Using r3Cseq package in R

For the r3Cseq, it is recommended to be using the latest version of the bioconductor package, the 3Cseq package itself, as well as RStudio to make the analysis easier.

Before starting, make sure that RStudio/R is within the desired working directory. Alternatively, setting up a new project with the option of creating a new directory, will also work.

Also do not forget to have the .bam files (or .sorted.bam files) to be within the desired working directory.

**The following are needed for using r3Cseq:**

- Forward primer sequence
- Reverse primer sequence
- Restriction Enzyme Name
- The chromosome of the Viewpoint (chromosome in which the bait is located)

**r3Cseq can be downloaded and installed using**

```r
BiocManager::install("r3Cseq")  OR

source("https://biooconductor.org/biocLite.R")

biocLite("r3Cseq")
```

**r3Cseq requires the BSgenome of H. sapiens (hg19). This can be downloaded in R through:**

```r
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")  OR

source("https://bioconductor.org/biocLite.R")

biocLite("BSgenome.Hsapiens.UCSC.hg19"
```

**load the r3Cseq package and the hg19 to R**

```r
>library(r3Cseq)
>library(BSgenome.Hsapiens.UCSC.hg19)
```


**Loading the data set (.bam files)**

###loading the experimental file###

```r
>expFile <- exp_file_name.bam
```

###loading the control file###

```r
contrFile <- contr_file_name.bam
```

**Creating the r3Cseq object**

```r
>object_name <- new("r3Cseq", organismName='hg19',
	alignedReadsBamExpFile=expFile, alignedReadsBamContrFile = contrFile,
	isControllInvolved=TRUE, viewpoint_chromosome='chr#',
	viewpoint_primer_forward='SEQUENCE_OF_FORWARD_PCR_PRIMER',
	viewpoint_primer_reverse='SEQUENCE_OF_REVERSE_PCR_PRIMER',
	expLabel="Name_for_Experimental",
	contrLabel="Name_for_Control",
	restrictionEnzyme='RE')
```


Note that r3Cseq can function without a control. r3Cseq can do the control and experimental file in tandem. It will find which interactions are significant within the samples (internally) through the RPM, p-value, q-value, etc.)

**However, if two samples are being used, it will not perform which interactions that are unique to the control and/or experimental, which interactions are shared. It also does not perform differential analysis of interactions i.e. if the interactions are shared but an interaction is expressed more in one sample than the other, it cannot do statistical analysis if the differences between these interactions are significant.**

For now, one sample can be loaded (no control). This can be done by deleting the control aspect of code. And assigning a FALSE for the argument is Control involved. Compare the code below to the ones above.

```r
>object_name <- new("r3Cseq", organismName='hg19',
	alignedReadsBamExpFile=expFile,
	isControllInvolved=FALSE, viewpoint_chromosome='chr#',
	viewpoint_primer_forward='SEQUENCE_OF_FORWARD_PCR_PRIMER',
	viewpoint_primer_reverse='SEQUENCE_OF_REVERSE_PCR_PRIMER',
	expLabel="Name_for_Experimental",
	restrictionEnzyme='RE')
```


**Analysis of the r3Cseq object**

This line loads the reads from the BAM file that was assigned earlier. Loading may take a while depending on speed of the CPU and RAMM capacity.

```r
>getRawReads(object_name)
```

r3Cseq has two options in terms of counting the reads. The first option is through restriction fragment. The line below performs the counting based on the regions between the restriction fragment

```r
>getReadCountPerRestrictionFragment(object_name)
```

What r3Cseq does is it fragments the hg19 based on the locations of the restriction fragment sites which was provided earlier when creating the r3Cseq object. After it fragments the genome, it counts how many reads (from the aligned bam file) are there in "per restriction fragment" fragment.

The second option is through a user defined window. This is where the user defines how many kilobases are there going to be per fragment. The code below can be used:

```r
>getReadCountPerWindow(object_name, windowSize=20e3)
```

This, in turn will fragment the hg19 every 20kb and proceeds to count how many reads are there in each window.


After counting the raw reads, the data is normalized by finding the RPM (reads per million) per the restriction (or window) fragment using the line below:

```r
>calculateRPM(object_name)
```

After data normalization, statistical analysis and finding interactions are done using the line below. Details about significant interaction calling for the samples are detailed in the Thongjuea 2012 paper which is referrenced at the refences section.

```r
>getInteractions
```

**r3Cseq results, data exportation, etc**

After getting the interactions, the data can now be viewed and exported.

To simply view the interactions, the following can lines can be used (omit control portion if control was not involved)

```r
>expResults <- expInteractionRegions(object_name)
>expReults

>contrResults <- contrInteractionRegions(object_name)
>contrResults
```

Raw reads and normalized reads can be extracted using:

```r
>export3CseqRawReads2bedGraph(object_name)

>export3Cseq2bedGraph(object_name)
```

To see all interactions along the genome (all chromosomes), use the following line:

```r
>plotOverviewInteractions(object_name)
```

To see the interactions near the bait/viewpoint:

```r
>plotInteractionsNearViewPoint(object_name)
```

To see Interactions per chromosome:

```r
>plotInteractionsPerChromosome(object_name,"chr#")
```

Finally, to export the interactions into a text file, use:

```r
>exportInteractions2txt(object_name)
```

This function automatically exports ALL interactions (significant or not) into a text file. Additionally, if there is a control file, this function prints the experimental and control interactions into two separate text files. The name is determined based on the name that was assigned for control and experimental labels when making the r3Cseq object.

Read the interactions .txt file generated

I am only interested in the interactions in the same chr of your bait, in this case it's **chr15**

So I will keep the interactions in that chr, filter out the non-significant ones (with a q-value > 0.1)

And then I will also going to filter interactions that are further away than 100kb from the promoter of RCCD1 (the bait)

And the interactions that are closer than 3kb from the promoter as well. (this is because of the technology of 4C, you can't tell if these interactions that are very close are artifacts or real ones, most likely self-ligation products which we are not interested)

In the end after you do this you'll have a list of interactions in the format of `chr` `start` `end`

Meaning the start and the end of where the promoter of your bait gene is interacting.

So in the end you should have a matrix with chr start end again, but now your start will be all the same because it's the promoter of the gene, and the end will be the middle of each interaction.

You can use the start of that gene as the bait.

Meaning the start and the end of where the promoter of your bait gene is interacting.

So in the end you should have a matrix with chr start end again, but now your start will be all the same because it's the promoter of the gene, and the end will be the middle of each interaction.

You can use the start of that gene as the bait.

```r
## Loading Libraries
library(r3Cseq)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(GenomicRanges)
library(Homo.sapiens)
library(dplyr)
library(rtracklayer)

## MCF interactions
mcf.int <- read.table("MCF7_DpnII.interaction.txt", header = T)
mcf_data <- mcf.int[mcf.int$chromosome == "chr15" & !mcf.int$q.value > 0.9, ] #Filtering the interactions within chr15 and q-value < 0.1

## Genomic Location of RCCD1 with entrez gene id
id <- "91433"
rccdGene <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene,
                  filter=list(gene_id=id))

## start and end of the gene (RCCD1)
rccdPromoter <- start(rccdGene)
mcf_end <- ((mcf_data$start+mcf_data$end)/2)

## Genomic Ranges interaction without the bait
mcf_test <- makeGRangesFromDataFrame(mcf_data)
rccd_test <- GRanges(seqnames = "chr15", ranges = IRanges(rccdPromoter, rccdPromoter+1))
mcf_int <- GenomicRanges::distance(rccd_test, mcf_test)
mcf_int2 <- mcf_test[GenomicRanges::distance(rccd_test, mcf_test) <= 70000 & GenomicRanges::distance(rccd_test, mcf_test) >= 3000] # 55 significant interactions

## Genomic Ranges interaction with bait
mcf <- cbind(rccdPromoter, mcf_end)
mcf_2 <- cbind(mcf, mcf_data$chromosome)
colnames(mcf_2) <- c("end", "start", "chromosome")
mcf_2 <- mcf_2[, c(3, 2, 1)]
mcf_3 <- as.data.frame(mcf_2)
mcf_3$start <- as.integer(mcf_3$start)
mcf_3$end <- as.integer(mcf_3$end)
sapply(mcf_3, class)

mcf_4 <- mcf_3
mcf_4$end <- pmax(mcf_3$start, mcf_3$end)
mcf_4$start <- pmin(mcf_3$start, mcf_3$end) # 97 interactions

mcf_4$diff <- (mcf_4$end - mcf_4$start)
mcf_5 <- mcf_4[mcf_4$diff >= 3000 & mcf_4$diff <= 70000,] #Removing interactions longer than 100kb and shorter than 3kb from the bait
mcf_5$diff <- NULL
mcf_loops <- makeGRangesFromDataFrame(mcf_5)  # 55 interactions

#####################################################################################
#####################################################################################

## BT549 interactions
bt.int <- read.table("BT_DpnII.interaction.txt", header = T)
bt_data <- bt.int[bt.int$chromosome == "chr15" & !bt.int$q.value > 0.9, ] #Filtering the interactions within chr15 and q-value < 0.1

## Genomic Location of RCCD1 with entrez gene id
id <- "91433"
rccdGene <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene,
                  filter=list(gene_id=id))

## start and end of the gene (RCCD1)
rccdPromoter <- start(rccdGene)
bt_end <- ((bt_data$start+bt_data$end)/2)

## Genomic Ranges interaction without the bait
bt_test <- makeGRangesFromDataFrame(bt_data)
rccd_test <- GRanges(seqnames = "chr15", ranges = IRanges(rccdPromoter, rccdPromoter+1))
bt_int <- GenomicRanges::distance(rccd_test, bt_test)
bt_int2 <- bt_test[GenomicRanges::distance(rccd_test, bt_test) <= 70000 & GenomicRanges::distance(rccd_test, bt_test) >= 3000] 

## Genomic Ranges interaction with bait
bt <- cbind(rccdPromoter, bt_end)
bt_2 <- cbind(bt, bt_data$chromosome)
colnames(bt_2) <- c("end", "start", "chromosome")
bt_2 <- bt_2[, c(3, 2, 1)]
bt_3 <- as.data.frame(bt_2)
bt_3$start <- as.integer(bt_3$start)
bt_3$end <- as.integer(bt_3$end)
sapply(bt_3, class)

bt_4 <- bt_3
bt_4$end <- pmax(bt_3$start, bt_3$end)
bt_4$start <- pmin(bt_3$start, bt_3$end)  # 96 interactions

bt_4$diff <- (bt_4$end - bt_4$start)
bt_5 <- bt_4[bt_4$diff >= 3000 & bt_4$diff <= 70000,] #Removing interactions longer than 100kb and shorter than 3kb from the bait
bt_5$diff <- NULL
bt_loops <- makeGRangesFromDataFrame(bt_5) # 53 interactions

## counting overlap interactions in Breast Cancer cell lines
breast_over <- GenomicRanges::intersect(mcf_int2, bt_int2) 
breast_over2 <- GenomicRanges::intersect(mcf_loops, bt_loops)

## Overlap interactions between KURA, UWB, MCF & BT549
overall_over <- GenomicRanges::intersect(breast_over, ovarian_over)
## Converting GRanges to bed file
rtracklayer::export.bed(breast_over, "Breast cell lines overlap interactions.bed")
rtracklayer::export.bed(overall_over, "All 4 cell lines overlap interactions.bed")

rtracklayer::export.bed(mcf_loops, "MCF_arc_interactions.bed")
rtracklayer::export.bed(bt_loops, "BT_arc_interactions.bed")
```









