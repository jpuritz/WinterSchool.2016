#B@G 2016 Mapping Exercise
Designed by Jon Puritz

#GOALS
1.	To learn how to map reads using the MEM algorithm of BWA
2.	To learn about the SAM/BAM file format 
3.	To learn how to use SAMtools to manipulate and view SAM/BAM files
4.	To learn how optimize read mappings by examining alignment statistics
5.	BONUS GOAL- learn how to use sam-stats from ea-utils.

#Exercise

Let's find our way back to your original working directory and make sure the dDocent module is loaded
```bash	
cd ~/D1W
module load dDocent/v2.18 
```
We are going to use the program BWA to learn about read mapping.  Specifically, we are going to be working with the
MEM algorithm of the program.  BWA-MEM is a fast, customizable, and accurate read mapping software.  Look at comparisons in
(http://arxiv.org/pdf/1303.3997.pdf) and (http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0090581).  
Let's take a look at bwa mem using some of our simulated data.
First, let's make sure we all use the same reference.
```bash	
bash remake_reference.sh 4 4 0.90 PE 8
```
This should return a reference of 1000 sequences.
First thing we need to do is to index the reference.fasta file, so that bwa search against it more quickly.
```bash	
bwa index reference.fasta
```
Now, let's align one of the samples to the reference, as simply as possible.  This means calling bwa mem and feeding it the reference and reads files
```bash	
bwa mem reference.fasta PopA_01.F.fq.gz PopA_01.R.fq.gz > PopA_01.sam
```
We saved the output in the file PopA_01.sam.  A Sequence Alignment Map or SAM file is the expanded the standard output for most mapping software.
A lot of information is stored in a SAM file.  We won't be able to cover it all today, so please be sure to read through the official description
(http://samtools.github.io/hts-specs/SAMv1.pdf).

Let's take a look at the SAM file we generated.
```bash	
head PopA_01.sam
```
```
	@SQ	SN:E1_L100	LN:206
	@SQ	SN:E2_L100	LN:206
	@SQ	SN:E3_L100	LN:206
	@SQ	SN:E4_L100	LN:206
	@SQ	SN:E5_L100	LN:206
	@SQ	SN:E6_L100	LN:206
	@SQ	SN:E7_L100	LN:206
	@SQ	SN:E8_L100	LN:206
	@SQ	SN:E9_L100	LN:206
	@SQ	SN:E10_L101	LN:207
```
This is the header section of the SAM file.  Each line begins with an @ and lists the referenece sequence name (SN:) and its length (LN:206).
The bulk of the data in a SAM file is in the alignment section though.  To see that, let's look at the first two lines that do not start with the @ character.
```bash	
mawk '!/@/' PopA_01.sam	| head -2
```
```
lane1_fakedata0_0	99	dDocent_Contig_702	2	60	94M	=	106	204	AATTCGGCTTGCAACGCAAGTGACGATTCCCACGGACATAACTGATCTAAGTAACTTCCAAATCTGGGAATGGGATTTCATAATTAAGGACTAT	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:0	MD:Z:94	AS:i:94	XS:i:0
lane1_fakedata0_0	147	dDocent_Contig_702	106	60	100M	=	2	-204	ACGACGAGCAATCCACAGACCTAGGCCCATCGAAGCGTCTTATGATTGATAACATCAGAGGGGGATGGGAGGTCCTGCTGTCGCATGGGAGAATACACCG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	NM:i:2	MD:Z:57A41C0	AS:i:94	XS:i:0
```
This output is a tab delimited text file, and each column holds a specific type of data.

* Col 1	=	the name of the read
* Col 2	=	This is a bitwise flag that describes the alignment.  The flag 99 means that this is the first read in a read pair, that this read is paired, is paired in a proper pair, and that it's mate is mapped on the one the reverse strand.
   These bitwise flags can be difficult to interpret with out help.  Checkout this website for a translator (http://broadinstitute.github.io/picard/explain-flags.html).  
Go to the website and translate what flag 147 means. 
* Col 3	=	The name of the reference sequence/contig that the read has mapped to.
* Col 4	=	The position of the left most base of the read on the reference sequence.  All the reference sequences in our example have a leading N, so 2 means it aligned to the beginning of the reference.
* Col 5	=	This is the map quality score.  It equals −10 log10 Pr{mapping position is wrong}, rounded to the nearest integer.  This is identical to a PHRED score, so a score of 60 means that there is a 0.000001 chance that the mapping is wrong.
* Col 6	=	The CIGAR string.  That's the Compact Idiosyncratic Gapped Alignment Report.  It's a single string of text that describe how the sequence aligns, including matches (M), insertions (I), deletions (D), mismatch (X), soft clipping (S), and hard clipping (H), etc.
   94M translates to 94 matching basepairs starting from the left most position (reported in Col 4).
* Col 7	=	The name of the reference sequence the read's mate maps to.  The equal sign here means the both map to the same reference contig.
* Col 8	=	The alignment position of the read's mate
* Col 9	=	The observed length of the total alignment.  If both reads map to the same contig it equals the number of bases from the leftmost mapped base to the rightmost mapped base.
   The leftmost segment has a plus sign and the rightmost has a minus sign.
* Col 10	=	The sequence of the read
* Col 11	=	The ASCII quality score for each basepair in the read.

This info is important to know and can be helpful for quick filtering.  However, going through an entire alignment in assessment would be tedious. 
In addition, all this data uses up a lot of disk space.  To remedy this, SAM files are often (and should be) converted into Binary Alignment Mao (BAM) Files.
The program SAMtools can be used to convert between the two formats.  SAMtools can also view, manipulate, and output statistics about SAM/BAM files.  Below, we are going to use it to examine 
the effect of different mapping parameters.  Let's use SAMtools to evaluate our mapping.  
```bash	
samtools view -Sbt reference.fasta PopA_01.sam | samtools flagstat -
```
```
[bam_header_read] EOF marker is absent. The input is probably truncated.
[samopen] SAM header is present: 1000 sequences.
37540 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
37370 + 0 mapped (99.55%:-nan%)
37540 + 0 paired in sequencing
18770 + 0 read1
18770 + 0 read2
25894 + 0 properly paired (68.98%:-nan%)
37200 + 0 with itself and mate mapped
170 + 0 singletons (0.45%:-nan%)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
The above code is converting our SAM file to a BAM file and piping it directly back to SAMtools to generate the stats you see.  The first two lines are errors from the piping, but don't affect the data.
Most of the above stats are self-explanatory, but let's key on a few of them.
The first line shows the total number of reads and line 3 shows how many of them mapped (and what percentage).  
You can also see that 0.71% of the mappings were singletons (only one read from a pair) and no reads had mates mapped to different reference contigs.
The last important figure is the number of reads (and percentage) that are properly paired.  Proper pairings are determined by the aligning software.
In general, the reads must be mated to the same reference contig, in the proper orientation, with high mapping quality, and have fragment size that meets expectations (usually based on the mean and SD of all mappings).
You can see that just with the defaults, the mapping performed very well.  Most of the reads mapped with no discordant mate mappings.  The one seemingly odd number is the low percentage of proper pairings.
Why might this be the case?

To remedy the problem, we can add a few extra parameters to our bwa code.  Basically, setting the limits for the proper insert size distribution.
```bash	
bwa mem reference.fasta PopA_01.F.fq.gz PopA_01.R.fq.gz -I 200,40 2>/dev/null | samtools view -SbT reference.fasta - | samtools flagstat -
```
```
37540 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
37540 + 0 mapped (100.00%:-nan%)
37540 + 0 paired in sequencing
18770 + 0 read1
18770 + 0 read2
37540 + 0 properly paired (100.00%:-nan%)
37540 + 0 with itself and mate mapped
0 + 0 singletons (0.00%:-nan%)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
The -I parameter sets the insert size at 200 and the standard deviation at 40.  Now the proper pairings jump up to 100% or all of the mapped reads.

Let's look at a second sample:
```bash
bwa mem reference.fasta PopA_02.F.fq.gz PopA_02.R.fq.gz -I 200,40 2>/dev/null | samtools view -SbT reference.fasta - | samtools flagstat -
```
```
38296 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
38295 + 0 mapped (100.00%:-nan%)
38296 + 0 paired in sequencing
19148 + 0 read1
19148 + 0 read2
38294 + 0 properly paired (99.99%:-nan%)
38294 + 0 with itself and mate mapped
1 + 0 singletons (0.00%:-nan%)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```	

There are still a very small percentage of reads that haven't mapped.  One way to optimize this is to change the scoring of alignments.
In short, alignments are scored by multiplying the number of matching bases by a match score (default 1), then the number of mismatching bases are multiplied by a mismatch penalty (default 4) and subtracted from the match score.  Next any gaps are penalized by a single gap opening penalty (default 6) and a gap extension penalty (default 1) multiplied the length of the gap.  The easiest way to make the alignment a little more liberal is to just decrease the mismatch and gap opening penalties by 1.
```bash	
bwa mem reference.fasta PopA_02.F.fq.gz PopA_02.R.fq.gz -I 200,40 -B 3 -O 5 2>/dev/null | samtools view -SbT reference.fasta - | samtools flagstat -
```
```
38296 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
38295 + 0 mapped (100.00%:-nan%)
38296 + 0 paired in sequencing
19148 + 0 read1
19148 + 0 read2
38294 + 0 properly paired (99.99%:-nan%)
38294 + 0 with itself and mate mapped
1 + 0 singletons (0.00%:-nan%)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

```	
Well, that didn't change anything.  These might be reads that don't match to any reference contig.  We can use mawk to quickly examine some of the reads.
```bash	
bwa mem reference.fasta PopA_02.F.fq.gz PopA_02.R.fq.gz -I 200,40 -B 3 -O 5 2>/dev/null |  mawk '/\*/' | head
```
```
lane1_fakedata837_6	133	dDocent_Contig_271	2	0	*	=	2	0	CGGGCAAATAGGCATGTGAACGTATTACCTCTCAGGCGCTTCTCTCGCGGTCGTTCAACCACTCAGTGATAAAAACGGTAAACAGGGCCTGTTAAGATTA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	AS:i:0	XS:i:0
```
Here, you can see that the reverse reads in this pair does not mapped to any contig.  This could mean that the second read is particularly divergent or that the reference is not completely representative.  However, at 99.99% we are doing well. 

To get a better idea of how to optimize mapping further, let's take a look at some real (messy) data.
```bash
mkdir realdata
cd realdata
ln -s /gdc_home5/groups/bag2016/monday/mapping2/* .
bwa mem reference.fasta JC_1119.R1.fq.gz JC_1119.R2.fq.gz -I 200,40 -t 8 2>/dev/null | samtools view -@8 -q 1 -SbT reference.fasta - | samtools flagstat -
```
```
2327545 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
2919 + 0 supplementary
0 + 0 duplicates
2327545 + 0 mapped (100.00% : N/A)
2324626 + 0 paired in sequencing
1172770 + 0 read1
1151856 + 0 read2
2157832 + 0 properly paired (92.82% : N/A)
2266708 + 0 with itself and mate mapped
57918 + 0 singletons (2.49% : N/A)
107084 + 0 with mate mapped to a different chr
81517 + 0 with mate mapped to a different chr (mapQ>=5)
```
If you're paying attention, you would have noticed that I added another parameter to the bwa mem command (-t).  This sets the number of processors to use.  This file is much larger.  I also added -@ and -q parameters to samtools
The -q 1 removes reads with 0 probability of mapping from being retained in the file and -@16 adds multithreading.  

Now, we have a baseline to work with for this individual.  Let's try to optimize it.  First, let's relax the mismatch and gap opening penalties.
```bash
bwa mem reference.fasta JC_1119.R1.fq.gz JC_1119.R2.fq.gz -I 200,40 -t 8 -B 3 -O 5 2>/dev/null | samtools view -@8 -q 1 -SbT reference.fasta - | samtools flagstat -
```
```
1803670 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
1803670 + 0 mapped (100.00%:-nan%)
1803670 + 0 paired in sequencing
911530 + 0 read1
892140 + 0 read2
1758193 + 0 properly paired (97.48%:-nan%)
1777652 + 0 with itself and mate mapped
26018 + 0 singletons (1.44%:-nan%)
19131 + 0 with mate mapped to a different chr
17882 + 0 with mate mapped to a different chr (mapQ>=5)
```
So, we've improved things somewhat.  We've mapped more reads and mapped more less with discordant mates, but decreased the proper pairings %, increased our singletons.
This discordant mappings, however, make up a very small percentage of the total mappings, around 1.5%.  Previously, it was 1.46%.  This probably is a good stopping point for
trying to increase the number of mappings.  The more we add at this point, the worse they are going to be.  Let's take some steps to reduce bad mappings.

We can add an additional flag to bwa that specifically aids with RAD mappings.  RADseqs tend to be more conserved at the 5' end simply because of the cut site.  In general, if the
beginning of a read needs to be trimmed off to match (especially after going through quality filtering), it most likely is not homologous to that locus and should be removed.
```bash
bwa mem reference.fasta JC_1119.R1.fq.gz JC_1119.R2.fq.gz -I 200,40 -t 8 -L 20,5 -B 3 -O 5 2>/dev/null | samtools view -@8 -q 1 -SbT reference.fasta - | samtools flagstat -
```
```
2329008 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
2314 + 0 supplementary
0 + 0 duplicates
2329008 + 0 mapped (100.00% : N/A)
2326694 + 0 paired in sequencing
1174336 + 0 read1
1152358 + 0 read2
2158731 + 0 properly paired (92.78% : N/A)
2268022 + 0 with itself and mate mapped
58672 + 0 singletons (2.52% : N/A)
107486 + 0 with mate mapped to a different chr
79551 + 0 with mate mapped to a different chr (mapQ>=5)
```
The improvement is slight but noticeable.  More reads are properly paired and there are less singletons and discordant mappings!

We can take this a step further by also filtering out reads that have a significant amount of clipping.  
Can you think of a way how to do this?

Remember that while the alignment is in SAM format, column 6 contains the CIGAR string.  If the read has been clipped at all (hard or soft), there will be an S or H
character proceed by a number (the number of bases clipped).  So, we can use awk to filter via a regular expression to remove reads with significant clipping (for example more than 10 bases).
```bash
bwa mem reference.fasta JC_1119.R1.fq.gz JC_1119.R2.fq.gz -I 200,40 -t 8 -L 20,5 -B 3 -O 5 2>/dev/null | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/'|samtools view -@8 -q 1 -SbT reference.fasta - | samtools flagstat -
``` 
```
2272415 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
2272415 + 0 mapped (100.00% : N/A)
2272415 + 0 paired in sequencing
1146942 + 0 read1
1125473 + 0 read2
2124931 + 0 properly paired (93.51% : N/A)
2221470 + 0 with itself and mate mapped
50945 + 0 singletons (2.24% : N/A)
94842 + 0 with mate mapped to a different chr
68401 + 0 with mate mapped to a different chr (mapQ>=5)
```
While this step definitely decreased the overall mappings by a small percentage it increased our properly paired percentage and decreased the singleton and discordant mappings.  

Let's try un relaxing the mapping parameters and compare the results with the previous
```bash
bwa mem reference.fasta JC_1119.R1.fq.gz JC_1119.R2.fq.gz -I 200,40 -t 3 -B 3 -O 5 -L 20,5 2>/dev/null | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/'| samtools view -@16 -q 1 -SbT reference.fasta - | samtools flagstat - > relaxed.stats
bwa mem reference.fasta JC_1119.R1.fq.gz JC_1119.R2.fq.gz -I 200,40 -t 3 -L 20,5 2>/dev/null | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/'| samtools view -@16 -q 1 -SbT reference.fasta - | samtools flagstat - > normal.stats
paste <(cut -f1 -d + relaxed.stats) <(cut -f1 -d + normal.stats) <(cut -f4-10 -d " " relaxed.stats) 
```
```
2272415 	2263752 	in total (QC-passed reads + QC-failed reads)
0 			0 			secondary
0 			0 			supplementary
0 			0 			duplicates
2272415 	2263752 	mapped (100.00% : N/A)
2272415 	2263752 	paired in sequencing
1146942 	1140984 	read1
1125473 	1122768 	read2
2124931 	2120013 	properly paired (93.51% : N/A)
2221470 	2215558 	with itself and mate mapped
50945 		48194 		singletons (2.24% : N/A)
94842 		93866 		with mate mapped to a different chr
68401 		69122 		with mate mapped to a different chr
```	
I cheated and added extra tabs to make the table better.  We can see that the relaxed stats lead to more mappings and  more proper pairings.
The only drawback to the more relaxed mappings is that there are more discordant mapping pairs and singletons.  These extra discordant mappings may be of very low quality and

Now on your own, try repeating the previous steps with a mapping quality cutoff of 5 and 10.  Does the same pattern hold true?
Your results should look like this:
```
2237573 	2227393 	in total (QC-passed reads + QC-failed reads)
0 			0		 	secondary
0 			0 			supplementary
0 			0		 	duplicates
2237573 	2227393 	mapped (100.00% : N/A)
2237573 	2227393 	paired in sequencing
1137155 	1129926 	read1
1100418 	1097467 	read2
2123913 	2114718 	properly paired (94.92% : N/A)
2192739 	2184255 	with itself and mate mapped
44834 		43138 		singletons (2.00% : N/A)
68401 		69122 		with mate mapped to a different chr
68401 		69122 		with mate mapped to a different chr 	

```
```

2208568 	2206988 	in total (QC-passed reads + QC-failed reads)
0		 	0 			secondary
0 			0 			supplementary
0 			0 			duplicates
2208568 	2206988 	mapped (100.00% : N/A)
2208568 	2206988 	paired in sequencing
1116047 	1112705 	read1
1092521 	1094283 	read2
2116636 	2112749 	properly paired (95.84% : N/A)
2171508 	2170128 	with itself and mate mapped
37060 		36860 		singletons (1.68% : N/A)
54457 		56970 		with mate mapped to a different chr
54457 		56970 		with mate mapped to a different chr
```
The benefit of relaxing the mapping parameters is even more pronounced at higher mapping qualities.
This looks like an optimized mapping setting for this sample.
Remember that you should trying using multiple samples from your data set before picking final mapping parameters

Bonus Material
Here are some other things you can try if you have extra time.

Another program from the ea-utils suite called sam-stats can be helpful to evaluate alignments as well.
You can quickly download and install it with the following commands:
```bash
curl -O https://ea-utils.googlecode.com/files/ea-utils.1.1.2-537.tar.gz
tar xvzf ea-utils.1.1.2-537.tar.gz
cd ea-utils.1.1.2-537
make
cp sam-stats ~/bin
cd ..
rm -rf ea-utils.1.1.2-537*
```
You can pipe results to this program, just like samtools, so you don't need to waste time writing files while optimizing.
Here is an example from one of our exercises above:
```bash
bwa mem reference.fasta JC_1119.R1.fq.gz JC_1119.R2.fq.gz -I 200,40 -t 8 -B 3 -O 5 -L 20,5 2>/dev/null | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/'| samtools view -q 1 -@16 -SbT reference.fasta - |  sam-stats -D -B 
```
```
reads	2272415
version	1.34.488
mapped reads	2272409
ambiguous	6
pct ambiguous	0.000251
max dup align	3
singleton mappings	121567
total mappings	2272415
mapped bases	219626385
library	paired-end
discordant mates	94842
pct forward	97.756
phred	33
forward	2221426
reverse	50989
len max	100
len mean	96.6489
len stdev	5.3856
mapq mean	56.8529
mapq stdev	11.1001
mapq Q1	60.00
mapq median	60.00
mapq Q3	60.00
snp rate	0.013529
ins rate	0.002009
del rate	0.001809
pct mismatch	51.7389
insert mean	206.5219
insert stdev	0.6871
insert Q1	206.00
insert median	206.00
insert Q3	207.00
base qual mean	37.1035
base qual stdev	4.2683
%A	28.4294
%C	22.3528
%G	20.5673
%T	28.6492
%N	0.0012
num ref seqs	40171
num ref aligned	26760
```
The useful info that sam-stats reports that SAMtools doesn't is the number (and percentage) of ambiguous mappings (6 in this case), and the number of reference sequences that the reads aligned to.  These can both be helpful for determining optimal mapping parameters.
Ideally, you would like to maximize the number of reference sequences mapped to and minimize the number of ambiguous mappings.  







