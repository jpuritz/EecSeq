# This code repeats the analyses done in Puritz and Lotterhos 2017


This code assumes that there are three subdirectories in this working directory
./RNA ./DNA which contain the raw sequencing files from the SRA BioProject PRJNA423022 and 
./Genome directory which contains files from the eastern oyster genome 
(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0)

## Setup genome
```bash
cd $WORKING_DIR/Genome
cd ./Assembled_chromosomes/seq
```
Assemble full reference from each assembled chromosome
```bash
zcat 6565_ref_C_virginica-3.0_chrMT.fa.gz 6565_ref_C_virginica-3.0_chr[0-9]*.fa.gz | sed 's/ref|//g' | sed 's/| Crass/  Crass/g' > genome.fasta
```
Create a "genome.file" to use for bedtools sorting
```bash
samtools faidx genome.fasta
mawk -v OFS='\t' {'print $1,$2'} genome.fasta.fai > genome.file
```
Switch to annotations section and create the needed bed files
```bash
cd $WORKING_DIR/Genome/GFF
ln -s $WORKING_DIR/Genome/Assembled_chromosomes/seq/genome.file .
```
This next step requires BEDOPS (https://github.com/bedops) to be installed on your system
```bash
gff2bed < ref_C_virginica-3.0_top_level.gff3 > ref_C_virginica-3.0_top.bed
```
Create file bed file with all exons
```bash
mawk '$8 ~ /exon/' ref_C_virginica-3.0_top.bed > ref3.0.exon.bed
```
Sort the bed file
```bash
bedtools sort -i ref3.0.exon.bed -faidx <(cut -f1 genome.file) > sorted.ref3.0.exon.bed
```
Remove duplicates and concatenate overlapping intervals from transcript variants
```bash
bedtools merge -i sorted.ref3.0.exon.bed > sorted.ref3.0.exon.sc.bed
```
Create bed file for all gene regions
```bash
mawk '$8 ~ /gene/' ref_C_virginica-3.0_top.bed > ref3.0.gene.bed
bedtools sort -i ref3.0.gene.bed -faidx <(cut -f1 genome.file) > sorted.ref3.0.gene.bed
bedtools merge -i sorted.ref3.0.gene.bed > sorted.ref3.0.gene.sc.bed
```
Create bed file for intergenic, intron, non-coding, and CDS regions
```bash
bedtools complement -i sorted.ref3.0.gene.bed -g genome.file -sorted | bedtools subtract -b sorted.ref3.0.exon.sc.bed -a - > cv.ref3.intergenic.bed
bedtools complement -i sorted.ref3.0.exon.sc.bed -g genome.file -sorted > cv.ref3.noncoding.bed
bedtools intersect -a cv.ref3.noncoding.bed -b sorted.ref3.0.gene.sc.bed -sorted > cv.ref3.intron.bed
mawk '$8 ~ /CDS/' ref_C_virginica-3.0_top.bed > ref3.0.CDS.bed
bedtools sort -i ref3.0.CDS.bed -faidx <(cut -f1 genome.file) > sorted.ref3.0.CDS.bed
bedtools merge -i sorted.ref3.0.CDS.bed > sorted.ref3.0.CDS.sc.b
bedtools sort -i sorted.ref3.0.CDS.sc.b -faidx <(cut -f1 genome.file) > sorted.ref3.0.CDS.sc.bed
```
Create bed file for untranslated regions (UTRs) of exons

```bash
bedtools subtract -a sorted.ref3.0.exon.bed -b sorted.ref3.0.CDS.bed -g genome.file -sorted > sorted.ref3.0.UTR.bed
bedtools merge -i <(bedtools sort -i sorted.ref3.0.UTR.bed -faidx <(cut -f1 genome.file))> sorted.ref3.0.UTR.sc.b
bedtools sort -i sorted.ref3.0.UTR.sc.b -faidx <(cut -f1 genome.file) > sorted.ref3.0.UTR.sc.bed
```
Create bed file for mtDNA
```bash
mawk '$1 ~ /NC_007175.2/' ref_C_virginica-3.0_top.bed > mtDNA.bed
```
## RNA trimming and custom adapter removal

`cd $WORKING_DIR/RNA`

Download custom adapter file

`wget https://raw.githubusercontent.com/jpuritz/EecSeq/master/Bioinformatics/TruSeq2-PE.fa`

Download custom dDocent v2.2.ed20 that looks for adapters in working directory

`wget https://raw.githubusercontent.com/jpuritz/EecSeq/master/Bioinformatics/dDocent`

Download the configuration file

`wget https://raw.githubusercontent.com/jpuritz/EecSeq/master/Bioinformatics/RNA.config`

Use nano to alter configuration file for your system

`nano RNA.config`

Run dDocent to trim files

`./dDocent RNA.config`

## RNA mapping to the genome

#### The next steps require STAR (https://github.com/alexdobin/STAR) to be installed 
```bash
ln -s $WORKING_DIR/Genome/Assembled_chromosomes/seq/virginica.ref3.0.fasta reference.fasta
ln -s $WORKING_DIR/Genome/GFF/ref_C_virginica-3.0_top_level.gff3 top.gff
mkdir genome
```
Create genome files with STAR
```bash
STAR --runMode genomeGenerate --runThreadN 64 --genomeDir ./genome --genomeFastaFiles reference.fasta --sjdbGTFfile top.gff --sjdbGTFtagExonParentTranscript Parent--sjdbOverhang 149
```
Run Star in dual-pass mode
```bash
STAR --runThreadN 32 --genomeDir ./genome/ --readFilesIn  RNASeq_1.R1.fq.gz RNASeq_1.R2.fq.gz  --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 --outSAMmapqUnique 40 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix RNA1m2  --readFilesCommand zcat
STAR --runThreadN 32 --genomeDir ./genome/ --readFilesIn  RNASeq_2.R1.fq.gz RNASeq_2.R2.fq.gz  --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 --outSAMmapqUnique 40 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix RNA2m2  --readFilesCommand zcat
STAR --runThreadN 32 --genomeDir ./genome/ --readFilesIn  RNASeq_3.R1.fq.gz RNASeq_3.R2.fq.gz  --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 --outSAMmapqUnique 40 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix RNA3m2  --readFilesCommand zcat
STAR --runThreadN 32 --genomeDir ./genome/ --readFilesIn  RNASeq_4.R1.fq.gz RNASeq_4.R2.fq.gz  --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 --outSAMmapqUnique 40 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix RNA4m2  --readFilesCommand zcat
```
```bash
STAR --runThreadN 32 --genomeDir ./genome/ --readFilesIn  RNASeq_1.R1.fq.gz RNASeq_1.R2.fq.gz  --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 --outSAMmapqUnique 40 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix RNA1m3  --readFilesCommand zcat --sjdbFileChrStartEnd ./RNA1m2SJ.out.tab ./RNA2m2SJ.out.tab ./RNA3m2SJ.out.tab ./RNA4m2SJ.out.tab
STAR --runThreadN 32 --genomeDir ./genome/ --readFilesIn  RNASeq_2.R1.fq.gz RNASeq_2.R2.fq.gz  --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 --outSAMmapqUnique 40 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix RNA2m3  --readFilesCommand zcat --sjdbFileChrStartEnd ./RNA1m2SJ.out.tab ./RNA2m2SJ.out.tab ./RNA3m2SJ.out.tab ./RNA4m2SJ.out.tab
STAR --runThreadN 32 --genomeDir ./genome/ --readFilesIn  RNASeq_3.R1.fq.gz RNASeq_3.R2.fq.gz  --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 --outSAMmapqUnique 40 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix RNA3m3  --readFilesCommand zcat --sjdbFileChrStartEnd ./RNA1m2SJ.out.tab ./RNA2m2SJ.out.tab ./RNA3m2SJ.out.tab ./RNA4m2SJ.out.tab
STAR --runThreadN 32 --genomeDir ./genome/ --readFilesIn  RNASeq_4.R1.fq.gz RNASeq_4.R2.fq.gz  --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 --outSAMmapqUnique 40 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix RNA4m3  --readFilesCommand zcat --sjdbFileChrStartEnd ./RNA1m2SJ.out.tab ./RNA2m2SJ.out.tab ./RNA3m2SJ.out.tab ./RNA4m2SJ.out.tab
```
Use samtools to merge each library into a single alignment file
```bash
samtools merge -@64 m4.merged.bam RNA1m3Aligned.sortedByCoord.out.bam RNA2m3Aligned.sortedByCoord.out.bam RNA3m3Aligned.sortedByCoord.out.bam RNA4m3Aligned.sortedByCoord.out.bam
```
Use samtools and mawk to filter out reads that did not uniquely map and reads that had significant hard or soft clipping >79 bp
```bash
samtools view -@64 -q4 -h -F 0x100 -F 0x400 m4.merged.bam| mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | samtools view -b > m4.q4.merged.bam
```

## Probe/RNA Coverage Analysis

Create symlinks to all bed files and the bedtools genome file
```bash
ln -s $WORKING_DIR/Genome/GFF/*.bed .
ln -s $WORKING_DIR/Genome/GFF/genome.file .
```
Use samtools to calculate the number of reads mapping, reads mapping to genes, and reads mapping to exons
```bash
paste <(samtools view -@32 -c m4.q4.merged.bam) <(samtools view -@32 -c -L sorted.ref3.0.gene.bed m4.q4.merged.bam) <(samtools view -@32 -c -L sorted.ref3.0.exon.sc.bed m4.q4.merged.bam)
```
Output:

`21990025	17234677	12059266`

Use bedtools to calculate per base pair coverage levels across various genomic regions
```bash
bedtools coverage -hist -b m4.q4.merged.bam -a cv.ref3.intron.bed -g genome.file -sorted  -split | grep ^all > AllRNAm4q4.hist.AllIntron.all.split.txt
bedtools coverage -hist -b m4.q4.merged.bam -a cv.ref3.intergenic.bed -g genome.file -sorted -split | grep ^all > AllRNAm4q4.hist.AllIntergenic.all.split.txt
bedtools coverage -hist -b m4.q4.merged.bam -a sorted.ref3.0.exon.sc.bed -g genome.file -sorted  -split | grep ^all > AllRNAm4q4.hist.AllExon.all.split.txt

bedtools coverage -hist -b m4.q4.merged.bam -a sorted.ref3.0.UTR.sc.bed -g genome.file -sorted  -split | grep ^all > AllRNAm4q4.hist.allUTR.all.split.txt
bedtools coverage -hist -b m4.q4.merged.bam -a sorted.ref3.0.CDS.sc.bed -g genome.file -sorted  -split | grep ^all > AllRNAm4q4.hist.AllCDS.all.split.txt
```

### R Code to recreate Figure 2 from above data

```R
setwd("$WORKING_DIR/RNA")
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)

print(files <- list.files(pattern="AllRNAm4q4.hist."))

labs <- c("CDS","Intergenic","Intron","UTR")


cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")

cov <- list()
for (i in 1:length(files)) {
  cov[[i]] <- read.table(files[i])[,c(2,5)]
  cov_cumul=1-cumsum(cov[[i]][,2])
  cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
  cov[[i]]$sample=labs[i]
}

cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")

pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7","#0072B2")

p <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(0,100)+
  scale_alpha(guide = 'none') +
  geom_line(size=1.5)+ 
  #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
  scale_color_manual(values=pcbPalette) +
  scale_fill_manual(values=pcbPalette) +
  ylab("% of Bases > Depth")+
  xlab("Depth")+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme(legend.title = element_blank()) + 
  theme(legend.position=c(0.75,0.75))

png(filename="Figure2.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
p
dev.off()
```
## DNA Adapter Trimming

Change to DNA directory
```bash
cd $WORKING_DIR/DNA
```
Symlink to Genome
```bash
ln -s $WORKING_DIR/Genome/Assembled_chromosomes/seq/virginica.ref3.0.fasta reference.fasta
```
Download same dDocent version, configuration files, and custom adapter filtering script
```bash
wget https://raw.githubusercontent.com/jpuritz/EecSeq/master/Bioinformatics/dDocent_ngs.sh
wget https://raw.githubusercontent.com/jpuritz/EecSeq/master/Bioinformatics/configDNA1
wget https://raw.githubusercontent.com/jpuritz/EecSeq/master/Bioinformatics/configDNA2
wget https://raw.githubusercontent.com/jpuritz/EecSeq/master/Bioinformatics/Adapter_filter.sh
```
Make both dDocent and filter script executable

```bash
chmod +x dDocent_ngs.sh
chmod +x Adapter_filter.sh
```
Run dDocent with first configuration file to run standard quality trimming and adapter removal

```bash
./dDocent_ngs.sh configDNA1
```
Run custom adapter filtering script
```bash
cat namelist | parallel ./Adapter_filter.sh {}
```
## Map genomic reads to genome
Run dDocent again with second configuration file for read mapping
```bash
./dDocent_ngs.sh configDNA2
```
Use Picard [(http://broadinstitute.github.io/picard/)](http://broadinstitute.github.io/picard/) to mark duplicates
```bash
java -Xms4g -jar /shared_lab/scripts/picard.jar MarkDuplicatesWithMateCigar I=ECI_1-RG.bam O=ECI_1-RGmd.bam M=ECI_1_dup_metrics.txt MINIMUM_DISTANCE=300 &> md.ECI1.log 
java -Xms4g -jar /shared_lab/scripts/picard.jar MarkDuplicatesWithMateCigar I=ECI_2-RG.bam O=ECI_2-RGmd.bam M=ECI_2_dup_metrics.txt MINIMUM_DISTANCE=300 &> md.ECI2.log 
java -Xms4g -jar /shared_lab/scripts/picard.jar MarkDuplicatesWithMateCigar I=ECI_3-RG.bam O=ECI_3-RGmd.bam M=ECI_3_dup_metrics.txt MINIMUM_DISTANCE=300 &> md.ECI3.log 
java -Xms4g -jar /shared_lab/scripts/picard.jar MarkDuplicatesWithMateCigar I=ECI_4-RG.bam O=ECI_4-RGmd.bam M=ECI_4_dup_metrics.txt MINIMUM_DISTANCE=300 &> md.ECI4.log 
java -Xms4g -jar /shared_lab/scripts/picard.jar MarkDuplicatesWithMateCigar I=ECI_7-RG.bam O=ECI_7-RGmd.bam M=ECI_7_dup_metrics.txt MINIMUM_DISTANCE=300 &> md.ECI7.log 
java -Xms4g -jar /shared_lab/scripts/picard.jar MarkDuplicatesWithMateCigar I=ECI_12-RG.bam O=ECI_12-RGmd.bam M=ECI_12_dup_metrics.txt MINIMUM_DISTANCE=300 &> md.ECI12.log 
```
Use SAMtools to remove duplicates, secondary alignments, mappings with a quality score less than ten, and reads with more than 80 bp clipped
```bash
samtools view -@32 -h -F 0x100 -q 10 -F 0x400 ECI_12-RGmd.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 64 -b > ECI_12.F.bam 
samtools view -@32 -h -F 0x100 -q 10 -F 0x400 ECI_1-RGmd.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 64 -b > ECI_1.F.bam
samtools view -@32 -h -F 0x100 -q 10 -F 0x400 ECI_4-RGmd.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 64 -b > ECI_4.F.bam 
samtools view -@32 -h -F 0x100 -q 10 -F 0x400 ECI_7-RGmd.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 64 -b > ECI_7.F.bam 
samtools view -@32 -h -F 0x100 -q 10 -F 0x400 ECI_3-RGmd.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 64 -b > ECI_3.F.bam 
samtools view -@32 -h -F 0x100 -q 10 -F 0x400 ECI_2-RGmd.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/'| samtools view -@ 64 -b > ECI_2.F.bam
```

## Generate data in Table 2

```bash
for i in `cat namelist`; do nom=$(samtools view -@32 $i.F.bam -c -L mtDNA.bed); denom=$(samtools view -@32 $i.F.bam -c); dup=$(mawk '/Unknown/' "$i"_dup_metrics.txt | cut -f9); paste <(echo $i) <(echo $(( `zcat $i.F.fq.gz | wc -l` /2 ))) <(echo $(( `zcat $i.R1a.fq.gz | wc -l` /2 ))) <(samtools view -@ 32 $i-RG.bam -c) <(python -c "print(round("$dup" * 100,2))") <(echo $denom) <(python -c "print(round("$nom"/"$denom" *100,2))"); done > data.table2
echo -e "Pool\tRaw_Reads\tFiltered_Reads\tMapped_Reads\t%_Duplicate\tFiltered_Mapped_Reads\t%_mapping_to_mitochondrial_genome" > header
cat header data.table2 > table2.txt
```

## Generate data for Figure 3

Merge all capture pools into one file
```bash
samtools merge -@64 filter.merged.bam ECI_1.F.bam ECI_2.F.bam ECI_3.F.bam ECI_4.F.bam ECI_7.F.bam ECI_12.F.bam
```
Get total coverage counts per exon
```bash
bedtools coverage -b filter.merged.bam -a sorted.ref3.0.exon.sc.bed -sorted -g genome.file -split -counts > cov.counts.filtered.merged.exon.stats
```
Change to RNA directory
```bash
cd $WORKING_DIR/RNA
```
Calculate total RNA coverage per exon
```bash
bedtools coverage -b m4.q4.merged.bam -a sorted.ref3.0.exon.sc.bed -sorted -g genome.file -counts -split > cov.m4q4.EiR.stats
```
Paste RNA and DNA data together and remove mtDNA data
```bash
paste cov.m4q4.EiR.stats <(cut -f4 $WORKING_DIR/DNA/cov.counts.filtered.merged.exon.stats) | mawk '!/NC_007175.2/' > RnD.cov.stats
```
Calculate lower 10th percentile of exon sizes
```bash
mawk '{print $3 -$2}' RnD.cov.stats |sort -g | perl -e '$d=.1;@l=<>;print $l[int($d*@l)]'
```
Resut: `59`

Calculate upper 10th percentile of exon sizes
```bash
mawk '{print $3 -$2}' RnD.cov.stats |sort -g | perl -e '$d=.9;@l=<>;print $l[int($d*@l)]'
```
Result: `517`

Mark exons into size classes based on size distribution and create data table
```bash
mawk '{if ( $3 -$2 > 517 ) print $0 "\tUpper"; else if ( $3 - $2 < 59 ) print $0 "\tLower"; else if ( $3 - $2 > 58 && $3 - $2 < 518) print $0 "\tMiddle" }' RnD.cov.stats > rnd.cov.stats.class
echo -e "Chrom\tStart\tEnd\tRNA_Coverage\tDNA_Coverage\tExon_Size_Class" > header
cat header rnd.cov.stats.class > ExonCoverage.txt
```

## R code for Figure 3

```R
library(MASS)
library(fields)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)

TotalExon <- read.table("./ExonCoverage.txt", header = TRUE)
TotalExon <-as.data.frame(TotalExon)

TotalExon$Exon_Size_Class <-factor(TotalExon$Exon_Size_Class, levels=c("Lower","Middle","Upper"))

TotalExon$Exon_Size_Class <- revalue(TotalExon$Exon_Size_Class, c("Lower"="Lower 10%", "Upper"="Upper 10%", "Middle"="Middle 80%"))

TX <- TotalExon
TX$RNA_Coverage <- log(TX$RNA_Coverage +1)
TX$DNA_Coverage <- log(TX$DNA_Coverage + 1)


TotalExon$density <- interp.surface(kde2d(TX$RNA_Coverage, TX$DNA_Coverage), TX[,c("RNA_Coverage", "DNA_Coverage")])

cbPalette <- c("#009E73","#D55E00","#56B4E9", "#0072B2","#E69F00", "#F0E442" , "#999999","#CC79A7")
t <- cbPalette[1]
cbPalette[1] <- cbPalette[2]
cbPalette[2] <- t

b <- ggplot(TotalExon, aes(x=RNA_Coverage+1,y=DNA_Coverage+1,alpha = 1/(density)),fill=Exon_Size_Class,color=Exon_Size_Class) + 
  geom_abline(intercept =0, slope =1) +
  geom_point(aes(color=TotalExon$Exon_Size_Class,fill=TotalExon$Exon_Size_Class,shape=TotalExon$Exon_Size_Class)) +
  scale_alpha_continuous(guide = "none",range = c(.3, .99)) + 
  scale_shape_manual(values=c(15,16,17), name="Exon Size Percentile") +
  scale_fill_manual(values=cbPalette , name="Exon Size Percentile")+
  scale_color_manual(values=cbPalette, name="Exon Size Percentile") +
  scale_x_log10(limits=c(1,150000),expand=c(0.02,0), breaks = c(0,1,10,100,1000,10000,100000),
                labels = c("0","0","10","100","1,000","10,000","100,000"))+
  scale_y_log10(limits=c(1,150000),expand=c(0.02,0), breaks = c(0,1,10,100,1000,10000,100000),
                labels = c("0","0","10","100","1,000","10,000","100,000"))+
  xlab("Total Number of RNA Reads per Exon")+
  ylab("Total Number of DNA Reads per Exon") +
  theme_bw() +
  theme(legend.position = c(0.85,0.25)) 

png(filename="Figure3.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
b
dev.off()    
```
## Calculating Sensitivity

First, we need to return to the RNA folder
```bash
cd $WORKING_DIR/RNA
```
Next step is to find exons with minimum thresholds of RNA coverage.  These will be our "target" sets along with confidence intervals.  Based on overal RNA coverage, we chose 35X as our "target" set and choose 15X boundaries around that number for confidence intervals.  We will create three `bed` files from our RNA exon coverage stats.

```bash
mawk 'BEGIN { FS = "\t" } ; $11 > 19' $WORKING_DIR/RNA/cov.m4q4.EiR.stats > m4q4.EiRc20.bed
mawk 'BEGIN { FS = "\t" } ; $11 > 49' $WORKING_DIR/RNA/cov.m4q4.EiR.stats > m4q4.EiRc50.bed
mawk 'BEGIN { FS = "\t" } ; $11 > 34' $WORKING_DIR/RNA/cov.m4q4.EiR.stats > m4q4.EiRc35.bed
```
### Calculating data for table 3

We will use a BASH function to automate this for us:

```bash
counts_per_target(){

#Calculate number of exons with more than 1X coverage
EXONC=$(bedtools coverage -b $1.F.bam -a sorted.ref3.0.exon.sc.bed -counts -sorted -g genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 35X targets with more than 1X coverage
X35XC=$(bedtools coverage -b $1.F.bam -a m4q4.EiRc35.bed -counts -sorted -g genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 20X targets with more than 1X coverage
X20XC=$(bedtools coverage -b $1.F.bam -a m4q4.EiRc20.bed -counts -sorted -g genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 50X targets with more than 1X coverage
X50XC=$(bedtools coverage -b $1.F.bam -a m4q4.EiRc50.bed -counts -sorted -g genome.file | mawk '$4 > 0' | wc -l) 

#Calculate the total number of targets for each set
EXON=$(cat sorted.ref3.0.exon.sc.bed | wc -l )
X35X=$(cat m4q4.EiRc35.bed | wc -l)
X20X=$(cat m4q4.EiRc20.bed | wc -l)
X50X=$(cat m4q4.EiRc50.bed | wc -l)

#Print results in pretty percentages
echo $1
echo `python -c "print(round("$EXONC"/"$EXON" * 100,1))"`"%"
echo `python -c "print(round("$X20XC"/"$X20X" * 100,1))"`"%"
echo `python -c "print(round("$X35XC"/"$X35X" * 100,1))"`"%"
echo `python -c "print(round("$X50XC"/"$X50X" * 100,1))"`"%"
  
}

export -f counts_per_target
```

Now with this function we can use `paste` and subshells to produce the table
```bash
paste <(echo -e "Targets\nAll Exons\n20XR Exons\n35XR Exons\n50XR Exons") <(counts_per_target ECI_2) <(counts_per_target ECI_4) <(counts_per_target ECI_7) <(counts_per_target ECI_1) <(counts_per_target ECI_3) <(counts_per_target ECI_12) > Table3.txt
```
View the results:

`cat Table3`
```bash
Targets	ECI_2	ECI_4	ECI_7	ECI_1	ECI_3	ECI_12
All Exons	88.0%	86.0%	85.8%	86.5%	87.9%	86.4%
20XR Exons	99.5%	99.4%	99.4%	99.4%	99.5%	99.4%
35XR Exons	99.6%	99.6%	99.6%	99.6%	99.6%	99.6%
50XR Exons	99.7%	99.7%	99.7%	99.7%	99.7%	99.7%
```

### Calculating data for supplemental table 5

This is effectively the same as table 3 with a simple change to the coverage variable

```bash
counts_per_target10(){

#Calculate number of exons with more than 1X coverage
EXONC=$(bedtools coverage -b $1.F.bam -a sorted.ref3.0.exon.sc.bed -counts -sorted -g genome.file | mawk '$4 > 9' | wc -l) 
#Calculate number of 35X targets with more than 1X coverage
X35XC=$(bedtools coverage -b $1.F.bam -a m4q4.EiRc35.bed -counts -sorted -g genome.file | mawk '$4 > 9' | wc -l) 
#Calculate number of 20X targets with more than 1X coverage
X20XC=$(bedtools coverage -b $1.F.bam -a m4q4.EiRc20.bed -counts -sorted -g genome.file | mawk '$4 > 9' | wc -l) 
#Calculate number of 50X targets with more than 1X coverage
X50XC=$(bedtools coverage -b $1.F.bam -a m4q4.EiRc50.bed -counts -sorted -g genome.file | mawk '$4 > 9' | wc -l) 

#Calculate the total number of targets for each set
EXON=$(cat sorted.ref3.0.exon.sc.bed | wc -l )
X35X=$(cat m4q4.EiRc35.bed | wc -l)
X20X=$(cat m4q4.EiRc20.bed | wc -l)
X50X=$(cat m4q4.EiRc50.bed | wc -l)

#Print results in pretty percentages
echo $1
echo `python -c "print(round("$EXONC"/"$EXON" * 100,1))"`"%"
echo `python -c "print(round("$X20XC"/"$X20X" * 100,1))"`"%"
echo `python -c "print(round("$X35XC"/"$X35X" * 100,1))"`"%"
echo `python -c "print(round("$X50XC"/"$X50X" * 100,1))"`"%"
  
}

export -f counts_per_target10
```

Now with this function we can use `paste` and subshells to produce the table
```bash
paste <(echo -e "Targets\nAll Exons\n20XR Exons\n35XR Exons\n50XR Exons") <(counts_per_target10 ECI_2) <(counts_per_target10 ECI_4) <(counts_per_target10 ECI_7) <(counts_per_target10 ECI_1) <(counts_per_target10 ECI_3) <(counts_per_target10 ECI_12) > SupTable5.txt
```

`cat SupTable5.txt`

```bash
Targets	ECI_2	ECI_4	ECI_7	ECI_1	ECI_3	ECI_12
All Exons	61.3%	54.4%	54.8%	55.1%	61.1%	56.7%
20XR Exons	98.5%	97.2%	97.5%	97.1%	98.4%	97.9%
35XR Exons	99.2%	98.8%	98.9%	98.7%	99.1%	99.0%
50XR Exons	99.4%	99.1%	99.2%	99.1%	99.3%	99.2%
```

## Generate data for Figure 4

Figure 4 is per bp sensitivity looking at coverage across our target sets, near targets (definied as 150 bp around the edge of targets, and off target (everything that is not near or on target).

First steps involve creating our different interval sets using bedtools.

```bash
bedtools flank -i m4q4.EiRc35.bed -b 150 -g genome.file | bedtools sort -faidx genome.file >  m4q4.EiRc35.neartarget.bed 
bedtools slop -i m4q4.EiRc35.bed -b 150 -g genome.file > m4q4.EiRc35.slop.bed 
bedtools complement -i m4q4.EiRc35.slop.bed -g genome.file > m4q4.EiRc35.offtarget.bed 

bedtools flank -i m4q4.EiRc20.bed -b 150 -g genome.file | bedtools sort -faidx genome.file >  m4q4.EiRc20.neartarget.bed 
bedtools slop -i m4q4.EiRc20.bed -b 150 -g genome.file > m4q4.EiRc20.slop.bed 
bedtools complement -i m4q4.EiRc20.slop.bed -g genome.file > m4q4.EiRc20.offtarget.bed 

bedtools flank -i m4q4.EiRc50.bed -b 150 -g genome.file | bedtools sort -faidx genome.file >  m4q4.EiRc50.neartarget.bed
bedtools slop -i m4q4.EiRc50.bed -b 150 -g genome.file > m4q4.EiRc50.slop.bed
bedtools complement -i m4q4.EiRc50.slop.bed -g genome.file > m4q4.EiRc50.offtarget.bed
```
With the target sets defined we again use bedtools to caculate coverage levels across these various genomic regions, and below we use GNU-parallel to speed things up.

```bash
ls *F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a m4q4.EiRc35.neartarget.bed -sorted -g genome.file  | grep ^all > {}.hist.EiRc35neartarget.all.txt' 
ls *F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a m4q4.EiRc35.offtarget.bed -sorted -g genome.file  | grep ^all > {}.hist.EiRc35nofftarget.all.txt' 
ls *F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a m4q4.EiRc35.bed -sorted -g genome.file  | grep ^all > {}.hist.EiRc35.all.txt' 

ls *F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a m4q4.EiRc20.neartarget.bed -sorted -g genome.file  | grep ^all > {}.hist.EiRc20neartarget.all.txt' 
ls *F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a m4q4.EiRc20.offtarget.bed -sorted -g genome.file  | grep ^all > {}.hist.EiRc20nofftarget.all.txt' 
ls *F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a m4q4.EiRc20.bed -sorted -g genome.file  | grep ^all > {}.hist.EiRc20.all.txt' 

ls *F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a m4q4.EiRc50.neartarget.bed -sorted -g genome.file  | grep ^all > {}.hist.EiRc50neartarget.all.txt' 
ls *F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a m4q4.EiRc50.offtarget.bed -sorted -g genome.file  | grep ^all > {}.hist.EiRc50nofftarget.all.txt' 
ls *F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a m4q4.EiRc50.bed -sorted -g genome.file  | grep ^all > {}.hist.EiRc50.all.txt' 
```
## R code to generate Figure 4

```R
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

make_graph= function(j){
  j1 <- gsub("ECI", "EC", j, ignore.case=FALSE, fixed=FALSE)
  print(files <- list.files(pattern=paste(j,".hist.*EiRc35*", sep = "")))
  labs <- c("On target","Near target", "Off target")
  
  cov <- list()
  for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])[,c(2,5)]
    cov_cumul=1-cumsum(cov[[i]][,2])
    cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
    cov[[i]]$sample=labs[i]
  }
  
  cov_df=do.call("rbind",cov)
  names(cov_df)[1:2]=c("depth","fraction")
  
  print(files <- list.files(pattern=paste(j,".hist.*EiRc20*", sep = "")))
  cov2 <- list()
  for (i in 1:length(files)) {
    cov2[[i]] <- read.table(files[i])[,c(2,5)]
    cov_cumul2=1-cumsum(cov2[[i]][,2])
    cov2[[i]]$cov_cumul <- c(1,cov_cumul2[-length(cov_cumul2)])
    cov2[[i]]$sample=labs[i]
  }
  cov2_df=do.call("rbind",cov2)
  names(cov2_df)[1:2]=c("depth","fraction")
  
  print(files <- list.files(pattern=paste(j,".hist.*EiRc50*", sep = "")))
  
  cov3 <- list()
  for (i in 1:length(files)) {
    cov3[[i]] <- read.table(files[i])[,c(2,5)]
    cov_cumul3=1-cumsum(cov3[[i]][,2])
    cov3[[i]]$cov_cumul <- c(1,cov_cumul3[-length(cov_cumul3)])
    cov3[[i]]$sample=labs[i]
  }
  cov3_df=do.call("rbind",cov3)
  names(cov3_df)[1:2]=c("depth","fraction")
  
  cov_df <- subset(cov_df, depth <101)
  cov2_df <- subset(cov2_df, depth <101)
  cov3_df <- subset(cov3_df, depth <101)
  
  cov2_df$high <-cov3_df$cov_cumul
  
  cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7")
  cbPalettep <- c("#0072B2", "#009E73","#D55E00" )
  
  cov_df$sample <- factor(cov_df$sample,levels=c("On target","Near target", "Off target"))
  cov2_df$sample <- factor(cov2_df$sample,levels=c("On target","Near target", "Off target"))
  
  p <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  +
    #xlim(0,200)+
    geom_ribbon(data=cov2_df,aes(ymin=cov_cumul,ymax=high, color=sample, fill=sample, alpha=0.4)) +
    scale_alpha(guide = 'none') +
    geom_line(size=1.5)+ 
    scale_color_manual(values=cbPalettep) +
    scale_fill_manual(values=cbPalettep) +
    ylab("% of Target Bases > Depth")+
    #scale_x_log10()+
    xlab("Depth")+
    ggtitle(eval(j1)) +
    theme_bw() +
    theme(plot.title = element_text(size = 10, face = "bold",hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.position="none")
  #theme(legend.position=c(0.92,0.88))
  
  return(p)
}

sample_names=c("ECI_2","ECI_4","ECI_7","ECI_1","ECI_3","ECI_12")

ECI_2 <- make_graph(sample_names[1])

ECI_4 <- make_graph(sample_names[2])
ECI_4 <- ECI_4 + theme(axis.title.y=element_text(color="transparent"))

ECI_7 <- make_graph(sample_names[3])
ECI_7 <- ECI_7 + theme(axis.title.y=element_text(color="transparent"))

ECI_1 <- make_graph(sample_names[4])

ECI_3 <- make_graph(sample_names[5])
ECI_3 <- ECI_3 + theme(axis.title.y=element_text(color="transparent"))

ECI_12 <- make_graph(sample_names[6])
ECI_12 <- ECI_12 + theme(axis.title.y=element_text(color="transparent"))



pdf(file="Figure4.pdf",width=14, height=6.5, bg="transparent")
multiplot(ECI_2,ECI_1,ECI_4,ECI_3,ECI_7,ECI_12, cols=3)

dev.off()

pdf(file="Figure4Legend.pdf",width=14, height=6.5, bg="transparent")

ECI_4 <- ECI_4 + theme(legend.position="bottom")
ECI_4
dev.off()
```

## Calculating Specificity

First let's calculate near and off-target intervals for all exons
```bash
bedtools flank -i sorted.ref3.0.exon.sc.bed -b 150 -g genome.file | bedtools sort -faidx genome.file >  sorted.ref3.0.exon.sc.neartarget.bed
bedtools slop -i sorted.ref3.0.exon.sc.bed -b 150 -g genome.file > sorted.ref3.0.exon.sc.slop.bed
bedtools complement -i sorted.ref3.0.exon.sc.slop.bed -g genome.file > sorted.ref3.0.exon.sc.offtarget.bed
```
Now we can create a specificity table for all exons and for expressed targets using a few more BASH functions

```bash
specExon(){

exon_reads=$(samtools view -@32 $1.F.bam -c -L sorted.ref3.0.exon.sc.bed)
exon_nearr=$(samtools view -@32 $1.F.bam -c -L sorted.ref3.0.exon.sc.neartarget.bed)
exon_nearo=$(samtools view $1.F.bam  -h -@32 -L sorted.ref3.0.exon.sc.bed | samtools view - -@32 -c -L sorted.ref3.0.exon.sc.neartarget.bed)
exon_offtr=$(samtools view -@32 $1.F.bam -c -L sorted.ref3.0.exon.sc.offtarget.bed)
exon_nearO=$(samtools view $1.F.bam  -h -@32 -L sorted.ref3.0.exon.sc.slop.bed | samtools view - -@32 -c -L sorted.ref3.0.exon.sc.offtarget.bed)
total=$(samtools view -@32 $1.F.bam -c)


echo -e $1"\t"`python -c "print(round("$exon_reads"/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"%\t"
  
}

export -f specExon
```
```bash
spec35X(){

exon_reads=$(samtools view -@32 $1.F.bam -c -L m4q4.EiRc35.bed)
exon_nearr=$(samtools view -@32 $1.F.bam -c -L m4q4.EiRc35.neartarget.bed )
exon_nearo=$(samtools view $1.F.bam  -h -@32 -L m4q4.EiRc35.bed | samtools view - -@32 -c -L m4q4.EiRc35.neartarget.bed )
exon_offtr=$(samtools view -@32 $1.F.bam -c -L m4q4.EiRc35.offtarget.bed)
exon_nearO=$(samtools view $1.F.bam  -h -@32 -L m4q4.EiRc35.neartarget.bed  | samtools view - -@32 -c -L m4q4.EiRc35.offtarget.bed)
total=$(samtools view -@32 $1.F.bam -c)


echo -e `python -c "print(round("$exon_reads"/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"%\t"
  
}

export -f spec35X
```
```bash
spec20X(){

exon_reads=$(samtools view -@32 $1.F.bam -c -L m4q4.EiRc20.bed)
exon_nearr=$(samtools view -@32 $1.F.bam -c -L m4q4.EiRc20.neartarget.bed )
exon_nearo=$(samtools view $1.F.bam  -h -@32 -L m4q4.EiRc20.bed | samtools view - -@32 -c -L m4q4.EiRc20.neartarget.bed )
exon_offtr=$(samtools view -@32 $1.F.bam -c -L m4q4.EiRc20.offtarget.bed)
exon_nearO=$(samtools view $1.F.bam  -h -@32 -L m4q4.EiRc20.neartarget.bed  | samtools view - -@32 -c -L m4q4.EiRc20.offtarget.bed)
total=$(samtools view -@32 $1.F.bam -c)


echo -e `python -c "print(round("$exon_reads"/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"%\t"
  
}

export -f spec20X
```
```bash
spec50X(){

exon_reads=$(samtools view -@32 $1.F.bam -c -L m4q4.EiRc50.bed)
exon_nearr=$(samtools view -@32 $1.F.bam -c -L m4q4.EiRc50.neartarget.bed )
exon_nearo=$(samtools view $1.F.bam  -h -@32 -L m4q4.EiRc50.bed | samtools view - -@32 -c -L m4q4.EiRc50.neartarget.bed )
exon_offtr=$(samtools view -@32 $1.F.bam -c -L m4q4.EiRc50.offtarget.bed)
exon_nearO=$(samtools view $1.F.bam  -h -@32 -L m4q4.EiRc50.neartarget.bed  | samtools view - -@32 -c -L m4q4.EiRc50.offtarget.bed)
total=$(samtools view -@32 $1.F.bam -c)


echo -e `python -c "print(round("$exon_reads"/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"%\t"
  
}
export -f spec50X
```
Now we use all the functions to create a table
```bash
echo -e "Pool\t%_in_Exons\t%_Near_Exons\t%Off_Target_Exons\t\t%_on_Target35X\t%_Near_Target35X\t%Off_Target35X\t\t%_on_Target20X\t%_Near_Target20X\t%Off_Target20X\t\t%_on_Target50X\t%_Near_Target50X\t%Off_Target50X" > Spec.Table

paste <(specExon ECI_2) <(spec35X ECI_2) <(spec20X ECI_2) <(spec50X ECI_2) >> Spec.Table
paste <(specExon ECI_4) <(spec35X ECI_4) <(spec20X ECI_4) <(spec50X ECI_4) >> Spec.Table
paste <(specExon ECI_7) <(spec35X ECI_7) <(spec20X ECI_7) <(spec50X ECI_7) >> Spec.Table
paste <(specExon ECI_1) <(spec35X ECI_1) <(spec20X ECI_1) <(spec50X ECI_1) >> Spec.Table
paste <(specExon ECI_3) <(spec35X ECI_3) <(spec20X ECI_3) <(spec50X ECI_3) >> Spec.Table
paste <(specExon ECI_12) <(spec35X ECI_12) <(spec20X ECI_12) <(spec50X ECI_12) >> Spec.Table
```
`cat Spec.Table`
```bash
Pool	%_in_Exons	%_Near_Exons	%Off_Target_Exons		%_on_Target35X	%_Near_Target35X	%Off_Target35X		%_on_Target20X	%_Near_Target20X	%Off_Target20X		%_on_Target50X	%_Near_Target50X	%Off_Target50X
ECI_2	49.1%	6.7%	44.2%		37.9%	3.5%	58.6%		42.3%	4.3%	53.3%		34.3%	3.0%	62.7%	
ECI_4	43.6%	6.9%	49.5%		33.3%	3.5%	63.1%		37.4%	4.4%	58.2%		30.1%	3.0%	67.0%	
ECI_7	48.6%	6.7%	44.6%		38.1%	3.6%	58.3%		42.3%	4.4%	53.2%		34.6%	3.1%	62.3%	
ECI_1	45.2%	7.2%	47.6%		34.5%	3.5%	62.0%		38.6%	4.4%	57.0%		31.2%	3.0%	65.9%	
ECI_3	48.8%	6.7%	44.5%		37.9%	3.6%	58.6%		42.3%	4.3%	53.4%		34.3%	3.0%	62.7%	
ECI_12	51.8%	6.6%	41.6%		40.7%	3.6%	55.7%		45.3%	4.4%	50.3%		36.9%	3.1%	60.0%
```

## Generating data for Figure 5

First, we need to break up exons into size classes.

```bash
mawk '$3 - $2 > 59' sorted.ref3.0.exon.sc.bed | mawk '$3 - $2 < 181' > exon_bin60to180.bed
mawk '$3 - $2 > 180' sorted.ref3.0.exon.sc.bed | mawk '$3 - $2 < 241' > exon_bin180to240.bed
mawk '$3 - $2 > 240' sorted.ref3.0.exon.sc.bed | mawk '$3 - $2 < 421' > exon_bin240to420.bed
mawk '$3 - $2 > 420' sorted.ref3.0.exon.sc.bed | mawk '$3 - $2 < 541' > exon_bin420to540.bed
```
Next, we use bedtools to divide each exon size class into discrete windows and remove mtDNA.
```bash
bedtools makewindows -b exon_bin60to180.bed -w 10 -i srcwinnum | mawk '$1 !~ /NC_007175.2/' > exon_bin60to180b10.bed 
bedtools makewindows -b exon_bin180to240.bed -w 15 -i srcwinnum| mawk '$1 !~ /NC_007175.2/' > exon_bin180to240b15.bed 
bedtools makewindows -b exon_bin240to420.bed -w 28 -i srcwinnum| mawk '$1 !~ /NC_007175.2/' > exon_bin240to420b28.bed 
bedtools makewindows -b exon_bin420to540.bed -w 30 -i srcwinnum| mawk '$1 !~ /NC_007175.2/' > exon_bin420to540b30.bed 
```

The next step is to use bedtools to calculate coverage across these descrete windows.
```bash
cat <(echo -e "Chrom\tStart\tEnd\tBin\tDepth") <(bedtools coverage -b filter.merged.bam -a exon_bin60to180b10.bed -mean -sorted -g genome.file) > Allcap.mean.exon_bin60to180b10.cov &
cat <(echo -e "Chrom\tStart\tEnd\tBin\tDepth") <(bedtools coverage -b filter.merged.bam -a exon_bin180to240b15.bed -mean -sorted -g genome.file) > Allcap.mean.exon_bin180to240b15.cov &
cat <(echo -e "Chrom\tStart\tEnd\tBin\tDepth") <(bedtools coverage -b filter.merged.bam -a exon_bin240to420b28.bed -mean -sorted -g genome.file) > Allcap.mean.exon_bin240to420b28.cov &
cat <(echo -e "Chrom\tStart\tEnd\tBin\tDepth") <(bedtools coverage -b filter.merged.bam -a exon_bin420to540b30.bed -mean -sorted -g genome.file) > Allcap.mean.exon_bin420to540b30.cov 
```

Finally, the bin numbering won't sort properly in R, so this short script will substitute the proper number
`echo sed.sh`

```bash
!#/bin/bash

sed 's/_1\t/_01\t/g' $1 | sed 's/_2\t/_02\t/g' | sed 's/_3\t/_03\t/g' | sed 's/_4\t/_04\t/g' | sed 's/_5\t/_05\t/g' | sed 's/_6\t/_06\t/g' | sed 's/_7\t/_07\t/g' | sed 's/_8\t/_08\t/g' | sed 's/_9\t/_09\t/g' > $1.1
mv $1.1 $1
```
```bash
ls Allcap.mean.exon_bin*b*.cov | parallel bash sed.sh {}
```
## R code for Figure 5
```R
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
#library(tidyquant)
library(zoo)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


Exon60to180 <- read.table("Allcap.mean.exon_bin60to180b10.cov", header = TRUE)

Exon60to180 <- as.data.frame(Exon60to180)
#Exon60to180$Depth <- Exon60to180$Depth/6


f <- function(x) {
  ans <- boxplot.stats(x)
  data.frame(ymin = ans$conf[1], ymax = ans$conf[2])
}


plot1<- ggplot(Exon60to180, aes(x=Bin, y=Depth)) + 
  stat_summary(fun.y="median", geom="bar", fill=cbPalette[1], alpha=0.9) +
  stat_summary(fun.data = f, geom = "errorbar") +
  theme_classic() +
  ylab("Median Coverage") +
  xlab("10 bp Windows Across Exons 60bp - 180bp in Length")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  coord_cartesian(ylim = c(0,45), expand = c(0,0))


Exon180to240 <- read.table("Allcap.mean.exon_bin180to240b15.cov", header = TRUE)

Exon180to240 <- as.data.frame(Exon180to240)
#Exon180to240$Depth <- Exon180to240$Depth/6


plot2<- ggplot(Exon180to240, aes(x=Bin, y=Depth)) + 
  stat_summary(fun.y="median", geom="bar", fill=cbPalette[2], alpha=0.9) +
  stat_summary(fun.data = f, geom = "errorbar") +
  #geom_abline(intercept = mean(Exon180to240$Depth), slope =0) +
  #scale_y_continuous(expand = c(0,0))+
  theme_classic() +
  ylab("Median Coverage") +
  xlab("15 bp Windows Across Exons 180bp - 240bp in Length")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  coord_cartesian(ylim = c(0,45), expand = c(0,0))

Exon240to420 <- read.table("Allcap.mean.exon_bin240to420b28.cov", header = TRUE)

Exon240to420 <- as.data.frame(Exon240to420)
#Exon240to420$Depth <- Exon240to420$Depth/6


plot3<- ggplot(Exon240to420, aes(x=Bin, y=Depth)) + 
  stat_summary(fun.y="median", geom="bar", fill=cbPalette[4], alpha=0.9) +
  stat_summary(fun.data = f, geom = "errorbar") +
  #geom_abline(intercept = mean(Exon240to420$Depth), slope =0) +
  #scale_y_continuous(expand = c(0,0))+
  theme_classic() +
  ylab("Median Coverage") +
  xlab("28 bp Windows Across Exons 240bp - 420bp in Length")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  coord_cartesian(ylim = c(0,45), expand = c(0,0))

Exon420to540 <- read.table("Allcap.mean.exon_bin420to540b30.cov", header = TRUE)

Exon420to540 <- as.data.frame(Exon420to540)
#Exon420to540$Depth <- Exon420to540$Depth/6


plot4<- ggplot(Exon420to540, aes(x=Bin, y=Depth)) + 
  stat_summary(fun.y="median", geom="bar", fill=cbPalette[5], alpha=0.9) +
  stat_summary(fun.data = f, geom = "errorbar") +
  #geom_abline(intercept = mean(Exon420to540$Depth), slope =0,linetype=dash) +
  #scale_y_continuous(expand = c(0,0))+
  theme_classic() +
  ylab("Median Coverage") +
  xlab("30 bp Windows Across Exons 420bp - 540bp in Length")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  coord_cartesian(ylim = c(0,45), expand = c(0,0))


png(filename="Figure5.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
multiplot(plot1, plot3,plot2,plot4, cols=2)
dev.off()
```

## Code to generate data for Figure 6

First, we use bedtools to calculate the GC content for each exon

```bash
bedtools nuc -fi reference.fasta -bed sorted.ref3.0.exon.sc.bed > gc.exon.stats
```

Next we can add the GC content as a column to a previous data table
paste ExonCoverage.txt <(mawk '!/NC_007175.2/' gc.exon.stats| cut -f5 ) | sed 's/5_pct_gc/GC_Content/' > ExonCoverageGC.txt

## R code to generate Figure 6
```R
library(MASS)
library(fields)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

cbPalette <- c("#009E73","#D55E00","#56B4E9", "#0072B2","#E69F00", "#F0E442" , "#999999","#CC79A7")


TotalExon <- read.table("./ExonCoverageGC.txt", header = TRUE)
TotalExon <-as.data.frame(TotalExon)

TotalExon$Exon_Size_Class <-factor(TotalExon$Exon_Size_Class, levels=c("Lower","Middle","Upper"))
TotalExon$Exon_Size_Class <- revalue(TotalExon$Exon_Size_Class, c("Lower"="Lower 10%", "Upper"="Upper 10%", "Middle"="Middle 80%"))

TX <- TotalExon
TX$RNA_Coverage <- log(TX$RNA_Coverage +1)
TX$DNA_Coverage <- log(TX$DNA_Coverage/6 + 1)

TotalExon$density2 <- interp.surface(kde2d(TX$GC_Content, TotalExon$DNA_Coverage), TX[,c("GC_Content", "DNA_Coverage")])
TotalExon$density2 <- TotalExon$density2 +0.000001


gc1 <- ggplot(data=TotalExon[which(TotalExon$Exon_Size_Class=="Middle 80%"),], aes(x=GC_Content,y=DNA_Coverage/6+1)) + 
  geom_point(color=cbPalette[2],fill=cbPalette[2],shape=15, alpha=0.1)+
  xlim(0,1)+
c("0","0","10","100","1,000","10,000","100,000"))+
  scale_y_log10(limits=c(1,30000),expand=c(0.02,0), breaks = c(0,1,10,100,1000,10000),labels = c("0","0","10","100","1,000","10,000"))+
  xlab("% of GC Content")+
  ylab("Mean Capture Depth") +
  theme_bw() +
  theme(legend.position = c(0.85,0.25)) 
gc2 <- ggplot(TotalExon, aes(x=GC_Content,y=DNA_Coverage/6+1)) + 
  geom_point(data=TotalExon[which(TotalExon$Exon_Size_Class=="Lower 10%"),], color=cbPalette[1],fill=cbPalette[1],shape=16, alpha=0.1)+
   xlim(0,1)+
c("0","0","10","100","1,000","10,000","100,000"))+
  scale_y_log10(limits=c(1,30000),expand=c(0.02,0), breaks = c(0,1,10,100,1000,10000),labels = c("0","0","10","100","1,000","10,000"))+
  xlab("% of GC Content")+
  ylab("Mean Capture Depth") +
  theme_bw() +
  theme(legend.position = c(0.85,0.25)) 
gc3 <- ggplot(TotalExon, aes(x=GC_Content,y=DNA_Coverage/6+1)) + 
  xlim(0,1)+
  geom_point(data=TotalExon[which(TotalExon$Exon_Size_Class=="Upper 10%"),], color=cbPalette[3],fill=cbPalette[3],shape=17, alpha=0.1)+
  #scale_x_log10(limits=c(1,150000),expand=c(0.02,0), breaks = c(0,1,10,100,1000,10000,100000),labels = c("0","0","10","100","1,000","10,000","100,000"))+
  scale_y_log10(limits=c(1,30000),expand=c(0.02,0), breaks = c(0,1,10,100,1000,10000),labels = c("0","0","10","100","1,000","10,000"))+
  xlab("% of GC Content")+
  ylab("Mean Capture Depth") +
  theme_bw() +
  theme(legend.position = c(0.85,0.25)) 

png(filename="Figure6.png", type="cairo",units="px", width=5600, 
    height=1500, res=400, bg="transparent")
multiplot(gc2,gc1,gc3,cols=3)
dev.off() 
```

## Comparing coverage in CDS and UTR in exons with both

First we can find UTRs by subtracting CDS from exons using bedtools.
```bash
bedtools subtract -a sorted.ref3.0.exon.bed -b sorted.ref3.0.CDS.bed -g genome.file -sorted > sorted.ref3.0.UTR.bed
bedtools merge -i sorted.ref3.0.UTR.bed > sorted.ref3.0.UTR.sc.bed
```
Next step is to find exons that have CDS and UTR.
```bash
bedtools intersect -a sorted.ref3.0.exon.sc.bed -b sorted.ref3.0.UTR.sc.bed -wa > sort.ref3.0.exonwithUTR.sc.bed
bedtools intersect -a sort.ref3.0.exonwithUTR.sc.bed -b sorted.ref3.0.CDS.sc.bed -wa > sort.ref3.0.exonwithbothUTRandCDS.bed
```
Next, we create CDS only and UTR only regions from those exons with both UTRs and CDS
```bash
bedtools intersect -a sorted.ref3.0.CDS.sc.bed -b sort.ref3.0.exonwithbothUTRandCDS.bed -wa > sort.comp.CDS.bed
bedtools intersect -a sorted.ref3.0.UTR.sc.bed -b sort.ref3.0.exonwithbothUTRandCDS.bed -wa > sort.comp.UTR.bed
```
Now, we use bedtools to calculate coverage across all captures for both regions
```bash
bedtools coverage -b filter.merged.bam -a sort.comp.CDS.bed -mean -sorted -g genome.file > AllCAP.mean.CDSwithUTR.cov
bedtools coverage -b filter.merged.bam -a sort.comp.UTR.bed -mean -sorted -g genome.file > AllCAP.mean.UTRwithCDS.cov
```
Create a combined data file
```bash
cat <(echo -e "Chrom\tStart\tEnd\tCoverage\tClass" ) <( mawk '{print $0"\tCDS"}' AllCAP.mean.CDSwithUTR.cov ) <( mawk '{print $0"\tUTR"}' AllCAP.mean.UTRwithCDS.cov ) > UTRvsCDS.txt
```
### R code for test
```R
CvU <- read.table("./UTRvsCDS.txt", header = TRUE)
CvU <-as.data.frame(CvU)

UTR <- CvU[which(CvU$Class=="UTR"),]
CDS <- CvU[which(CvU$Class=="CDS"),]

UTR$Coverage <- UTR$Coverage/6
CDS$Coverage <- CDS$Coverage/6

t.test(CDS$Coverage,UTR$Coverage)
```
Results:
```R
	Welch Two Sample t-test

data:  CDS$Coverage and UTR$Coverage
t = 40.063, df = 135580, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 5.822934 6.421991
sample estimates:
mean of x mean of y 
 17.70905  11.58659 
```

Now change to RNA directory

`cd $Working_DIR/RNA`

Repeat bedtools steps
```bash
bedtools subtract -a sorted.ref3.0.exon.bed -b sorted.ref3.0.CDS.bed -g genome.file -sorted > sorted.ref3.0.UTR.bed
bedtools merge -i sorted.ref3.0.UTR.bed > sorted.ref3.0.UTR.sc.bed
bedtools intersect -a sorted.ref3.0.exon.sc.bed -b sorted.ref3.0.UTR.sc.bed -wa > sort.ref3.0.exonwithUTR.sc.bed
bedtools intersect -a sort.ref3.0.exonwithUTR.sc.bed -b sorted.ref3.0.CDS.sc.bed -wa > sort.ref3.0.exonwithbothUTRandCDS.bed
bedtools intersect -a sorted.ref3.0.CDS.sc.bed -b sort.ref3.0.exonwithbothUTRandCDS.bed -wa > sort.comp.CDS.bed
bedtools intersect -a sorted.ref3.0.UTR.sc.bed -b sort.ref3.0.exonwithbothUTRandCDS.bed -wa > sort.comp.UTR.bed

bedtools coverage -b m4.q4.merged.bam -a sort.comp.UTR.bed -mean -sorted -g genome.file -split  > m4q4.mean.UTR.cov 
bedtools coverage -b m4.q4.merged.bam -a sort.comp.CDS.bed -mean -sorted -g genome.file -split > m4q4.mean.CDSwithUTR.cov

cat <(echo -e "Chrom\tStart\tEnd\tCoverage\tClass" ) <( mawk '{print $0"\tCDS"}' m4q4.mean.CDSwithUTR.cov ) <( mawk '{print $0"\tUTR"}' m4q4.mean.UTR.cov ) > UTRvsCDSRNA.txt

```
### R Code
```R
CvUR <- read.table("./UTRvsCDSRNA.txt", header = TRUE)
CvUR <-as.data.frame(CvUR)

UTRR <- CvUR[which(CvUR$Class=="UTR"),]
CDSR <- CvUR[which(CvUR$Class=="CDS"),]


t.test(CDSR$Coverage,UTRR$Coverage)
```
Result:
```
	Welch Two Sample t-test

data:  CDSR$Coverage and UTRR$Coverage
t = 22.677, df = 129300, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 4.931675 5.864832
sample estimates:
mean of x mean of y 
13.653255  8.255002 
```

## Generate data for Figure 7

Change back to DNA directory

`cd $WORKING_DIR/Genome`

The first step is to generate depth per bp data for each capture pool.
```bash
ls *.F.bam | sed 's/.F.bam//g' | parallel "samtools depth -aa {}.F.bam > {}.genome.depth"
```

Next, we extract the region of interest.
```bash
mawk '$1 ~ /NC_035780.1/ && $2 > 32736205' ECI_1.genome.depth | mawk '$2 < 32866205' > ECI_1.graph.depth 
mawk '$1 ~ /NC_035780.1/ && $2 > 32736205' ECI_2.genome.depth | mawk '$2 < 32866205' > ECI_2.graph.depth 
mawk '$1 ~ /NC_035780.1/ && $2 > 32736205' ECI_4.genome.depth | mawk '$2 < 32866205' > ECI_4.graph.depth 
mawk '$1 ~ /NC_035780.1/ && $2 > 32736205' ECI_3.genome.depth | mawk '$2 < 32866205' > ECI_3.graph.depth 
mawk '$1 ~ /NC_035780.1/ && $2 > 32736205' ECI_7.genome.depth | mawk '$2 < 32866205' > ECI_7.graph.depth 
mawk '$1 ~ /NC_035780.1/ && $2 > 32736205' ECI_12.genome.depth | mawk '$2 < 32866205' > ECI_12.graph.depth
```

Next, we add a column for the pool identifier and concatenate into single data table.
```bash
sed -i 's/$/\tEC_1/g' ECI_1.graph.depth &
sed -i 's/$/\tEC_2/g' ECI_2.graph.depth &
sed -i 's/$/\tEC_4/g' ECI_4.graph.depth &
sed -i 's/$/\tEC_3/g' ECI_3.graph.depth &
sed -i 's/$/\tEC_7/g' ECI_7.graph.depth &
sed -i 's/$/\tEC_12/g' ECI_12.graph.depth & 

echo -e "Contig\tbp\tDepth\tSample" > header
cat header ECI_*.graph.depth > TotalCov.txt
```

For the graph, we all need the annotations for the gene regions, the exons, and CDS
```bash
mawk '$1 ~ /NC_035780.1/ && $4 > 32736205' ~/j.puritz_remote/new_genome/GFF/ref_C_virginica-3.0_top_level.gff3 | mawk '$5 < 32866205' | mawk '$3 == "exon"' | cut -f1,4,5,9 | uniq -w 30 | sed 's/ID=.*product=//g' | sed 's/;trans.*//g' | sed 's/%.*//g' > exons
cat <(echo -e "Contig\tStart\tEnd\tTreatment") exons > exon.list

mawk '$1 ~ /NC_035780.1/ && $4 > 32736205' ~/j.puritz_remote/new_genome/GFF/ref_C_virginica-3.0_top_level.gff3 | mawk '$5 < 32866205' | mawk '$3 == "mRNA"' | cut -f1,4,5,9 | uniq -w 30 | sed 's/ID=.*product=//g' | sed 's/;trans.*//g' | sed 's/%.*//g' > genes
cat <(echo -e "Contig\tStart\tEnd\tTreatment") genes > genes.list


mawk '$1 ~ /NC_035780.1/ && $4 > 32736205' ~/j.puritz_remote/new_genome/GFF/ref_C_virginica-3.0_top_level.gff3 | mawk '$5 < 32866205' | mawk '$3 == "CDS"' | cut -f1,4,5,9 | uniq -w 30 | sed 's/ID=.*product=//g' | sed 's/;trans.*//g' | sed 's/%.*//g' > CDS
cat <(echo -e "Contig\tStart\tEnd\tTreatment") CDS > CDS.list
```
We also need to perform similar steps for the RNA data
```bash
cd $Working_DIR/RNA

samtools view -@64 -q4 -h -F 0x100 -F 0x400 ../RNA1m3Aligned.sortedByCoord.out.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | samtools view -b > RNA1.bam
samtools view -@64 -q4 -h -F 0x100 -F 0x400 ../RNA2m3Aligned.sortedByCoord.out.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | samtools view -b > RNA2.bam
samtools view -@64 -q4 -h -F 0x100 -F 0x400 ../RNA3m3Aligned.sortedByCoord.out.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | samtools view -b > RNA3.bam
samtools view -@64 -q4 -h -F 0x100 -F 0x400 ../RNA4m3Aligned.sortedByCoord.out.bam | mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | samtools view -b > RNA4.bam


bedtools coverage -b RNA1.bam -a <(echo -e "NC_035780.1\t32740000\t32755000") -d -split -sorted -g genome.file | mawk 'BEGIN { FS = "\t" }{t = $2+$4 -1; print $1"\t"t"\t"$5}'> RNA1.exon.cov.stats
bedtools coverage -b RNA2.bam -a <(echo -e "NC_035780.1\t32740000\t32755000") -d -split -sorted -g genome.file | mawk 'BEGIN { FS = "\t" }{t = $2+$4 -1; print $1"\t"t"\t"$5}' > RNA2.exon.cov.stats
bedtools coverage -b RNA3.bam -a <(echo -e "NC_035780.1\t32740000\t32755000") -d -split -sorted -g genome.file | mawk 'BEGIN { FS = "\t" }{t = $2+$4 -1; print $1"\t"t"\t"$5}' > RNA3.exon.cov.stats
bedtools coverage -b  RNA4.bam -a <(echo -e "NC_035780.1\t32740000\t32755000") -d -split -sorted -g genome.file | mawk 'BEGIN { FS = "\t" }{t = $2+$4 -1; print $1"\t"t"\t"$5}'> RNA4.exon.cov.stats

sed -i 's/$/\tRNA_1/g' RNA1.exon.cov.stats 
sed -i 's/$/\tRNA_2/g' RNA2.exon.cov.stats 
sed -i 's/$/\tRNA_3/g' RNA3.exon.cov.stats 
sed -i 's/$/\tRNA_4/g' RNA4.exon.cov.stats 

echo -e "Contig\tbp\tDepth" > header

cat header *.exon.cov.stats > TotalRNACov1.txt
cd $Working_DIR/DNA
ln -s $Working_DIR/RNA/TotalRNACov1.txt
```

## R code for Figure 7
```R
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)

cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")
DepC <- read.table("TotalCov.txt", header = TRUE)
DepR <- read.table("TotalRNACov1.txt", header = TRUE)

DepC <- as.data.frame(DepC)
DepC$Sample <- factor(DepC$Sample,levels=c("EC_2","EC_4","EC_7","EC_1","EC_3","EC_12"))
DepR <- as.data.frame(DepR)
DepR$Sample <- factor(DepR$Sample,levels=c("RNA_1","RNA_2","RNA_3","RNA_4"))

exons <- read.table("exon.list", header = TRUE, sep = "\t")
exons <- as.data.frame(exons)

genes <- read.table("genes.list", header = TRUE, sep = "\t")
genes <- as.data.frame(genes)

cds <- read.table("cds.list", header = TRUE, sep = "\t")
cds <- as.data.frame(cds)

subDepC <-subset(DepC, bp <32755000 & bp > 32739000)
subDepR <-subset(DepR, bp <32755000 & bp > 32739000)
subexons <-subset(exons, End <32755205 & End > 32740205)
subgenes <-subset(genes, End <32800757 & Start < 32754201)
subcds <-subset(cds, End <32800757 & Start < 32755000)
subDepR$Depth <- subDepR$Depth / -1
submean.cov <- ddply(subDepC, .(Contig,bp), summarize,  Depth=mean(Depth))
submeanR.cov <- ddply(subDepR, .(Contig,bp), summarize,  Depth=mean(Depth))
subgenes$End[4] <- 32755000


pie(rep(1, length(cbPalette)), labels = sprintf("%d (%s)", seq_along(cbPalette), 
                                                cbPalette), col = cbPalette)
redcol <-"#940000"

cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")
cbPalettedd <- c( "#009E73","#D55E00", "#E69F00")


dd <- ggplot(subDepC, aes(x= bp, y=Depth)) +
  geom_area(aes(group=Sample),position = "identity",color=alpha("grey30",0.25),fill=cbPalette[4], alpha=0.1, linetype="dotted")+  
  geom_line(data=submean.cov,aes(y=rollmean(Depth, 100, na.pad=TRUE)),colour=cbPalette[4], size =1.0, alpha=0.9)  +
  geom_line(data=submeanR.cov,aes(y=rollmean(Depth, 100, na.pad=TRUE)),colour=redcol, size =1.0, alpha=0.9)  +
  geom_area(data=subDepR, aes(group=Sample),position = "identity",color=alpha("grey30",0.25),fill=redcol, alpha=0.1, linetype="dotted")+
  scale_color_manual(values=cbPalettedd) +
  geom_segment(data=subgenes, aes(x = Start, y = 715, xend = End, yend = 715), size = 6,color=cbPalette[9], alpha=1)+
  geom_segment(data=subexons,aes(x = Start, y = 715, xend = End, yend = 715, color=Treatment),size = 4, alpha=1) +
  geom_segment(data=subcds,aes(x = Start, y = 715, xend = End, yend = 715),size = 1, color="grey90", alpha=1) +
  theme_bw()+
  coord_cartesian(xlim = c(32740000,32755000))+
  xlim(32740000,32755000) +
  scale_y_continuous(limits=c(-415,735),labels=c("250","0","500"), breaks=c(-250,0,500),expand=c(0.01,0)) +
  theme(legend.position="none")

png(filename="Figure7.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
dd
dev.off()
```
## Code for SNP calling and statistics
First, use dDocent_ngs to call raw variants
```bash
wget https://raw.githubusercontent.com/jpuritz/EecSeq/master/Bioinformatics/SNPConfig
bash ./dDocent_ngs.sh SNPConfig
```
Next decompose raw variants into SNPs and InDel calls
```bash
vcfallelicprimitives TotalRawSNPs.vcf --keep-info --keep-geno > decomp.raw.vcf
```
Next, remove InDel calls
```bash
mawk '$0 ~ /#/ || $0 ~ /TYPE\=snp/' decomp.raw.vcf > snp.raw.vcf
```
Next, filter using VCFtools to generate stats

Total SNPs: `vcftools --vcf snp.raw.vcf `

Total SNPs with quality score higher than 20: `vcftools --vcf snp.raw.vcf --minQ 20`

Total Exome SNPS: `vcftools --vcf snp.raw.vcf --bed sorted.ref3.0.exon.sc.bed --minQ 20`

Exome SNPs with minimum mean 16X coverage: `vcftools --vcf snp.raw.vcf --minQ 20 --bed sorted.ref3.0.exon.sc.bed --min-meanDP 16`

Exome SNPs with minimum mean 32X coverage: `vcftools --vcf snp.raw.vcf --minQ 20  --bed sorted.ref3.0.exon.sc.bed --min-meanDP 32`

Exome SNPs with minimum mean 48X coverage: `vcftools --vcf snp.raw.vcf --minQ 20  --bed sorted.ref3.0.exon.sc.bed --min-meanDP 48`

Exome SNPs with minimum mean 80X coverage: `vcftools --vcf snp.raw.vcf --minQ 20  --bed sorted.ref3.0.exon.sc.bed --min-meanDP 80`

SNPs with minimum mean 80X coverage outside of exome: `vcftools --vcf snp.raw.vcf --minQ 20  --exclude-bed sorted.ref3.0.exon.sc.bed --min-meanDP 80`


## Code for analysis of hyRAD-X

Change to GFF directory

`cd $WORKING_DIR/Genome/GFF`

Get mRNA intervals

`mawk '$8 == "mRNA"' ref_C_virginica-3.0_top.bed > ref3.0_mRNA.bed`

Take only 1 transcript variant per mRNA

`bedtools merge -i ref3.0_mRNA.bed  -c 10 -o "first" | mawk '{ print $1"\t"$2"\t"$3"\tmRNA="NR}'  > ref3.0_mRNA.sc.bed`

Remove introns from bed file and then convert into fasta format
```bash
bedtools subtract -a ref3.0_mRNA.sc.bed -b cv.ref3.intron.sc.bed > ref3.0_mRNA_nointrons.bed
bedtools getfasta -name -tab -fi ~/j.puritz_remote/EecSeq2/reference.fasta -bed ref3.0_mRNA_nointrons.bed | sed 's/::.*\t/\t/g' | mawk 'BEGIN {m=0;seq=0} {if ($1 == m) seq= seq$2;else if ($1 != m) {print ">"m; print seq; m=$1; seq = $2}}' | mawk '!/>0/ && !/^0/' > mRNA.nointrons.fasta
```
Count number of mRNAs and the number of mRNAs with two restriction sites
```bash
mawk '/>/' mRNA.nointrons.fasta | wc -l
grep "TTAA.*TTAA" mRNA.nointrons.fasta | wc -l
```
Result:
```bash
31383
29555
```
### R code for SimRAD
```R
source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")
library(SimRAD)

RE1_1 <- "T"
RE1_2 <- "TAA"
RE2_1 <- "T"
RE2_2 <- "TAA"

transcripts <- ref.DNAseq("mRNA.nointrons.fasta", subselect.contigs = FALSE)
digest <-insilico.digest(transcripts,RE1_1,RE1_2,RE2_1,RE2_2, verbose=TRUE)
adapters <- adapt.select(digest, type= "AB",RE1_1,RE1_2,RE2_1,RE2_2)
wid.digest.adapt <- size.select(adapters, min.size = 100, max.size = 100000, graph=TRUE, verbose=TRUE)
```
Results:
```R
Number of type AA fragments:440881
220184 fragments between 100 and 1e+05 bp 
```
