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

```bashbedtools subtract -a sorted.ref3.0.exon.bed -b sorted.ref3.0.CDS.bed -g genome.file -sorted > sorted.ref3.0.UTR.bed
bedtools merge -i <(bedtools sort -i sorted.ref3.0.UTR.bed -faidx <(cut -f1 genome.file))> sorted.ref3.0.UTR.sc.b
bedtools sort -i sorted.ref3.0.UTR.sc.b -faidx <(cut -f1 genome.file) > sorted.ref3.0.UTR.sc.bed
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
```bashsamtools view -@64 -q4 -h -F 0x100 -F 0x400 m4.merged.bam| mawk '$6 !~/[8-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | samtools view -b > m4.q4.merged.bam
```
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

