#!/bin/bash

###Information
#reference filter workflow: https://github.com/gokul-r/RNAediting_RNAonly/blob/master/README_pipeline
#reference GATK code: http://people.duke.edu/~ccc14/duke-hts-2017/Statistics/08032017/GATK-pipeline-sample.html
#reference results: https://csibioinfo.nus.edu.sg/csingsportal/docs/#rnaediting


cd /exports/eddie/scratch/s1796194/foolab/sra

group="DCM" #either PPCM, DCM, HCM, Control
sample=`sed -n "$SGE_TASK_ID p" ../code/Failed_list_$group.txt`
ncore=4

module load java
module load python/2.7.10
module load igmm/apps/STAR/2.7.1a
#module load igmm/apps/picard/2.17.11
module load igmm/apps/GATK/4.0.10.1
module load igmm/apps/samtools/1.10
module load igmm/apps/bcftools/1.9
module load igmm/apps/BEDTools/2.27.1
module load igmm/apps/vcftools/0.1.13
module load igmm/apps/stringtie/1.3.5
module load igmm/apps/bowtie/2.3.5.1
module load igmm/apps/sratoolkit/2.10.8
export PATH=$PATH:/exports/cmvm/eddie/sbms/groups/young-lab/yiru/apps/bamUtil/bin
export PATH=$PATH:~/bin

# File directions
ref_dir="/exports/cmvm/eddie/sbms/groups/young-lab/yiru/refdata"
perl_dir="/exports/cmvm/eddie/sbms/groups/young-lab/yiru/perlCode"
apps_dir="/exports/cmvm/eddie/sbms/groups/young-lab/yiru/apps"
#output_dir="/exports/cmvm/eddie/sbms/groups/young-lab/yiru/output/$group"
final_output_dir="/exports/cmvm/eddie/sbms/groups/young-lab/yiru/final/$group"
annovar_dir=$apps_dir"/annovar"

# Program file
mypicard=$apps_dir"/picard.jar"
mygatk=$apps_dir"/gatk.jar"
myAEI=$apps_dir"/RNAEditingIndexer/RNAEditingIndex"

# Splicing junction database
ucsc_anno=$ref_dir"/hg19_UCSC_genes_20201219.gtf"
refseq_anno=$ref_dir"/hg19_RefSeq_curated_genes_20201219.gtf"
genecode_anno=$ref_dir"/hg19_GENCODE_basic_genes_V27lift37_20201219.gtf"
ensembl_anno=$ref_dir"/hg19_Ensembl_genes_20201219.gtf"
merged_anno=$ref_dir"/hg19.merged.gene.anno.txt"

# Annotations
refGenome=$ref_dir"/ucsc.hg19.fasta"
g1000Anno=$ref_dir"/1000G_omni2.5.hg19.sites.vcf"
exonAnno=$ref_dir"/hg19.GRCh38.wesp.snps_indels.chr.sorted.vcf"
dbSnpAnno=$ref_dir"/dbsnp_human_9606_b150_GRCh37p13.chr.vcf"
aluAnno=$ref_dir"/hg19_alu_anno.bed"
repeatMasker=$ref_dir"/hg19_repeatMasker.bed"
simpleRepeatsAnno=$ref_dir"/simpleRepeat.bed"
combinedSNP=$ref_dir"/combined.snp.chr.txt"
g1000Indel=$ref_dir"/1000G_phase1.indels.hg19.sites.vcf"
millsIndel=$ref_dir"/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"


######################## Sample analysis ##########################
time=$(date "+%Y-%m-%d %H:%M:%S")
echo $time $sample "data download start..\n"

prefetch $sample --output-directory /exports/eddie/scratch/s1796194/foolab/sra
cd /exports/eddie/scratch/s1796194/foolab/sra/$sample
fastq-dump --split-files --gzip $sample.sra
#cp -r /exports/cmvm/eddie/sbms/groups/young-lab/yiru/refdata ./
#cp -r /exports/cmvm/eddie/sbms/groups/young-lab/yiru/perlCode ./
rm $sample.sra

sampleR1=$sample"_1.fastq.gz"
sampleR2=$sample"_2.fastq.gz"

time=$(date "+%Y-%m-%d %H:%M:%S")
echo $time $sample "data download done\n"


##################### Mapping of RNAseq reads ####################

time=$(date "+%Y-%m-%d %H:%M:%S")
echo $time "STAR alignment for $sample started...\n"

# Step 1. first pass alignment by STAR
## star 1-pass index
mkdir star_idx1
STAR --runThreadN $ncore --runMode genomeGenerate \
	--genomeDir ./star_idx1 \
	--genomeFastaFiles $refGenome \
	--sjdbGTFfile $ucsc_anno $refseq_anno $genecode_anno $ensembl_anno
## star 1-pass align
mkdir star_1pass
STAR --runThreadN $ncore --genomeDir ./star_idx1 \
	--readFilesIn $sampleR1 $sampleR2 \
	--readFilesCommand zcat \
	--outFileNamePrefix ./star_1pass/$sample.
rm -r star_idx1

# Step 2. second pass alignment by STAR with recognized splicing junctions
## star 2-pass index
mkdir star_idx2
STAR --runThreadN $ncore --runMode genomeGenerate \
        --genomeDir ./star_idx2 \
        --genomeFastaFiles $refGenome \
        --sjdbFileChrStartEnd ./star_1pass/$sample.SJ.out.tab
rm -r star_1pass
## star 2-pass align
mkdir star_2pass
STAR --runThreadN $ncore --genomeDir ./star_idx2 \
        --readFilesIn $sampleR1 $sampleR2 \
        --readFilesCommand zcat \
        --outFileNamePrefix ./star_2pass/$sample.
rm -r star_idx2
rm $sampleR1 $sampleR2

# Step 3. sort & add read groups
java -jar $mypicard AddOrReplaceReadGroups \
	I=./star_2pass/$sample.Aligned.out.sam \
	O=$sample.star.sort.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true \
	RGID=foo RGSM=$sample RGLB=$sample RGPL=Illumina RGPU=$sample
rm -r star_2pass

#Step 4. use Picard MarkDuplicates to remove duplicate reads
java -jar $mypicard MarkDuplicates \
	I=$sample.star.sort.bam \
	O=$sample.star.sort.rmdup.bam \
	M=$sample.star.sort.rmdup.metric \
	ASSUME_SORT_ORDER=coordinate \
	REMOVE_DUPLICATES= false \
	VALIDATION_STRINGENCY=LENIENT \
	CREATE_INDEX=true

#Step 5. delete intron reads and MAPQ change for GATK
java -jar $mygatk -T SplitNCigarReads \
	-R $refGenome \
	-I $sample.star.sort.rmdup.bam \
	-o $sample.star.sort.rmdup.split.bam \
	-rf ReassignOneMappingQuality \
	-RMQF 255 \
	-RMQT 60 \
	-U ALLOW_N_CIGAR_READS
rm $sample.star.sort.rmdup.ba* $sample.star.sort.rmdup.metric

#Step 6. filter out unmapped reads and reads with mapping quality < 20 with samtools
samtools view -h -b -q20 $sample.star.sort.rmdup.split.bam > $sample.star.sort.rmdup.split.q20.bam
rm $sample.star.sort.rmdup.split.ba*

#Step 7. Index this alignment file with samtools
samtools index $sample.star.sort.rmdup.split.q20.bam

#Step 8. Perform indel realignment in GATK
java -jar $mygatk -R $refGenome \
	-T RealignerTargetCreator \
	-I $sample.star.sort.rmdup.split.q20.bam \
	-o $sample.realign.intervals \
	-known $dbSnpAnno \
	-known $exonAnno \
	-known $g1000Anno
java -jar $mygatk -T IndelRealigner \
	-R $refGenome \
	-targetIntervals $sample.realign.intervals \
	-I $sample.star.sort.rmdup.split.q20.bam \
	-o $sample.star.sort.rmdup.split.q20.realign.bam
java -jar $mypicard BuildBamIndex \
	I=$sample.star.sort.rmdup.split.q20.realign.bam
rm $sample.realign.intervals $sample.star.sort.rmdup.split.q20.bam $sample.star.sort.rmdup.split.q20.bam.bai

#Step 9. Perform base quality score recalibration in GATK
mkdir log
java -jar $mygatk -T BaseRecalibrator \
	-R $refGenome \
	-I $sample.star.sort.rmdup.split.q20.realign.bam \
	-o $sample.recal1.grp \
	-knownSites $g1000Indel \
	-knownSites $millsIndel \
	-knownSites $dbSnpAnno \
	-knownSites $exonAnno  \
	-knownSites $g1000Anno \
	-log ./log/$sample.recal1.grp.log
java -jar $mygatk -T PrintReads -R $refGenome \
	-I $sample.star.sort.rmdup.split.q20.realign.bam \
	-o $sample.star.sort.rmdup.split.q20.realign.recal.bam \
	-BQSR $sample.recal1.grp \
	-log ./log/$sample.star.sort.rmdup.realign.recal.bam.log
java -jar $mypicard BuildBamIndex \
	I=$sample.star.sort.rmdup.split.q20.realign.recal.bam
rm $sample.star.sort.rmdup.split.q20.realign.ba* $sample.recal1.grp

time=$(date "+%Y-%m-%d %H:%M:%S")
echo $time "STAR alignment is done."

##################### Variant calling and filtering ####################

time=$(date "+%Y-%m-%d %H:%M:%S")
echo $time "SNP calling started..."

#Step 1. Initial variant calling with GATK UnifiedGenotyper
samtools mpileup -DSugf $refGenome $sample.star.sort.rmdup.split.q20.realign.recal.bam -q 20 -Q 25 | bcftools call -mv > $sample.samtools.raw.vcf
#min mismatch ??
#bcftools view -Ncvg - ??

#Step 2.1. Convert VCF format to our custom variant format
perl $perl_dir/Convert_VCF.pl $sample.samtools.raw.vcf $sample.samtools.raw.txt
rm $sample.samtools.raw.vcf

#Step 2.2. Filter out candidates that overlap with dbSNP/1000 genomes/UW exome calls
bedtools intersect -a $sample.samtools.raw.txt -b $combinedSNP -v > $sample.samtools.unique.txt
rm $sample.samtools.raw.txt

#Step 3. Remove mismatches in first 6 bp of reads
perl $perl_dir/Remove_mismatch_first6bp.pl $sample.samtools.unique.AGCT.txt $sample.star.sort.rmdup.split.q20.realign.recal.bam $sample.samtools.unique.AGCT.rm6bp.txt
rm $sample.samtools.unique.AGCT.txt

#Step 4. Separate candidates in Alu and non-Alu regions of the genome
bedtools intersect -a $sample.samtools.unique.AGCT.rm6bp.txt -b $aluAnno -wa > $sample.samtools.unique.AGCT.rm6bp.Alu.txt
bedtools intersect -a $sample.samtools.unique.AGCT.rm6bp.txt -b $aluAnno -v > $sample.samtools.unique.AGCT.rm6bp.nonAlu.txt
rm $sample.samtools.unique.AGCT.rm6bp.txt

time=$(date "+%Y-%m-%d %H:%M:%S")
echo $time "SNP calling is done."


##################### Further filtering of nonAlu ####################
time=$(date "+%Y-%m-%d %H:%M:%S")
echo $time "nonAlu filtering started..."

#Step 1. Use bedtools to filter candidates that are in simple repeats
bedtools intersect -a $sample.samtools.unique.AGCT.rm6bp.nonAlu.txt -b $simpleRepeatsAnno -v > $sample.nonAlu.rmSR.txt

#Step 2. Filter intronic candidates that are within 4 bp of splicing junctions
perl $perl_dir/Filter_intron_near_splicejuncts.pl $sample.nonAlu.rmSR.txt $merged_anno > $sample.nonAlu.rmSR.rm4bp.txt

#Step 3. Filter candidates in homopolyer runs
perl $perl_dir/RemoveHomoNucleotides.pl $sample.nonAlu.rmSR.rm4bp.txt $refGenome $sample.nonAlu.rmSR.rm4bp.rmhp.txt

##Optional: remove chrM to accelerate speed
#awk '$1=="chrM"{print}' $sample.nonAlu.rmSR.rm4bp.rmhp.txt > $sample.nonAlu.rmSR.rm4bp.rmhp.chrM.txt
#awk '$1!="chrM"{print}' $sample.nonAlu.rmSR.rm4bp.rmhp.txt > $sample.nonAlu.rmSR.rm4bp.rmhp.noM.txt

#Step 4. Use BLAT to ensure unique mapping
perl $perl_dir/BLAT_candidates.pl $ncore $sample.nonAlu.rmSR.rm4bp.rmhp.txt $sample.star.sort.rmdup.split.q20.realign.recal.bam $refGenome $sample.nonAlu.rmSR.rm4bp.rmhp.blat.txt

#Step 5. Use bedtools to separate out candidates in repetitive non-Alu and nonrepetitive regions of the genome
bedtools intersect -a $sample.nonAlu.rmSR.rm4bp.rmhp.blat.txt -b $repeatMasker -wa > $sample.nonAlu.repetitive.txt
bedtools intersect -a $sample.nonAlu.rmSR.rm4bp.rmhp.blat.txt -b $repeatMasker -v > $sample.nonAlu.unique.txt

#Step 6. change to standard format and annotate
awk '{print $1,$2,$3,$5,$6,"nonAlu",$7,$4}' $sample.nonAlu.repetitive.txt | sed 's/ /\t/g;s/,/\t/g' > $sample.candidate.nonAlu.txt
awk '{print $1,$2,$3,$5,$6,"Alu",$7,$4}' $sample.samtools.unique.rm6bp.Alu.txt | sed 's/ /\t/g;s/,/\t/g' > $sample.candidate.Alu.txt
awk '{print $1,$2,$3,$5,$6,"Uniq",$7,$4}' $sample.nonAlu.unique.txt | sed 's/ /\t/g;s/,/\t/g' > $sample.candidate.Uniq.txt
cat $sample.candidate.*.txt > $sample.candidate.all.txt

time=$(date "+%Y-%m-%d %H:%M:%S")
echo $time "nonAlu filtering is done."


##################### ANNOVAR annotation ####################
time=$(date "+%Y-%m-%d %H:%M:%S")
echo $time "ANNOVAR annotation started..."

perl $annovar_dir/table_annovar.pl $sample.candidate.all.txt $annovar_dir/humandb/ -buildver hg19 -out $sample.all -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -polish -xref $annovar_dir/example/gene_xref.txt

cut -f '1-10' $sample.all.hg19_multianno.txt > anno.tmp
echo "#"`sed -n '1p' anno.tmp` "editing_site Mut_Freq Coverage Mut_read_count" | sed 's/ /\t/g' > colname.tmp
sed -i '1d' anno.tmp
cut -f '6-9' $sample.candidate.all.txt > call.tmp
paste anno.tmp call.tmp | sed 's/ /_/g'> $sample.all.tmp
awk -F '\\.\t' '{print $0"\t"NF-1}'  $sample.all.tmp | awk '$15!=5 {$15=null; print}' | sed 's/ /\t/g'> $sample.annotated.tmp
bedtools sort -i $sample.annotated.tmp > $sample.annotated.sorted.tmp
cat colname.tmp $sample.annotated.sorted.tmp > $sample.editing_candidates_annotated.txt
rm *.tmp

awk -F "\t" '$4=="A"{print}' $sample.editing_candidates_annotated.txt > $sample.editing_candidates_annotated.AG.txt
awk -F "\t" '$5=="C"{print}' $sample.editing_candidates_annotated.txt > $sample.editing_candidates_annotated.CT.txt


time=$(date "+%Y-%m-%d %H:%M:%S")
echo $time "ANNOVAR annotation is done."


##################### Save Files ####################
cp $sample.candidate.nonAlu.txt $sample.candidate.Alu.txt  $sample.candidate.Uniq.txt $sample.editing_candidates_annotated.AG.txt $$sample.editing_candidates_annotated.CT.txt $final_output_dir

cd ../
rm -r /exports/eddie/scratch/s1796194/foolab/sra/$sample