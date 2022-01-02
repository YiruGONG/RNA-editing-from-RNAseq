#!/bin/bash

###################### software download #########################
#picard
wget https://github.com/broadinstitute/picard/releases/download/2.24.1/picard.jar
#GATK
#https://storage.cloud.google.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

#BLAT
wget https://genome-test.gi.ucsc.edu/~kent/src/blatSrc36.zip
unzip blatSrc36.zip
cd blatSrc
uname -a
export MACHTYPE=x86_64
mkdir -p ~/bin/x86_64
make
cd ~
vim .bashrc      (export PATH=/home/s1796194/bin/x86_64:$PATH)
source .bashrc

#ANNOVAR
cd /exports/cmvm/eddie/sbms/groups/young-lab/yiru/
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
tar xvfz annovar.latest.tar.gz
cd annovar
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a humandb/

#RNAEditingIndexer (AEI calculation)
git clone https://github.com/a2iEditing/RNAEditingIndexer.git
git clone https://github.com/statgen/bamUtil.git   ##for banUtils required in this software


######################## data download ##########################

####### sample data #########

wget https://www.encodeproject.org/files/ENCFF855TSM/@@download/ENCFF855TSM.fastq.gz
mv ENCFF855TSM.fastq.gz ENCFF855TSM_1.fastq.gz
wget https://www.encodeproject.org/files/ENCFF971MYL/@@download/ENCFF971MYL.fastq.gz
mv ENCFF971MYL.fastq.gz ENCFF971MYL_2.fastq.gz

####### result reference #########
#CSI portal
cd /exports/cmvm/eddie/sbms/groups/young-lab/yiru/RNA-Editing


####### reference data #########

##1) GATK annotations: hg19 reference seq, 1000 genome variant calling sites (SNP) and GATK known indels
lftp ftp.broadinstitute.org -u gsapubftp-anonymous
cd /bundle/hg19
mget 1000G_omni2.5.hg19.sites.vcf* 1000G_phase1.indels.hg19.sites.vcf* ucsc.hg19* Mills_and_1000G_gold_standard.indels.hg19.sites.vcf*
quit
for i in $(ls 1000G_omni2.5.*.gz 1000G_phase1.*.gz Mills_*.gz ucsc.hg19*.gz);do gunzip $i;done

##make refseq index
samtools faidx ucsc.hg19.fasta
#bowtie2-build ucsc.hg19.fasta ucsc.hg19.fasta
#bwa index ucsc.hg19.fasta

##2) dbsnp from UCSC (all/common/clinVar/Mult/BadCoords)
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/archive/common_all_20170403.vcf.gz
#wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/archive/common_all_20170403.vcf.gz.md5
#wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/archive/common_all_20170403.vcf.gz.tbi
gunzip common_all_20170403.vcf.gz
mv common_all_20170403.vcf dbsnp_human_9606_b150_GRCh37p13.vcf
awk '{ 
        if($0 !~ /^#/) 
            print "chr"$0;
        else if(match($0,/(##contig=<ID=)(.*)/,m))
            print m[1]"chr"m[2];
        else print $0 
      }' dbsnp_human_9606_b150_GRCh37p13.vcf > dbsnp_human_9606_b150_GRCh37p13.chr.vcf

##3) WESP known indel:
wget http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz
mv ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz GRCh38.wesp.snps_indels.vcf.tar.gz
tar xvf GRCh38.wesp.snps_indels.vcf.tar.gz
mkdir GRCh38.wesp.snps_indels
mv ESP6500SI-V2-SSA137.GRCh38-liftover.* ./GRCh38.wesp.snps_indels
cd GRCh38.wesp.snps_indels
vcf-concat *.vcf | vcf-sort > hg19.GRCh38.wesp.snps_indels.vcf
cp hg19.GRCh38.wesp.snps_indels.vcf ../
cd ../
awk '{ 
        if($0 !~ /^#/) 
            print "chr"$0;
        else if(match($0,/(##contig=<ID=)(.*)/,m))
            print m[1]"chr"m[2];
        else print $0 
      }' hg19.GRCh38.wesp.snps_indels.vcf > hg19.GRCh38.wesp.snps_indels.chr.vcf
java -jar $mypicard SortVcf \
    I=hg19.GRCh38.wesp.snps_indels.chr.vcf \
    O=hg19.GRCh38.wesp.snps_indels.chr.sorted.vcf \
    SEQUENCE_DICTIONARY=ucsc.hg19.dict

##4) combined known SNP data from dbSNP/1000genome/UW ExonSeq
sed 's/##contig=<ID=chr/##contig=<ID=/g' 1000G_omni2.5.hg19.sites.vcf | sed 's/^chr//g'  > 1000G_omni2.5.hg19.sites.nochr.vcf
vcf-concat dbsnp_human_9606_b150_GRCh37p13.vcf 1000G_omni2.5.hg19.sites.nochr.vcf hg19.GRCh38.wesp.snps_indels.vcf | vcf-sort > combined.snp.vcf

###add "chr"
awk '{ 
        if($0 !~ /^#/) 
            print "chr"$0;
        else if(match($0,/(##contig=<ID=)(.*)/,m))
            print m[1]"chr"m[2];
        else print $0 
      }' combined.snp.vcf > combined.snp.chr.vcf
##transfer vcf file to bed-like file
perl ./perlCode/Convert_ref_VCF.pl combined.snp.chr.vcf combined.snp.chr.txt

#5) other annotation data download from UCSC browser: 
#hg19_alu_anno.txt
#simpleRepeat.bed
#hg19_repeatMasker.bed
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.ensGene.gtf.gz | gunzip
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.ncbiRefSeq.gtf.gz | gunzip

##splicing Junction files (known genes)
ucsc_anno="hg19_UCSC_genes_20201219.gtf"
refseq_anno="hg19_RefSeq_curated_genes_20201219.gtf"
genecode_anno="hg19_GENCODE_basic_genes_V27lift37_20201219.gtf"
ensembl_anno="hg19_Ensembl_genes_20201219.gtf"

#merged known genes
stringtie --merge $ucsc_anno $refseq_anno $genecode_anno $ensembl_anno -m 0 -o hg19.merged.gene.anno.txt

#splicing junction sequence
wget https://github.com/gokul-r/RNAediting_RNAonly/releases/download/1/hg19_junctionseqs_75bp.fa
cat ucsc.hg19.fasta hg19_junctionseqs_75bp.fa > hg19.withJunctions75bp.fa
samtools faidx hg19.withJunctions75bp.fa
bwa index hg19.withJunctions75bp.fa
