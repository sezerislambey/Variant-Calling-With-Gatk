module load java
module load picard_tools
module load gatk
module load samtools
module load bcftools
module load bwa

# #download reference genome
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human*
# gunzip human_g1k_v37.fasta.gz
# java -jar /home/u/anaconda3/share/picard-2.27.4-0/picard.jar CreateSequenceDictionary R=human_g1k_v37.fasta O=human_g1k_v37.dict
# 
# #download modern human known-sites (hg19)
# 1kg snps and indel
# wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.*
# wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps*
# 
# #other database
# wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5*
# wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf*
# wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3.hg19.sites.vcf*
# 
# #download practice bam files (chr20, GBR, 1000 genomes)
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00096/alignment/HG00096.chrom20*
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00100/alignment/HG00100.chrom20*
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00108/alignment/HG00108.chrom20*
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00111/alignment/HG00111.chrom20*

# 

# vcf file view example:
zless Mills_and_1000G_gold_standard.indels.b37.vcf.gz


# Mapping
bwa-mem2 mem -t 1 -R "@RG\tID:sample\tSM:sample\tPL:platform" ref.fa /sample_1.fastq.gz sample_2.fastq.gz > sample.mapped.sam


# sam > bam
samtools view -Sb sample.mapped.sam > sample.mapped.bam
samtools view -h sample.mapped.bam | less -S


# Basic syntax
samtools view HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20101123.bam | head -n 10
samtools view HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20101123.bam 20:60000-60010
samtools index HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20101123.bam
gatk --java-options -Xmx4G BaseRecalibrator


# mark duplicates
java -jar /home/u/anaconda3/share/picard-2.27.4-0/picard.jar MarkDuplicates INPUT=HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20101123.bam OUTPUT=dedup_HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20101123.bam M=HG00096.chrom20_metrics.txt


# add read group
java -jar /home/u/anaconda3/share/picard-2.27.4-0/picard.jar AddOrReplaceReadGroups I=dedup_HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20101123.bam O=RG_dedup_HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20101123.bam RGID=GBR RGLB=lib1 RGPL=ILLUMINA SORT_ORDER=coordinate RGPU=unit1 RGSM=HG00096


# 3. base quality score recalibration (BQSR)
gatk --java-options -Xmx4G BaseRecalibrator -R ../ref_genome/human_g1k_v37.fasta -I RG_dedup_test.bam --known-sites ../known_sites/1000G_phase1.indels.hg19.sites.clean.vcf.gz --known-sites ../known_sites/1000G_phase1.snps.high_confidence.hg19.sites.clean.vcf.gz -O test_BQSR.txt


# 4. apply BQSR
gatk --java-options -Xmx4G ApplyBQSR -R ../ref_genome/human_g1k_v37.fasta -I RG_dedup_test.bam --bqsr-recal-file test_BQSR.txt -O BQSR_RG_dedup_test.bam
# or
java -jar apps/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar ApplyBQSR -R ../ref_genome/human_g1k_v37.fasta -I RG_dedup_test.bam --bqsr-recal-file test_BQSR.txt -O BQSR_RG_dedup_test.bam


# 5. call genotypes and get gvcfs #repeat this step for each sample
gatk --java-options -Xmx4G HaplotypeCaller -R ../ref_genome/human_g1k_v37.fasta -I BQSR_RG_dedup_test.bam -O ../vcf/BQSR_RG_dedup_test.g.vcf.gz -bamout BQSR_RG_dedup_test.bamout.bam -ploidy 2 -ERC GVCF -L 20


# 6. consolidate gvcfs
gatk --java-options -Xmx4G GenomicsDBImport -R ../ref_genome/human_g1k_v37.fasta -V BQSR_RG_dedup_test.g.vcf.gz -V BQSR_RG_dedup_test2.g.vcf.gz -V BQSR_RG_dedup_test3.g.vcf.gz --genomicsdb-workspace-path chr20_test --intervals 20


# 7. jointly call variants
gatk --java-options -Xmx4G GenotypeGVCFs -R ../ref_genome/human_g1k_v37.fasta -V gendb://chr20_test -O test.vcf


# 8. separate snps and indels
gatk --java-options -Xmx4G SelectVariants -V test.vcf -select-type SNP -O test.snp.vcf
gatk --java-options -Xmx4G SelectVariants -V test.vcf -select-type INDEL -O test.indel.vcf


# 9a. hard filtering
# snp
gatk --java-options -Xmx4G VariantFiltration -V test.snp.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O {}
# indel
gatk --java-options -Xmx4G VariantFiltration -V test.indel.vcf -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O {}

# 9b.1. VQSR
gatk --java-options -Xmx4G VariantRecalibrator -V test.snp.vcf --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP -mode SNP --max-gaussians 6 --resource:1000G,known=false,training=true,truth=true,prior=10.0 ../known_sites/1000G_phase1.snps.high_confidence.hg19.sites.clean.vcf.gz -O test.snp.recal.txt --tranches-file test.snp.tranch.txt --rscript-file test.snp.R.plot

# 9b.2 apply VQSR
gatk -java-options -Xmx4G ApplyVQSR -V test.snp.vcf -O test.snp.VQSR.vcf --truth-sensitivity-filter-level 99.7 --tranches-file test.snp.tranch.txt --recal-file test.snp.recal.txt -mode SNP --create-output-variant-index true












