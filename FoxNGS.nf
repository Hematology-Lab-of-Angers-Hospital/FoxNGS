fastqSetup
bcl + sasmplesheet
bcl2fastq -> fastq

process qualityControl {
    input:
    path fastq from fastqChannel

    output:
    publishDir "${params.outDir}/QC", mode: 'copy'

    script:
    """
    fastqc 
    """

}

process buildIndex {
    publishDir "${params.publish_dir}/"

    input:
    path

    output:
    path sam_filt into samChannel

    script:
    """
    bwa mem -t params.bwa.index.threads -R params.bwa.index.header $reference_genome $fastq_r1 $fastq_r2 > sam_file
    """
}

proccess bamSetup {

}

samtools view
samtools sort
samtools index
samtools view -F ox4
samtools index
samtools coverage
samtools intersect
samtools flagstat
picard MarkDuplicates
samtools index -> sort.dupmark.bam

process qCReport {
    input:
    path bam from bamChannel

    output:
    path qc_report into ???
}

bedtools coverage > cover_analysis
R -e "rmarkdown::render('${RSCRIPT}', params = list(directory='$(pwd)',file='${name}_couverture_analyse.bed',agent='${agent}',pipeline='${0}',output='${COUV}/Statistic_couverture.csv',output_gene='/media/$USER/Elements/Result_NGS/Stat_gene/Statistic_couverture_gene.csv',ratio_library='${int_ratio}'),output_file='$(pwd)/${name}_couverture_analyse.bed.html')" -> _couverture_analysis.bed.html

/*
-> sort.dupmark.bam
samtools pileup
java varscan mpileup2cns -> .varsan.vcf
*/
process varscanVariantCaller {
    input:
    path bam from bamChannel

    output:
    path varscanVariants into variantChannel

    script:
    """
    samtools mpileup -Q 13 -q 0 -A -B -d 100000 -f $refFasta $bam > mpileup
    java -jar src/varscan mpileup2cns mpileup --min-coverage 50 --min-reads2 8 --min-avg-qual 30 --min-var-freq 0.02 --p-value 0.1 --strand-filter 0 --output-vcf --variants > varscanVariants
    """
}

process gatkVariantCaller {
    input:
    path bam from bamChannel

    output:
    path gatk_variants into variantChannel

    script:
    """
    java -jar src/GATK BaseRecalibrator -I $bam -r $refFasta --known-sites -o recal_data_table
    java -jar src/GATK ApplyBQRS -I recal_data_table -r $refFasta -o bqsr_bam
	java -jar src/GATK HaplotypeCaller -I bqsr_bam -R $refFasta -O gatk_output --min-base-quality-score 30 --minimum-mapping-quality 20 --dont-use-soft-clipped-bases true
    awk '{gsub(",Date=[^>]+>",">");}1' gatk_output
    """
}

/*
-> sort.dupmark.bam + RefGen
GATK Mutec2
gunzip ->  mutec2.vcf
*/
process gatkVariantCaller {
    input:
    file bam from bamChannel

    output:
    file mutec2Variants into variantChannel

    script:
    """
    java -jar src/GATK Mutect2 -I bam -R $refFasta --min-base-quality-score 30 --dont-use-soft-clipped-bases true -O zippedMutec2Variant
    gunzip zippedMutec2Variants > mutec2Variants"
    """
}

/*
-> sort.dupmark.bam + RefGen + config_file_pindel.txt
pindel
pindel -> pindel.vcf
*/
process pindelVariantCaller {
    input:
    path bam from bamChannel

    output:
    path pindelVariants into variantChannel

    script:
    """
    src/pindel -f $BWA_FASTA -i $bam -j $BED_PINDEL -o ITD
    src/pindel -P ITD -r $BWA_FASTA  -R x -d 00000000 -G -v pindelVariants
    """-y
}

.gatk.vcf, .mutec2.vcf, .varscan.vcf, .pindel.vcf
python3 annotations_processing -> fusion_annotation_simplified_$name_.csv