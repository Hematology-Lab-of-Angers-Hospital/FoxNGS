/// ----- F O X   N G S   M O D U L E S -----



// ___PROCESSES___


process FASTQ_SETUP {
    input:
        path(bcl)
        path(singularity_image)
        path(samplesheet)

    output:
        path()

    script:
    """
    singularity exec ${singularity_image} bcl2fastq --runfolder ${bcl_runfolder} --output-dir $REPERTORY/ --sample-sheet ${samplesheet} --use-bases-mask Y151,I8,I8,Y151 --no-lane-splitting -r $task.cpus -p $task.cpus -w $task.cpus
    """
}


process QUALITY_CONTROL {
/*
Runs FastQC on fastq input to check read quality

WARNING: This quality control step will NOT stop the analysis
    even if the reads are too low quality for significant results.
    
    Check quality output in the results/$patient_id/QC directory before any
    other analysis result

INPUT
    < fastq_r1 and fastq_r2: NGS reads 
        (both files if doing paired-end, only fastq_r1 if doing single-end)

OUTPUT
    > QC report file in the 
        (.html files in the 'QC' folder of the patient directory)

INPUT FROM
    <- setup_paired_end_reads channel

OUTPUT INTO
    -> results/$patient_id/QC directory
*/
    cpus 2
    publishDir "$params.results/$sampleId/QC", mode: 'copy'

    input:
        tuple val(sampleId), path(fastq_r1), path(fastq_r2)

    output:
        path '*'

    script:
    """
    fastqc \
        -t ${task.cpus} \
        ${fastq_r1} \ 
        ${fastq_r2}
    """
}

process FAI_SETUP {
/*
Runs samtools faidx to make the reference genome index (.fai file)
if that is necessary

TIP: This process is not as long as INDEX_SETUP but the pipeline has parameters
to fetch the .fai file from the data/reference folder; consider copying the
.fai output file and set param.indexed_genome to 'true' to skip this step.

INPUT
    < reference_genome: path to the reference genome (.fna file)

OUTPUT
    > .fai genome index file

INPUT FROM 
    <- reference_genome channel

OUTPUT INTO
    -> indexed_genome channel
*/
    input:
        path(reference_genome)

    output:
        path("${reference_genome}.fai"), emit: indexed_genome

    script:
    """
    samtools faidx ${reference_genome} > ${reference_genome}.fai
    """
}

process REFERENCE_DICT_SETUP {
/*
Runs samtools faidx to make the reference genome index (.fai file)
if that is necessary

TIP: This process is not as long as INDEX_SETUP but the pipeline has parameters
to fetch the .fai file from the data/reference folder; consider copying the
.fai output file and set param.indexed_genome to 'true' to skip this step.

INPUT
    < reference_genome: path to the reference genome (.fna file)
    < gatk            : path to the GATK toolkit java archive (.jar) file

OUTPUT
    > .dict genome dictionnnary file

INPUT FROM 
    <- reference_genome channel
    <- gatk channel

OUTPUT INTO
    -> reference_dict channel
*/
    input:
        path(gatk)
        path(reference_genome)
        
    output:
        path("*.dict")

    script:
    """
    java -jar ${gatk} CreateSequenceDictionary -R ${reference_genome}
    """
}

process INDEX_SETUP {
/*
Runs bwa index on the reference genome

TIP: This is a long process. Consider fetching the output files in the work
directory, moving them to data/reference and setting the makeIndex parameter
to 'false' to skip this process 

INPUT
    < reference_genome: path to the reference genome (.fna file)

OUTPUT
    (linked)
    > amb, ann, bwt, pac and sa genome index files 

INPUT FROM
    <- reference_genome channel

OUTPUT INTO
    -> 
*/

    input:
        path(reference_genome)

    output:
        path("*")

    script:
    """
    bwa index ${reference_genome}
    """
}

process SAM_SETUP {
    cpus 4
/*
Sets up the sam file by calling bwa mem to align the fastq on 
the reference genome.
If the path to index files has not been passed to parameters or the index files
are absent, the pipeline builds the index files with the INDEX_SETUP beforehand

INPUTS
(linked)
    < sampleId             : patient ID value
    < fastq_r1 and fastq_r2: NGS reads 
        (both files if doing paired-end, only fastq_r1 if doing single-end)
    
    < reference_genome     : path to the reference genome (.fna file)
    
(linked)
    < .amb, .ann, .bwt, .pac and .sa genome index files 

OUTPUT
(linked)
    < sampleId: patient ID value
    < sam     : path to the sam (Sequence Alignment Map, aligned reads) file 

INPUT FROM
    <- reference_genome          : reference_genome channel

    <- amb, ann, bwt, pac, sa:
    reference_index channel 
        if the reference genome has already been indexed)
    OR
    INDEX_SETUP process
        (if the genome index files are not found)

    <- sampleId,fastq_r!,fastq_r2: setup_reads channel

OUTPUT INTO
    -> BAM_SETUP process
*/
    input:
        tuple val(sampleId), path(fastq_r1), path(fastq_r2)
        path(reference_genome)
        tuple path(amb), path(ann), path(bwt), path(pac), path(sa)
    
    output:
        tuple val(sampleId), path("${sampleId}.sam")

    script:
    """
    bwa mem -t $task.cpus -R "@RG\\tID:C5-${sampleId}\\tPL:illumina\\tPU:HXXX\\tLB:Solexa\\tSM:C5-${sampleId}" ${reference_genome} ${fastq_r1} ${fastq_r2} -o ${sampleId}.sam
    """
}

process BAM_SETUP {
/*
Sets up the bam (binary ve

rsion of the reads mapped to the genome) 
by calling samtools view and samtools sort

INPUTS
(linked)
    < sampleId: patient ID value
    < sam     : path to the sam (Sequence Alignment Map, aligned reads) file 
OUTPUT
(linked)
    > sampleId  : patient ID value
    > sorted_bam: path to the sorted bam file

INPUT FROM
    <- SAM_SETUP process

OUTPUT INTO
    -> BAM_SETUP process (sorted_bam)
    -> STDOUT (mapping stats)
*/
    input:
        tuple val(sampleId), path("${sampleId}.sam")

    output:
        tuple val(sampleId), path("${sampleId}_sorted.bam"), emit: sorted_bam
        tuple env(mapped_reads), env(unmapped_reads), env(chimeric_reads), emit: mapping_stats
        
    script:
    """
    samtools view -@ $task.cpus -Sh ${sampleId}.sam -bo ${sampleId}.bam
    samtools sort -@ $task.cpus ${sampleId}.bam -o ${sampleId}_sorted.bam

    mapped_reads=`samtools view -h -c ${sampleId}_sorted.bam`
    unmapped_reads=`samtools view -f 0x4 -h -@ $task.cpus -c -b ${sampleId}_sorted.bam`
    chimeric_reads=`samtools view  -f 0x800 -h -@ $task.cpus -c -b ${sampleId}_sorted.bam`
    """
}

process BAM_MAPPING {
/*
Filters the reads not mapped to the reference genome (ox4 tag) with samtools 
Sets up the bam index (bai) with samtools index

INPUTS
(linked)
    < sampleId  : patient ID value
    < sorted.bam: path to the sorted bam file

OUTPUT
(linked)
    > sampleId             : patient ID value
    > mapped_sorted.bam    : path to file containing 
                         the reads correctly mapped on the genome 
    > mapped_sorted.bam.bai: path to the mapped_sorted_bam index

INPUT FROM
    <- BAM_SETUP process

OUTPUT INTO
    -> DUPMARK_BAM_SETUP, COVERAGE_ANALYSIS processes
*/
    input:
        tuple val(sampleId), path("${sampleId}_sorted.bam")

    output:
        tuple val(sampleId), path("${sampleId}_mapped_sorted.bam"), path("${sampleId}_mapped_sorted.bam.bai")

    script:
    """
    samtools view -F ox4 -h -@ $task.cpus -b ${sampleId}_sorted.bam > ${sampleId}_mapped_sorted.bam
    samtools index -@ $task.cpus ${sampleId}_mapped_sorted.bam > ${sampleId}_mapped_sorted.bam.bai
    """
}

process ON_TARGET_MAPPING {
/*


INPUTS
(linked)
    < sampleId             : patient ID value
    < mapped_sorted_bam    : path to file containing 
                             the reads correctly mapped on the genome 
    < mapped_sorted_bam.bai: path to the mapped_sorted_bam index

    < bed_hemato       : coordinates of targets relevant to the analysis

OUTPUT
(linked)
    > sampleId     : patient ID value
    > on_target_bam: path to binary mapped with reads that are sorted and 
                     off target reads filtered
    > on_target_bam.bai: path to the on_target bam index

INPUT FROM
    <- BAM_MAPPING process

OUTPUT INTO
    -> DUPMARK_BAM_SETUP, COVERAGE_ANALYSIS processes
*/
    input:
        tuple val(sampleId), path("${sampleId}_mapped_sorted.bam"), path("${sampleId}_mapped_sorted.bam.bai")
        path(bed_hemato)

    output:
        tuple val(sampleId), path("${sampleId}_on_target.bam"), path("${sampleId}_on_target.bam.bai")

    script:
    """
    bedtools intersect -a ${sampleId}_mapped_sorted.bam -b ${bed_hemato} > ${sampleId}_on_target.bam
    samtools index -@ $task.cpus ${sampleId}_on_target.bam > ${sampleId}_on_target.bam.bai
    """
}

process DUPMARK_BAM_SETUP {
/*
Uses picard to filter duplicates in the on_target bam file

INPUTS
    (linked)
    < sampleId         : patient ID value
    < on_target_bam    : path to binary mapped with reads that are sorted and 
                         off target reads filtered
    < on_target_bam.bai: path to the on_target bam index

    < picard           : 

OUTPUT
(linked)
    > sampleId       : patient ID value
    > dupmark_bam    : path to binary mapped reads that are sorted, 
                       off target reads filtered, with marked suplicates
    > dupmark_bam.bai: path to the dupmark bam index

INPUT FROM
    <- sampleId,on_target_bam,on_target_bam.bai: ON_TARGET_MAPPING process
    <- picard                                  : picard channel    

OUTPUT INTO
    -> 
*/
    input:
        tuple val(sampleId), path("${sampleId}_on_target.bam"), path("${sampleId}_on_target.bam.bai") 
        path(picard)

    output:
        tuple val(sampleId), path("${sampleId}_dupmark.bam"), path("${sampleId}_dupmark.bam.bai")

    script:
    """
    java -jar ${picard} MarkDuplicates -I ${sampleId}_on_target.bam -M ${sampleId}.marked_dup.metrics.txt -O ${sampleId}_dupmark.bam
    samtools index -@ $task.cpus ${sampleId}_dupmark.bam > ${sampleId}_dupmark.bam.bai
    """
}

process COVERAGE_ANALYSIS {
/*
Uses samtools make a coverage analysis
    - Mapped reads on reference genome
    - Mapped reads that are on/off target
    - Stastistics on QC pass/fail
    - Coverage of gene panel exons

INPUTS
    (linked)
    < sampleId             : patient ID value
    < mapped_sorted_bam    : path to file containing 
                             the reads correctly mapped on the genome 
    < mapped_sorted_bam.bai: path to the mapped_sorted_bam index

    (linked)
    < sampleId             : patient ID value
    < on_target_bam        : path to binary mapped with reads that are sorted 
                             and off target reads filtered
    < on_target_bam.bai    : path to the on_target bam index
    
    (linked)
    < sampleId             : patient ID value
    < dupmark_bam          : path to binary mapped reads that are sorted, 
                             off target reads filtered, with marked duplicates
    < dupmark_bam.bai      : path to the dupmark bam index

    < report_rscript       :
    < bed_hemato           : coordinates of targets relevant to the analysis
    < bed_exon             : coordinates of exons relevant to the analysis

OUTPUT
    (linked)
    > coverage_bed
    > on_target_bed
    > bam_sort_stats
    > off_target_bam
    > coverage

INPUT FROM
    <- sampleId,mapped_sorted_bam,mapped_sorted_bam.bai: BAM_MAPPING process
    <- sampleId,on_target_bam,on_target_bam.bai        : ON_TARGET_MAPPING process
    <- sampleId,dupmark_bam,dupmark_bam.bai            : DUPMARK_BAM_SETUP process
    <- bed_hemato                                      : bed_hemato channel
    <- bed_exon                                        : bed_exon channel

OUTPUT INTO
    $params.results/$sampleId/Coverage; coverage analysis directory
*/
    publishDir "$params.results/$sampleId/Coverage", mode: 'copy'

    input:
        tuple val(sampleId), path("${sampleId}_mapped_sorted.bam"), path("${sampleId}_mapped_sorted.bam.bai")
        tuple val(sampleId), path("${sampleId}_on_target.bam"), path("${sampleId}_on_target.bam.bai")
        tuple val(sampleId), path("${sampleId}_dupmark.bam"), path("${sampleId}_dupmark.bam.bai")
        tuple env(mapped_reads), env(unmapped_reads), env(chimeric_reads)
//        path(report_rscript)
        path(bed_hemato)
        path(bed_exon)

    output:
        path "*"

    script:
    """
    samtools coverage ${sampleId}_mapped_sorted.bam -m > ${sampleId}_coverage.bed
    samtools coverage ${sampleId}_on_target.bam -m > ${sampleId}_on_target.bed
    samtools flagstat -@ $task.cpus ${sampleId}_on_target.bam > ${sampleId}_bam_sort_stats
    bedtools intersect -a ${sampleId}_mapped_sorted.bam -b ${bed_hemato} -v > ${sampleId}_off_target.bam
    bedtools coverage -a ${bed_exon} -b ${sampleId}_dupmark.bam -d > ${sampleId}_coverage
    """
//    R -e "rmarkdown::render('${report_rscript}', params = list(file='coverage',user='${USER}',pipeline='FoxNGS V1.0A',output='statistics_coverage.csv',output_gene='/media/t-chu-027/Elements/Result_NGS/Stat_gene/Statistic_couverture_gene.csv',ratio_library='${int_ratio}'),output_file='${sampleId}_couverture_analyse.bed.html')"
}

process VARSCAN {
/*
process description

INPUT
    (linked)
    < sampleId        : patient ID value
    < dupmark.bam     : path to binary mapped reads that are sorted, 
                        off target reads filtered, with marked duplicates
    < dupmark.bai     : path to the dupmark bam index

    < reference_genome: path to the reference genome (.fna file)
    < varscan         : path to the varscan software java archive (jar) file


OUTPUT
    > mutect2.vcf: variant call format (vcf) file from the mutect2 software

INPUT FROM
    <-

OUTPUT INTO
    ->
*/
    input:
        tuple val(sampleId), path("${sampleId}_dupmark.bam"), path("${sampleId}_dupmark.bam.bai")
        path(reference_genome)
        path(varscan)

    output:
        path("${sampleId}_varscan.vcf")

    script:
    """
    samtools mpileup -Q 13 -q 0 -A -B -d 100000 -f ${reference_genome} ${sampleId}_dupmark.bam -o mpileup
    java -jar ${varscan} mpileup2cns mpileup --min-coverage 50 --min-reads2 8 --min-avg-qual 30 --min-var-freq 0.02 --p-value 0.1 --strand-filter 0 --output-vcf --variants > ${sampleId}_varscan.vcf
    """
}

process MUTECT2 {
/*
process description

INPUT
    (linked)
    < sampleId        : patient ID value
    < dupmark.bam     : path to binary mapped reads that are sorted, 
                        off target reads filtered, with marked duplicates
    < dupmark.bai     : path to the dupmark bam index

    < reference_genome: path to the reference genome (.fna file)
    < indexed genome  : path to the reference genome index (.fna.fai file)
    < gatk            : path to the GATK toolkit java archive (.jar) file

OUTPUT
    > mutect2.vcf: variant call format (vcf) file from the mutect2 software

INPUT FROM
    <-

OUTPUT INTO
    ->
*/
    input:
        tuple val(sampleId), path("${sampleId}_dupmark.bam"), path("${sampleId}_dupmark.bam.bai")
        path(reference_genome)
        path(indexed_genome)
        path(reference_dict)
        path(gatk)

    output:
        path("${sampleId}_mutect2.vcf")

    script:
    """
	java -jar ${gatk} Mutect2 -R ${reference_genome} -I ${sampleId}_dupmark.bam --min-base-quality-score 30 --dont-use-soft-clipped-bases true --native-pair-hmm-threads $task.cpus -O ${sampleId}_mutect2.vcf.gz
    gunzip ${sampleId}_mutect2.vcf.gz
    """
}


process HAPLOTYPECALLER {
/*
process description

INPUT
    (linked)
    < sampleId        : patient ID value
    < dupmark.bam     : path to binary mapped reads that are sorted, 
                        off target reads filtered, with marked duplicates
    < dupmark.bai     : path to the dupmark bam index

    < dbsnp           : path to the dbsnp archive of human signe nucleotide 
                        variations and associated annotations 
    < reference_genome: path to the reference genome (.fna file)
    < gatk            : path to the GATK toolkit java archive (.jar) file

OUTPUT
    > haplotypecaller.vcf: variant call format (vcf) file from the 
                            HaplotypeCaller function in the GATK toolkit

INPUT FROM
    <- sampleId, dupmark.bam, dupmark.bai: DUPMARK_BAM_SETUP process
    <- dbsnp                             : dbsnp channel
    <- reference_genome                  : reference_genome channel
    <- gatk                              : gatk channel

OUTPUT INTO
    -> haplotypecaller.vcf: 
*/
    input:
        tuple val(sampleId), path("${sampleId}_dupmark.bam"), path("${sampleId}_dupmark.bam.bai")
        path(dbsnp)
        path(reference_genome)
        path(gatk)

    output:
        path("${sampleId}_haplotypecaller.vcf")

    script:
        """
        java -jar ${gatk} BaseRecalibrator \
            -I ${sampleId}_dupmark.bam \
            -R ${reference_genome} \
            --known-sites ${dbsnp} \
            -O ${sampleId}_recal.table


        java -jar ${gatk} ApplyBQSR \
            -I ${sampleId}_dupmark.bam \
            -R ${reference_genome} \
            -bqsr ${sampleId}_recal.table s\
            -O ${sampleId}_bqsr.bam


        java -jar ${gatk} HaplotypeCaller  \
            -I ${sampleId}_bqsr.bam \
            -R ${reference_genome} \
            --native-pair-hmm-threads $task.cpus \
            --min-base-quality-score 30 \
            --minimum-mapping-quality 20 \
            --dont-use-soft-clipped-bases true
            -O ${sampleId}_haplotypecaller.vcf \

        awk '{gsub(",Date=[^>]+>",">");}1' ${sampleId}_haplotypecaller.vcf
        """
}

/*
process PINDEL {
    input:
        tuple val(sampleId), path('dupmark.bam'), path(dupmark_bam.bai)
        path(reference_genome)
        path(bed_pindel)

    output:
        path('pindel.vcf')

    script:
        """

        """
}
*/

process ANNOVAR {
    input:
    path(vcf)
    path(annovar)
    path(humandb)
    path(python_annot)
    path(python_vcf_to_csv)
    value()

    output:
    path(variant_table)

    script:
    """
	${annovar}/table_annovar.pl ${vcf} ${annovar_db} -buildver hg${reference_version} -out annotated.vcf -remove -protocol refGene,cytoBand,cosmic92,cosmic89,avsnp138,gnomad211_exome,clinvar_20200316,dbnsfp35a,IARC,icgc21 -operation gx,r,f,f,f,f,f,f,f,f -nastring . -thread $task.cpus -polish -vcfinput -xref ${annovar}/humandb/hg19_refGene.txt

	python3 ${python_vcf_to_csv} -inVCF annotated.vcf -toType table -out variation_annotation.csv

	python3 ${python_annot} -f variation_annotation.csv -o processed_annotation.csv -m ${method}
    """
}