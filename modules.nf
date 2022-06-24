/// ----- F O X   N G S   M O D U L E S -----



// ___PROCESSES___

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
    > QC report file 

INPUT FROM
    <- setup_paired_end_reads channel

OUTPUT INTO
    -> results/$patient_id/QC directory
*/
    cpus 2
    publishDir "$params.results/$sampleId/QC"

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

process INDEX_SETUP {
/*
Runs bwa index on the reference genome

TIP: This is a long process. Consider fetching the output files in the work
directory, moving them to data/reference and setting the makeIndex parameter
to 'false' to skip this process 

INPUT
    < reference_genome

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
/*
Sets up the sam file by calling bwa mem to align the fastq on 
the reference genome.
If the path to index files has not been passed to parameters or the index files
are absent, the pipeline builds the index files with the INDEX_SETUP beforehand

INPUTS
(linked)
    < sampleID: patient ID
    < fastq_r1 and fastq_r2: NGS reads 
        (both files if doing paired-end, only fastq_r1 if doing single-end)
    
    < reference genome
    
(linked)
    < .amb, .ann, .bwt, .pac and .sa genome index files 

OUTPUT
(linked)
    < sampleId: patient ID
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
        tuple val(sampleId), path(sam)

    script:
        """
        bwa mem -t $task.cpus -R "@RG\\tID:C5-${sampleId}\\tPL:illumina\\tPU:HXXX\\tLB:Solexa\\tSM:C5-${sampleId}" ${reference_genome} ${fastq_r1} ${fastq_r2} -o sam
        """
}

process BAM_SETUP {
/*
Sets up the bam (binary version of the reads mapped to the genome) 
by calling samtools view and samtools sort

INPUTS
(linked)
    < sampleId: patient ID
    < sam     : path to the sam (Sequence Alignment Map, aligned reads) file 
OUTPUT
(linked)
    > sampleId  : patient ID
    > sorted_bam: path to the sorted bam file

INPUT FROM
    <- SAM_SETUP process

OUTPUT INTO
    -> BAM_SETUP process (sorted_bam)
    -> STDOUT (mapping stats)
*/
    input:
        tuple val(sampleId), path(sam)

    output:
        tuple val(sampleId), path(sorted_bam), emit: sorted_bam
        stdout emit: mapped_stats
        
    script:
        """
        samtools view -@ $task.cpus -Sh ${sam} -bo bam
        samtools sort -@ $task.cpus bam -o sorted_bam

        echo "Mapped: `samtools view -h -c sorted_bam` reads"
        echo "Unmapped: `samtools view -f 0x4 -h -@ $task.cpus -c -b sorted_bam` reads"
        echo "Chimeric: `samtools view  -f 0x800 -h -@ $task.cpus -c -b sorted_bam` reads" 
        """
}

process BAM_MAPPING {
/*
Filters the reads not mapped to the reference genome (ox4 tag) with samtools 
Sets up the bam index (bai) with samtools index

INPUTS
(linked)
    < sampleId  : patient ID
    < sorted_bam: path to the sorted bam file

OUTPUT
(linked)
    > sampleId             : patient ID
    > mapped_sorted_bam    : path to file containing 
                         the reads correctly mapped on the genome 
    > mapped_sorted_bam.bai: path to the mapped_sorted_bam index

INPUT FROM
    <- BAM_SETUP process

OUTPUT INTO
    -> DUPMARK_BAM_SETUP, COVERAGE_ANALYSIS processes
*/
    input:
        tuple val(sampleId), path(sorted_bam)

    output:
        tuple val(sampleId), path(mapped_sorted_bam), path('mapped_sorted_bam.bai')

    script:
        """
        samtools view -F ox4 -h -@ $task.cpus -b ${sorted_bam} > mapped_sorted_bam
        samtools index -@ $task.cpus mapped_sorted_bam > mapped_sorted_bam.bai
        """
}

process ON_TARGET_MAPPING {
/*


INPUTS
(linked)
    < sampleId             : patient ID
    < mapped_sorted_bam    : path to file containing 
                             the reads correctly mapped on the genome 
    < mapped_sorted_bam.bai: path to the mapped_sorted_bam index

    < bed_hemato       : coordinates of targets relevant to the analysis

OUTPUT
(linked)
    > sampleId     : patient ID
    > on_target_bam: path to binary mapped with reads that are sorted and 
                     off target reads filtered
    > on_target_bam.bai: path to the on_target bam index

INPUT FROM
    <- BAM_MAPPING process

OUTPUT INTO
    -> DUPMARK_BAM_SETUP, COVERAGE_ANALYSIS processes
*/
    input:
        tuple val(sampleId), path(mapped_sorted_bam), path('mapped_sorted_bam.bai')
        path(bed_hemato)

    output:
        tuple val(sampleId), path(on_target_bam), path('on_target_bam.bai')

    script:
        """
        bedtools intersect -nonamecheck -a ${mapped_sorted_bam} -b ${bed_hemato} > on_target_bam
        samtools index -@ $task.cpus on_target_bam > on_target_bam.bai
        """
}

process DUPMARK_BAM_SETUP {
/*
Uses picard to filter duplicates in the on_target bam file

INPUTS
(linked)
    < sampleId         : patient ID
    < on_target_bam    : path to binary mapped with reads that are sorted and 
                         off target reads filtered
    < on_target_bam.bai: path to the on_target bam index

OUTPUT
(linked)
    > sampleId       : patient ID
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
    tuple val(sampleId), path(on_target_bam), path('on_target_bam.bai') 
    path(picard)

    output:
    tuple val(sampleId), path(dupmark_bam), path('dupmark_bam.bai')

    script:
        """
        java -jar ${picard} MarkDuplicates -I ${on_target_bam} -M ${sampleId}.marked_dup.metrics.txt -O dupmark_bam
        samtools index -@ $task.cpus dupmark_bam > dupmark_bam.bai
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
    < sampleId             : patient ID
    < mapped_sorted_bam    : path to file containing 
                             the reads correctly mapped on the genome 
    < mapped_sorted_bam.bai: path to the mapped_sorted_bam index

(linked)
    < sampleId             : patient ID
    < on_target_bam        : path to binary mapped with reads that are sorted and 
                             off target reads filtered
    < on_target_bam.bai    : path to the on_target bam index
    
(linked)
    < sampleId             : patient ID
    < dupmark_bam          : path to binary mapped reads that are sorted, 
                             off target reads filtered, with marked duplicates
    < dupmark_bam.bai      : path to the dupmark bam index

    < bed_hemato       : coordinates of targets relevant to the analysis
    < bed_exon         : coordinates of exons relevant to the analysis

OUTPUT
(linked)
    > coverage_bed
    > on_target_bed
    > bam_sort_stats
    > off_target_bam
    > coverage

INPUT FROM
    <- sampleID,mapped_sorted_bam,mapped_sorted_bam.bai: BAM_MAPPING process
    <- sampleId,on_target_bam,on_target_bam.bai        : ON_TARGET_MAPPING process
    <- sampleId,dupmark_bam,dupmark_bam.bai            : DUPMARK_BAM_SETUP process
    <- bed_hemato                                      : bed_hemato channel
    <- bed_exon                                        : bed_exon channel

OUTPUT INTO
    $params.results/$sampleId/Coverage; coverage analysis directory
*/
    publishDir "$params.results/$sampleId/Coverage"

    input:
        tuple val(sampleId), path(mapped_sorted_bam), path('mapped_sorted_bam.bai')
        tuple val(sampleId), path(on_target_bam), path('on_target_bam.bai')
        tuple val(sampleId), path(dupmark_bam), path('dupmark_bam.bai')
//        path(report_rscript)
        path(bed_hemato)
        path(bed_exon)

    output:
        path "*"

    script:
        """
        samtools coverage ${mapped_sorted_bam} -m > coverage_bed
        samtools coverage ${on_target_bam} -m > on_target_bed
        samtools flagstat -@ $task.cpus ${on_target_bam} -O bam_sort_stats
        bedtools intersect -nonamecheck -a ${mapped_sorted_bam} -b ${bed_hemato} -v > off_target_bam
        bedtools coverage -a ${bed_exon} -b ${dupmark_bam} -d > coverage
        """
}
//      R -e "rmarkdown::render('${report_rscript}', params = list(directory='$(pwd)',file='${sampleId}_couverture_analyse.bed',user='${USER}',pipeline='${0}',output='statistics_coverage.csv',output_gene='/media/t-chu-027/Elements/Result_NGS/Stat_gene/Statistic_couverture_gene.csv',ratio_library='${int_ratio}'),output_file='$(pwd)/${sampleId}_couverture_analyse.bed.html')"

process VARSCAN {
    input:
        tuple val(sampleId), path(dupmark_bam), path('dupmark_bam.bai')
        path(varscan)
        path(reference_genome)

    output:
        path(varscan.vcf.gz)

    script:
        """
        samtools mpileup -Q 13 -q 0 -A -B -d 100000 -f ${reference_genome} ${dupmark_bam} > mpileup
        java -jar ${varscan} mpileup2cns \
        mpileup \
        --min-coverage 50 \
        --min-reads2 8 \
        --min-avg-qual 30 \
        --min-var-freq 0.02 \
        --p-value 0.1 \
        --strand-filter 0 \
        --output-vcf \
        --variants \
        > varscan.vcf.gz
        """
}

process MUTECT2 {
    input:
        tuple val(sampleId), path(dupmark_bam), path('dupmark_bam.bai')
        path(reference_genome)
        path(gatk)

    output:
        path(mutect2.vcf.gz)

    script:
        """
	    java -jar ${gatk} Mutect2 \
            -R ${reference_genome} \
            -I ${dupmark_bam} \
            --min-base-quality-score 30 \
            --dont-use-soft-clipped-bases true \
            --native-pair-hmm-threads 16 \
            -O mutect2.vcf.gz
        """
}

/*
process HAPLOTYPECALLER {
    input:
        tuple val(sampleId), path(dupmark_bam), path(dupmark_bam.bai)
        path(reference_genome)
        path(gatk)

    output:
        path(haplotypecaller.vcf.gz)

    script:
        """
        java -jar ${gatk} BaseRecalibrator \
            -I ${dupmark_bam} \
            -R ${reference_genome} \
            --known-sites $DBSNP \
            -O recal_data_table

        java -jar ${gatk} ApplyBQSR \
            -I ${dupmark_bam} \
            -R ${reference_genome} \
            -bqsr recal_data_table \
            -O bqsr_bam

        java -jar ${gatk} HaplotypeCaller  \
            -I recal_data_table \
            -R ${reference_genome} \
            --native-pair-hmm-threads 16 \
            --min-base-quality-score 30 \
            --minimum-mapping-quality 20 \
            --dont-use-soft-clipped-bases true
            -O haplotypecaller.vcf.gz \
        """
}


process PINDEL {
    input:
        tuple val(sampleId), path(dupmark_bam), path(dupmark_bam.bai)
        path(reference_genome)
        path(bed_pindel)

    output:

    script:
        """

        """
}
*/

// ___FUNCTIONS___
def setup_paired_end_reads(fastq_input) {
    fastq_input
        .map{ left, right ->  
            def sampleId = left
            def fastq_r1 = right[0]
            def fastq_r2 = right[1]
            tuple(sampleId, fastq_r1, fastq_r2)
        }
}

def setup_index(reference_index) {
    reference_index
        .map{ index ->  
            def amb = index[0]
            def ann = index[1]
            def bwt = index[2]
            def pac = index[3]
            def sa = index[4]
            tuple(amb, ann, bwt, pac, sa)
        }
}