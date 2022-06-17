PROCESSES
process QUALITY_CONTROL {
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
    input:
        tuple val(sampleId), path(sorted_bam)

    output:
        tuple val(sampleId), path(mapped_sorted_bam), path(mapped_sorted_bai)

    script:
        """
        samtools view -F ox4 -h -@ $task.cpus -b ${sorted_bam} > mapped_sorted_bam
        samtools index -@ $task.cpus mapped_sorted_bam > mapped_sorted_bai
        """
}

process ON_TARGET_MAPPING {
    input:
        tuple val(sampleId), path(mapped_sorted_bam), path(mapped_sorted_bai)
        path(bed_hemato)

    output:
        tuple val(sampleId), path(on_target_bam), path(on_target_bai)

    script:
        """
        bedtools intersect -a ${mapped_sorted_bam} -b ${bed_hemato} > on_target_bam
        samtools index -@ $task.cpus > on_target_bai
        """
}

process DUPMARK_BAM_SETUP {
    input:
    tuple val(sampleId), path(on_target_bam), path(on_target_bai) 
    path(picard)

    output:
    tuple val(sampleId), path(dupmark_bam), path(dupmark_bai)

    script:
    """
    java -jar ${picard} MarkDuplicates I=${on_target_bam} O=dupmark_bam M=${sampleId}.marked_dup.metrics.txt
    samtools index -@ $task.cpus dumpark_bam > dupmark_bai
    """${dupmark_bam}
}

process COVERAGE_ANALYSIS {
    publishDir "$params.results/$sampleId/Coverage"

    input:
        tuple val(sampleId), path(mapped_sorted_bam), path(mapped_sorted_bai)
        tuple val(sampleId), path(on_target_bam), path(on_target_bai)
        tuple val(sampleId), path(dupmark_bam), path(dupmark_bai)
        path(bed_hemato)
        path(bed_exon)

    output:
        path "*"

    script:
        """
        samtools coverage ${mapped_sorted_bam} -m > coverage_bed
        samtools coverage ${on_target_bam} -m > on_target_bed
        samtools flagstat $task.cpus ${on_target_bam} > bam_sort_stats
        bedtools intersect -a ${mapped_sorted_bam} -b ${bed_hemato} -v > off_target_bam
        bedtools coverage -a ${bed_exon} -b ${dupmark_bam} -d > couverture_
        """
}


// MODULES
def setup_reads_channel(fastq_input) {
    fastq_input
        .map{ left, right ->  
            def sampleId = left
            def fastq_r1 = right[0]
            def fastq_r2 = right[1]
            tuple(sampleId, fastq_r1, fastq_r2)
        }
}

def setup_index_channel(reference_index) {
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