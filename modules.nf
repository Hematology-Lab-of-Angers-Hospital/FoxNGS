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
    < sampleId             : patient ID value
    < fastq_r1 and fastq_r2: NGS reads 
        (both files if doing paired-end, only fastq_r1 if doing single-end)

OUTPUT
    (linked)
    > sampleId          : patient ID value
    > R1_001.fastqc.html: quality control metrics web page for the first end of paired-end reads
    > R2_001.fastqc.html: quality control metrics web page for the other end of paired-end reads

INPUT FROM
    <- setup_paired_end_reads channel

OUTPUT INTO
    -> results/$patient_id/QC directory
*/
    tag("${sampleId}")

    cpus 2
    publishDir "${params.results}/${sampleId}/QC", mode: 'copy'

    input:
        tuple val(sampleId), path(fastq_r1), path(fastq_r2)

    output:
        path('*')

    script:
    """
    fastqc -t ${task.cpus} ${fastq_r1} ${fastq_r2}
    """
}

process FAI_SETUP {
/*
Runs samtools faidx to make the reference genome index (.fai file)
if that is necessary

TIP: This process is not as long as REFERENCE_INDEX_SETUP but the pipeline has parameters
to fetch the .fai file from the data/reference folder; consider copying the
.fai output file and set param.indexed_genome to 'true' to skip this step.

INPUT
    < reference_genome: path to the reference genome (.fna file)

OUTPUT
    > path to the reference genome index (.fna.fai file)

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

TIP: This process is not as long as REFERENCE_INDEX_SETUP but the pipeline has parameters
to fetch the .fai file from the data/reference folder; consider copying the
.fai output file and set param.indexed_genome to 'true' to skip this step.

INPUT
    < reference_genome: path to the reference genome (.fna file)
    < gatk            : path to the GATK toolkit java archive (.jar) file

OUTPUT
    > path to the reference genome dictionary (.dict file)

INPUT FROM 
    <- reference_genome: reference_genome value channel
    <- gatk            : gatk value channel

OUTPUT INTO
    -> reference_dict channel
*/
    cpus 2

    input:
        path(gatk)
        path(reference_genome)
        
    output:
        path("*.dict")

    script:
    """
    java -jar ${gatk} CreateSequenceDictionary \
        -R ${reference_genome}
    """
}

process BED_INTERVAL_SETUP {
/*
Uses the Picard BedToInterval function to convert a bed file 
to an interval list that serves an an input for the CollectHsMetrics in Picard 

INPUT
    < picard        : path to the Picard java archive file (.jar)
    < bed           : path to a bed file
    < reference_dict: path to the reference genome dictionary (.dict file)


OUTPUT
    > path to the bed converted to an interval list

INPUT FROM
    <- picard        : picard value channel
    <- bed           : bed_exon OR bed_bait value channel
    <- reference_dict: REFERENCE_DICT_SETUP process

OUTPUT INTO
    -> COLLECT_HS_METRICS (The process is called twice as BAIT_BED_INTERVAL_SET and EXON_BED_INTERVAL_SET
    Both processes outputs are passed to COLLECT_HS_METRICS)
*/
    input:
        path(picard)
        path(bed)
        val(filetype)
        path(reference_dict)

    output:
        path("${filetype}_list.interval_list")

    script:
    """
    java -jar ${picard} BedToIntervalList \
        -I ${bed} \
        -O ${filetype}_list.interval_list \
        -SD ${reference_dict}
    """
}

process REFERENCE_INDEX_SETUP {
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
    <- reference_genome value channel

OUTPUT INTO
    -> REFERENCE_INDEX_SETUP process (if make_index parameter is set on true)
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

process VARIATION_INDEX_SETUP {
/*
process description

INPUT
    < gatk  : path to the GATK toolkit java archive (.jar) file
    < dbsnp : path to the dbsnp archive of human signe nucleotide 
              variations and associated annotations

OUTPUT
    > path to the .idx dbsnp index file

INPUT FROM
    <- gatk : gatk value channel
    <- dbsnp: dbsnp value channel

OUTPUT INTO
    -> HAPLOTYPECALLER process
*/
    cpus 2

    input:
        path(gatk)
        path(dbsnp)

    output:
        path('*.idx')

    script:
    """
    java -jar ${gatk} IndexFeatureFile \
        -I ${dbsnp}
    """
}

process ANNOTATION_DICTIONNARY_SETUP {
/*
Calls annotation_process to produce a .json dictionnary file of known
artifacts and favorite transcripts. Allows to tailor
the filtering process, removing known sequencing artifacts 
and selecting favorite transcripts.

Used in StageInMode 'copy' to alleviate for file opening issues
with symlinks in the script.

INPUT
    < python_annot    : Path to the python script building a .json dictionary from
                        an artifact and a transcript base
    < transcript_base : File containing known transcripts
    < artifact_base   : File containing known sequencing artifacts in the lab
    < pipeline_version: current version of the pipelin, used to name the dictionary

OUTPUT
    > variant_dictionary_routine_v${pipeline_version}.json: path to the .idx dbsnp index file

INPUT FROM
    <- python_annot    : python_annot value channel
    <- transcript_base : transcript_base value channel
    <- artifact_base   : artifact_base value channel
    <- pipeline_version: pipeline_version value channel

OUTPUT INTO
    -> MERGE_ANNOVAR_FILES process
*/
    stageInMode "copy"

    input:
        path(python_annot)
        path(artifact_base)
        path(transcript_base)
        val(pipeline_version)

    output:
        path("variant_dictionary_routine_v${pipeline_version}.json")

    script:
    """
    python ${python_annot} -a ${artifact_base} -t ${transcript_base} -c variant_dictionary_routine_v${pipeline_version}.json
    """
}

process BAM_SETUP {
/*
Sets up the bam (binary version of the reads mapped to the genome) 
by calling samtools view and samtools sort. Both are piped together
to greatly optimize computing time.

Now uses 16 CPUs per process to compensate for alignment bias of bwa-mem
when CPU allocation changes. 
Legacy pipeline used 16 CPUs, so will FoxNGS.

INPUTS
    (linked)
    < sampleId: patient ID value
    < fastq_r1: path to the r1 fastq.gz file
    < fastq_r2: path to the r2 fastq.gz file

OUTPUT
    (linked)
    > sampleId  : patient ID value
    > sorted_bam: path to the bam (Binary Alignment Map, aligned reads in binary setup) file 
                  with reads sorted by starting position

INPUT FROM
    <- fastq_r1, fastq_r2 : fastq_input channel

OUTPUT INTO
    -> sampleId, sorted.bam: MAPPED_BAM_SETUP process (sorted_bam)
    -> mapping_stats       : STDOUT
*/
    tag "${sampleId}"
    cpus 16

    input:
        tuple val(sampleId), path(fastq_r1), path(fastq_r2)
        path(reference_genome)
        tuple path(amb), path(ann), path(bwt), path(pac), path(sa)

    output:
        tuple val(sampleId), path("${sampleId}_sorted.bam")
        
    script:
    """
    bwa mem -t $task.cpus -j -R "@RG\\tID:C5-${sampleId}\\tPL:illumina\\tPU:HXXX\\tLB:Solexa\\tSM:C5-${sampleId}" ${reference_genome} ${fastq_r1} ${fastq_r2} | samtools sort -@ ${task.cpus} -O BAM -o ${sampleId}_sorted.bam

    mapped_reads=`samtools view -h -c ${sampleId}_sorted.bam`
    unmapped_reads=`samtools view -f 0x4 -h -@ $task.cpus -c -b ${sampleId}_sorted.bam`
    chimeric_reads=`samtools view  -f 0x800 -h -@ $task.cpus -c -b ${sampleId}_sorted.bam`
    """
}

process MAPPED_BAM_SETUP {
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
    tag "${sampleId}"
    cpus 4
    publishDir "${params.results}/${sampleId}", mode: 'copy'

    input:
        tuple val(sampleId), path("${sampleId}_sorted.bam")

    output:
        tuple val(sampleId), path("${sampleId}_mapped_sorted.bam"), path("${sampleId}_mapped_sorted.bam.bai")

    script:
    """
    samtools view -F 0x4 -h -@ $task.cpus -b ${sampleId}_sorted.bam > ${sampleId}_mapped_sorted.bam

    samtools index -@ $task.cpus ${sampleId}_mapped_sorted.bam > ${sampleId}_mapped_sorted.bam.bai
    """
}

process IN_TARGET_MAPPING {
/*
Filters reads by position (relative to reference geneme mapping), keeping only
position that are clinically relevant
Sets up the bam index (bai) with samtools index

INPUTS
(linked)
    < sampleId             : patient ID value
    < mapped_sorted_bam    : path to file containing 
                             the reads correctly mapped on the genome 
    < mapped_sorted_bam.bai: path to the mapped_sorted_bam index

    < bed_bait              : coordinates of sequences relevant to the analysis

OUTPUT
(linked)
    > sampleId     : patient ID value
    > on_target_bam: path to binary mapped sorted reads,
                     with off target reads filtered
    > on_target_bam.bai: path to the on_target bam index

INPUT FROM
    <- MAPPED_BAM_SETUP process

OUTPUT INTO
    -> DUPMARK_BAM_SETUP, COVERAGE_ANALYSIS processes
*/
    tag "${sampleId}"
    cpus 2

    input:
        tuple val(sampleId), path(mapped_sorted_bam), path("${sampleId}_mapped_sorted.bam.bai")
        path(bed_bait)

    output:
        tuple val(sampleId), path("${sampleId}_in_target.bam"), path("${sampleId}_in_target.bam.bai")

    script:
    """
    bedtools intersect -a ${sampleId}_mapped_sorted.bam -b ${bed_bait} > ${sampleId}_in_target.bam

    samtools index -@ $task.cpus ${sampleId}_in_target.bam > ${sampleId}_in_target.bam.bai
    """
}

process DUPMARK_BAM_SETUP {
/*
Uses picard to filter duplicates in the on_target bam file.

As long as GATK variant calling tools do not make a difference 
in duplicate marking (either modify hexadecimal flag of duplicates or 
consider additional duplicate flag in duplicate filtering) between optical duplicates 
and sequencing duplicates due to over amplification in targeted panels, 
variant callers will have to ignore duplicate marking results in order 
to increase reference/mutant reporting.

INPUTS
    (linked)
    < sampleId         : patient ID value
    < on_target_bam    : path to binary mapped sorted reads,
                         with off target reads filtered
    < on_target_bam.bai: path to the on_target bam index

    < picard           : path to the Picard java archive file (.jar)

OUTPUT
(linked)
    > sampleId       : patient ID value
    > dupmark_bam    : path to binary mapped reads that are sorted, 
                       off target reads filtered, with marked duplicates
    > dupmark_bam.bai: path to the dupmark_bam index

INPUT FROM
    <- sampleId,on_target_bam,on_target_bam.bai: ON_TARGET_MAPPING process
    <- picard                                  : picard value channel

OUTPUT INTO
    -> COVERAGE_ANALYSIS, VARSCAN, MUTECT2, HAPLOTYPECALLER, PINDEL processes
*/
    tag "${sampleId}"
    cpus 4
    publishDir "${params.results}/${sampleId}/Coverage", mode: 'copy', pattern: "*.marked_dup.metrics.txt"

    input:
        tuple val(sampleId), path("${sampleId}_in_target.bam"), path("${sampleId}_in_target.bam.bai") 
        path(picard)

    output:
        tuple val(sampleId), path("${sampleId}_dupmark.bam"), path("${sampleId}_dupmark.bam.bai"), emit: bam
        path("${sampleId}.marked_dup.metrics.txt"), emit: metrics

    script:
    """
    java -Xmx4g -jar ${picard} MarkDuplicates \
        -I ${sampleId}_in_target.bam \
        -M ${sampleId}.marked_dup.metrics.txt \
        --REMOVE_SEQUENCING_DUPLICATES true \
        -O ${sampleId}_dupmark.bam
    
    samtools index -@ $task.cpus ${sampleId}_dupmark.bam > ${sampleId}_dupmark.bam.bai
    """
}

process BQSR_BAM_SETUP {
/*
Calls GATK to build a model of base quality covariation, 
to subsequently apply that model to recalibrate quality
in the bam file

Uses .vcf of known sites (such as dbsnp) to not affect these positions.

DO NOT USE ON TARGETED PANELS, as these have not enough variation to
provide a covariate model that has enough relevance for base revalirbation, 
uselessly driving down the average base quality.

INPUT
    < gatk            : path to the GATK toolkit java archive (.jar) file

    (linked)
    < sampleId        : patient ID value
    < dupmark_bam     : path to binary mapped reads that are sorted, 
                        off target reads filtered, with marked duplicates
    < dupmark_bam.bai : path to the dupmark_bam index

    < reference_genome: path to the reference genome (.fna file)
    < indexed_genome  : path to the reference genome index (.fna.fai file)
    < reference_dict  : path to the reference genome dictionary (.dict file)
    < dbsnp           : path to the reference database of SNP
    < dbsnp_idx       : path to the dbsnp index

OUTPUT
    (linked)
    > sampleId     : patient ID value
    > bqsr.bam     : path to the base recalibrated bam file
    > bqsr.bam.bam : path to the bqsr.bam index

    (linked)
    > sampleId    : patient ID value
    > recal_table : table of base recalibrations, produced by the 
                    BaseRecalibrator GATK function

INPUT FROM
    <- sampleId, dupmark.bam, dupmark.bam.bai : DUPMARK_BAM_SETUP process
    <- indexed_genome                         : FAI_SETUP process
    <- reference_dict                         : REFERENCE_DICT_SETUP process
    <- dbsnp_idx                              : VARIATION_INDEX_SETUP process

    <- gatk                                   : gatk channel
    <- reference_genome                       : reference_genome channel
    <- dbsnp                                  : dbsnp channel

OUTPUT INTO
    -> sampleId, bqsr.bam, bqsr.bam.bai : HAPLOTYPECALLER process
    -> sampleId, recal_table            : PATIENT_REPORT process
*/
    input:
    path(gatk)
    tuple val(sampleId), path("${sampleId}_dupmark.bam"), path("${sampleId}_dupmark.bam.bai")
    path(reference_genome)
    path(indexed_genome)
    path(reference_dict)
    path(dbsnp)
    path(dbsnp_idx)

    output:
    tuple val(sampleId), path("${sampleId}_bqsr.bam"), path("${sampleId}_bqsr.bam.bai"), emit: bam
    tuple val(sampleId), path("${sampleId}_recal.table"), emit: base_recalibration

    script:
    """
    java -Xmx4g -jar ${gatk} BaseRecalibrator \
        -I ${sampleId}_dupmark.bam \
        -R ${reference_genome} \
        --known-sites ${dbsnp} \
        -O ${sampleId}_recal.table

    java -Xmx2g -jar ${gatk} ApplyBQSR \
        -I ${sampleId}_dupmark.bam \
        -R ${reference_genome} \
        -bqsr ${sampleId}_recal.table \
        -O ${sampleId}_bqsr.bam

    samtools index -@ $task.cpus ${sampleId}_bqsr.bam > ${sampleId}_bqsr.bam.bai
    """
}

process MPILEUP_SETUP{
/*
Calls samtools on the bam file (with marked duplicates) 
to produce a mpileup file

INPUT
    (linked)
    < sampleId       : patient ID value
    < dupmark_bam    : path to binary mapped reads that are sorted, 
                       off target reads filtered, with marked duplicates
    < dupmark_bam.bai: path to the dupmark bam index


OUTPUT
    (linked)
    > sampleId: patient ID value
    > mpileup : 

INPUT FROM
    <- sampleId, dupmark.bam, dupmark.bam.bai: DUPMARK_BAM_SETUP process

    <- reference_genome : reference_gnome channel

OUTPUT INTO
    -> sampleId, mpileup : VARSCAN process
*/
    input:
        tuple val(sampleId), path("${sampleId}_dupmark.bam"), path("${sampleId}_dupmark.bam.bai")
        path(reference_genome)

    output:
        tuple val(sampleId), path("${sampleId}.mpileup")

    script:
    """
    samtools mpileup -Q 13 -q 0 -A -B -d 100000 -f ${reference_genome} ${sampleId}_dupmark.bam -o ${sampleId}.mpileup
    """
// ????????
}

process COLLECT_HS_METRICS {
/*
Uses the CollectHSMetrics funtion from the Picard suite to gather multiple
metrics on coverage, mapping... The output is used in the multiqc report
The NEAR_DISTANCE is set to 0 so that the selected_reads percentage
metrics in the output reflects an "in bait"/"on target" metric.

INPUT
    < picard               : path to the Picard java archive file (.jar)\

(linked)
    < sampleId             : patient ID value
    < mapped_sorted.bam    : path to file containing 
                             the reads correctly mapped on the genome 
    < mapped_sorted.bam.bai: path to the mapped_sorted_bam index

    < reference_genome     : path to the reference genome (.fna file)
    < reference_fai        : path to the reference genome index (.fna.fai file)
    < design_interval_list : path to the bait bed converted to an interval list
    < exon_interval_list   : path to the target exon bed converted to an interval list

OUTPUT
    > output_hs_metrics  : table of metrics collected by the CollectHsMetrics tool

INPUT FROM
    <- sampleId, mapped_sorted.bam, mapped_sorted.bam.bai: MAPPED_BAM_SETUP process
    <- reference_fai                                     : FAI_SETUP process
    <- design_interval_list, exon_interval_list          : BED_INTERVAL_SETUP process

    <- picard                                            : picard value channel
    <- reference_genome                                  : reference_genome value channel

OUTPUT INTO
    -> $params.results/$sampleId/Coverage: coverage analysis directory
*/
    tag("${sampleId}")
    publishDir "${params.results}/${sampleId}/Coverage", mode: 'copy'

    input:
        path(picard)
        tuple val(sampleId), path("${sampleId}_mapped_sorted.bam"), path("${sampleId}_mapped_sorted.bam.bai")
        path(reference_genome)
        path(reference_fai)
        path(design_interval_list)
        path(exon_interval_list)

    output:
        tuple val(sampleId), path("${sampleId}_output_hs_metrics.txt")

    script:
    """
    java -jar ${picard} CollectHsMetrics \
        I=${sampleId}_mapped_sorted.bam \
        O=${sampleId}_output_hs_metrics.txt \
        R=${reference_genome} \
        NEAR_DISTANCE=0 \
        BAIT_INTERVALS=${design_interval_list} \
        TARGET_INTERVALS=${exon_interval_list}
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

    < bed_bait             : coordinates of targets relevant to the analysis
    < bed_exon             : coordinates of exons relevant to the analysis

OUTPUT
    (linked)
    > 
    _in_target

INPUT FROM
    <- sampleId,mapped_sorted_bam,mapped_sorted_bam.bai: MAPPED_BAM_SETUP process
    <- sampleId,on_target_bam,on_target_bam.bai        : ON_TARGET_MAPPING process
    <- sampleId,dupmark_bam,dupmark_bam.bai            : DUPMARK_BAM_SETUP process
    <- bed_bait                                        : bed_bait value channel
    <- bed_exon                                        : bed_exon value channel

OUTPUT INTO
    $params.results/$sampleId/Coverage; coverage analysis directory
*/  
    tag "${sampleId}"
    cpus 2
    publishDir "${params.results}/${sampleId}/Coverage", mode: 'copy', pattern: "*.mosdepth.*"

    input:
        tuple val(sampleId), path("${sampleId}_mapped_sorted.bam"), path("${sampleId}_mapped_sorted.bam.bai"), path("${sampleId}_in_target.bam"), path("${sampleId}_in_target.bam.bai")
        path(bed_bait)
        path(bed_exon)

    output:
        path '*'
        tuple val(sampleId), path("${sampleId}_exon.per-base.bed.gz"), path("${sampleId}_in_target.mosdepth.region.dist.txt"), path("${sampleId}_exon.mosdepth.region.dist.txt"), path("${sampleId}_exon.regions.bed.gz"), emit: mosdepth_stats

    script:
    """
    samtools flagstat -@ $task.cpus ${sampleId}_in_target.bam > ${sampleId}_in_target_stats
    
    mosdepth ${sampleId}_in_target -b ${bed_bait} ${sampleId}_mapped_sorted.bam -t $task.cpus

    mosdepth ${sampleId}_exon -b ${bed_exon} ${sampleId}_in_target.bam -t $task.cpus
    """
}


process PATIENT_REPORT {
/*
Intersects a mosdepth exon.per_base coverage file with a hotspots postion file
to check if these hotspots are covered above 200x. Issues a warning in the file
named hotspots_warning.txt if any base in the coverage bed has a hospot position
with a coverage value under 200x.

Used in StageInMode 'copy' to alleviate for file opening issues
with symlinks when using MultiQC.

Runs in a few seconds thanks to low level scripting.

INPUT
    (linked)
    > sampleId             : patient ID value
    > exon.per-base.bed.gz : path to the file containing coverage of exons with
                             1 bp resolution

OUTPUT
    > hotspot_warnings.txt

INPUT FROM
    <- sampleId, exon.per-base.bed.gz: COVERAGE_ANALYSIS process 
    <- bed_hotspots                  : bed_hotspots value channel

OUTPUT INTO
    -> $params.results/$sampleId/Coverage; coverage analysis directory
*/
    tag "${sampleId}"
    publishDir "${params.results}/${sampleId}", mode: 'copy'
    stageInMode 'copy'

    input:
        tuple val(sampleId), path("${sampleId}_exon.per-base.bed.gz"), path("${sampleId}_in_target.mosdepth.region.dist.txt"), path("${sampleId}_exon.mosdepth.region.dist.txt"), path("${sampleId}_exon.regions.bed.gz"), path("${sampleId}_output_hs_metrics.txt")
        path(bed_hotspots)
        path(bed_exon)
        path(exon_template)
        path(hotspot_template)
        path(patient_report_config)

    output:
        path("${sampleId}_report.html")

    script:
    """
    cat ${hotspot_template} > ${sampleId}_hotspot_coverage_mqc.csv
    cat ${exon_template} > ${sampleId}_exon_coverage_mqc.csv

    gunzip ${sampleId}_exon.per-base.bed.gz

    bedtools intersect -a ${bed_hotspots} -b ${sampleId}_exon.per-base.bed -wb | cut -f 2,3,4,8 | awk -F '\t' -v OFS=',' '{ if(\$4 < 200) print NR,\$3,\$1,\$2,\$4 }' | sed 's/,/-/3' >> ${sampleId}_hotspot_coverage_mqc.csv
    bedtools intersect -a ${sampleId}_exon.per-base.bed -b ${bed_exon} -wb | cut -f 2,3,4,8 | awk -F '\t' -v OFS=',' '{if(\$3 < 200) print \$4}' | uniq -c | sed -e 's/^[ \t]*//' | sed -e "s/ /,/g" | awk -F ',' -v OFS=',' '{print \$2,\$1}' >> ${sampleId}_exon_coverage_mqc.csv

    multiqc . -c ${patient_report_config} -n ${sampleId}_report.html
    """
}


process CONTROL_FREEC {
/*
INPUT
<

OUTPUT
>

INPUT FROM
<-

OUTPUT INTO
->
*/
    tag "$sampleId"
    publishDir "${params.results}/${sampleId}/Variation/CNV", mode: 'copy'

    input:
        tuple val(sampleId), path("${sampleId}.mpileup")
        path(freec_config)
        path(R_freec_graphics)

    output:
        path("*")

    script:
    """
    freec -conf ${freec_config} -sample ${sampleId}.mpileup
    cat ${R_freec_graphics} | R --slave --args ${sampleId}.mpileup_ratio.txt ${sampleId}.mpileup_BAF.txt $sampleId
    """
}


process VARSCAN {
/*
Calls genome variation in the .bam files (comparing it to the reference genome)
using varscan

INPUT
    (linked)
    < sampleId        : patient ID value
    < dupmark.bam     : path to binary mapped reads that are sorted, 
                        off target reads filtered, with marked duplicates
    < dupmark.bai     : path to the dupmark bam index

    < reference_genome: path to the reference genome (.fna file)
    < varscan         : path to the varscan software java archive (.jar) file


OUTPUT
    (linked)
    > sampleId   : patient ID value
    > varscan.vcf: variant call format (.vcf) file from the varscan software

INPUT FROM
    <- sampleId,dupmark.bam,dupmark.bai: DUPMARK_BAM_SETUP process
    <- reference_genome                : reference_genome value channel
    <- varscan                         : varscan value channel

OUTPUT INTO
    -> ANNOVAR process
*/  
    tag "${sampleId}"
    cpus 2

    input:
        tuple val(sampleId), path("${sampleId}.mpileup")
        path(varscan)

    output:
        tuple val(sampleId), val("varscan"), path("${sampleId}_varscan.vcf"), emit: varscan_variation

    script:
    """
    java -jar ${varscan} mpileup2cns \
        ${sampleId}.mpileup \
        --min-coverage 50 \
        --min-reads2 8 \
        --min-avg-qual 30 \
        --min-var-freq 0.0175 \
        --p-value 0.1 \
        --strand-filter 0 \
        --output-vcf \
        --variants > ${sampleId}_varscan.vcf
    """
}


process MUTECT2 {
/*
Calls genome variation in the .bam files (comparing it to the reference genome)
using Mutect2 from the GATK software suite.
Duplicate reads filter disabled, see DUPMARK_BAM_SETUP process
description for more information.

INPUT
    < gatk            : path to the GATK toolkit java archive (.jar) file

    (linked)
    < sampleId        : patient ID value
    < dupmark.bam     : path to binary mapped reads that are sorted, 
                        off target reads filtered, with marked duplicates
    < dupmark.bai     : path to the dupmark bam index

    < reference_genome: path to the reference genome (.fna file)
    < indexed genome  : path to the reference genome index (.fna.fai file)
    < reference_dict  : path to the reference genome dictionary (.dict file)

OUTPUT
    (linked)
    > sampleId   : patient ID value
    > mutect2.vcf: variant call format (.vcf) file from the mutect2 software

INPUT FROM
    <- sampleId,dupmark.bam,dupmark.bam.bai: DUPMARK_BAM_SETUP process
    <- reference_genome                    : reference_genome value channel
    <- indexed_genome                      : REFERENCE_INDEX_SETUP
    <- reference_dict                      : REFERENCE_DICT_SETUP process

OUTPUT INTO
    -> ANNOVAR process OR VEP process
*/
    tag "${sampleId}"
    cpus 2

    input:
        path(gatk)
        tuple val(sampleId), path("${sampleId}_bqsr.bam"), path("${sampleId}_bqsr.bam.bai")
        path(reference_genome)
        path(indexed_genome)
        path(reference_dict)

    output:
        tuple val(sampleId), val("mutect2"), path("${sampleId}_mutect2.vcf")

    script:
    """
    java -jar ${gatk} Mutect2 \
        -I ${sampleId}_bqsr.bam \
        -R ${reference_genome} \
        --min-base-quality-score 30 \
        --dont-use-soft-clipped-bases true \
        --disable-read-filter NotDuplicateReadFilter \
        --native-pair-hmm-threads 16 \
        -O ${sampleId}_mutect2.vcf.gz

    java -jar ${gatk} FilterMutectCalls \
        -V ${sampleId}_mutect2.vcf.gz \
        -R ${reference_genome} \
        -O ${sampleId}_mutect2.vcf

	java -jar ${gatk} Mutect2 \
        -I ${sampleId}_dupmark.bam \
        -R ${reference_genome} \
        --min-base-quality-score 30 \
        --dont-use-soft-clipped-bases true \
        --disable-read-filter NotDuplicateReadFilter \
        --native-pair-hmm-threads 16 \
        -L chr5:170837531-170837888 \
        -O ${sampleId}_NPM1_mutect2.vcf.gz

    gunzip ${sampleId}_NPM1_mutect2.vcf.gz

    awk '!/^#+/' ${sampleId}_NPM1_mutect2.vcf >> ${sampleId}_mutect2.vcf
    """
}


process HAPLOTYPECALLER {
/*
Calls genome variation in the .bam files (comparing it to the reference genome)
using HaplotypeCaller from the GATK software suite. Uses dbsnp to filter SNP
that are reported as non-pathogenic.
Duplicate reads filter disabled, see DUPMARK_BAM_SETUP process
description for more information.

INPUT
    < gatk            : path to the GATK toolkit java archive (.jar) file

    (linked)
    < sampleId        : patient ID value
    < dupmark.bam     : path to binary mapped reads that are sorted, 
                        off target reads filtered, with marked duplicates
    < dupmark.bai     : path to the dupmark bam index

    < reference_genome: path to the reference genome (.fna file)
    < indexed_genome  : path to the reference genome index (.fna.fai file)
    < reference_dict  : path to the reference genome dictionary (.dict file)
    < dbsnp           : path to the dbsnp archive of human signe nucleotide 
                        variations and associated annotations
    < dbsnp_idx       : path 


OUTPUT
    (linked)
    > sampleId           : patient ID value
    > haplotypecaller.vcf: variant call format (.vcf) file from the 
                           HaplotypeCaller function in the GATK toolkit

INPUT FROM
    <- sampleId, dupmark.bam, dupmark.bai: DUPMARK_BAM_SETUP process
    <- reference_genome                  : reference_genome value channel
    <- indexed_genome                    : FAI_SETUP process
    <- reference_dict                    : REFERENCE_DICT_SETUP process
    <- dbsnp                             : dbsnp value channel
    <- dbsnp_idx                         : VARIATION_INDEX_SETUP process

OUTPUT INTO
    -> ANNOVAR process OR VEP process
*/
    tag "${sampleId}"
    publishDir "${params.results}/${sampleId}/Coverage", mode: 'copy', pattern: "*_recal.table"
    cpus 2

    input:
        path(gatk)
        tuple val(sampleId), path("${sampleId}_bqsr.bam"), path("${sampleId}_bqsr.bam.bai")
        path(reference_genome)
        path(indexed_genome)
        path(reference_dict)

    output:
        tuple val(sampleId), val("gatk"), path("${sampleId}_haplotypecaller.vcf"), emit: gatk_variation

    script:
    """
    java -jar ${gatk} HaplotypeCaller \
        -I ${sampleId}_bqsr.bam \
        -R ${reference_genome} \
        --min-base-quality-score 30 \
        --minimum-mapping-quality 20 \
        --native-pair-hmm-threads $task.cpus \
        --disable-read-filter NotDuplicateReadFilter \
        --dont-use-soft-clipped-bases true \
        -O ${sampleId}_haplotypecaller.vcf

    awk '{gsub(",Date=[^>]+>",">");}1' ${sampleId}_haplotypecaller.vcf
    """
}


process PINDEL {
/*
Calls genome variation in the .bam files (comparing it to the reference genome)
using pindel to find specific InDel in CALR and FLT3 (hence the bed input).

Used in StageInMode 'copy' to alleviate for file opening issues
with symlinks in the Pindel config file.

INPUT
    (linked)
    < sampleId             : patient ID value
    < dupmark.bam          : path to binary mapped reads that are sorted, 
                             off target reads filtered, with marked duplicates
    < dupmark.bai          : path to the dupmark bam index

    (linked)
    > sampleId             : patient ID value
    > mapped_sorted.bam    : path to file containing 
                             the reads correctly mapped on the genome 
    > mapped_sorted.bam.bai: path to the mapped_sorted_bam index

OUTPUT
    (linked)
    > sampleId  : patient ID value
    > pindel.vcf: variant call format (vcf) file from the pindel variant caller

INPUT FROM
    <- sampleId, dupmark.bam, dupmark.bai            : DUPMARK_BAM_SETUP process
    <- sampleId, mapped_sorted.bam, mapped_sorted.bai: MAPPED_BAM_SETUP process
    <- reference_genome                              : reference_genome value channel
    <- pindel                                        : pindel value channel
    <- bed_pindel                                    : bed_pindel value channel

OUTPUT INTO
    -> ANNOVAR process
*/
    tag "${sampleId}"
    stageInMode "copy"
    afterScript 'rm *.bam*; rm *.fna; rm *_ITD'

    input:
        tuple val(sampleId), path("${sampleId}_bqsr.bam"), path("${sampleId}_bqsr.bam.bai")
        tuple val(sampleId), path("${sampleId}_mapped_sorted.bam"), path("${sampleId}_mapped_sorted.bam.bai")
        path(reference_genome)
        path(reference_fai)
        path(pindel)
        path(bed_pindel)

    output:
        tuple val(sampleId), val("pindel"), path("${sampleId}_pindel.vcf")

    script:
    """
    echo -e "${sampleId}_bqsr.bam 800 Duplicate_mark\n${sampleId}_mapped_sorted.bam 800 All_read" > config_file_pindel.txt
    
    ${pindel}/pindel -f ${reference_genome} -i config_file_pindel.txt -j ${bed_pindel} -T $task.cpus -o ${sampleId}_ITD

    ${pindel}/pindel2vcf -P ${sampleId}_ITD -r ${reference_genome} -R x -d 00000000 -G -v ${sampleId}_pindel.vcf
    """
}


process MERGE_VCF {
    publishDir "${params.results}/${sampleId}/Variation", mode: 'copy', pattern: 'Annotation_*.csv'

    input:
        tuple val(sampleId), val(mutect2), path("${sampleId}_mutect2.vcf"), val("gatk"), path("${sampleId}_gatk.vcf")
        path(picard)
        path(reference_genome)
        path(reference_dict)

    output:
        tuple val(sampleId), path("${sampleId}_merged_dedup.vcf") 

    script:
    """
    java -jar ${picard} MergeVcfs \
        -I ${sampleId}_mutect2.vcf \
        -I ${sampleId}_gatk.vcf \
        -R ${reference_genome} \
        -D ${reference_dict} \
        -O ${sampleId}_merged.vcf
    
    LC_ALL=C sort -t \$'\t' -k1,1 -k2,2n -k4,4 ${sampleId}_merged.vcf | \
    awk -F '\t' '/^#/ {print;prev="";next;} {key=sprintf("%s\t%s\t%s",\$1,\$2,\$4);if(key==prev) next;print;prev=key;}' > ${sampleId}_merged_dedup.vcf
    """
}


process ANNOVAR {
/*
Annotates the vcf files from the variant calling processes using ANNOVAR and the following databases:
    refGene, cytoBand, cosmic92, cosmic89, avsnp138, gnomad211_exome, clinvar_20200316, dbnsfp35a, IARC, icgc21
Simplifies the ANNOVAR output from a vcf file to a csv, 
Applies filtering to get clinically significant variants 

INPUT
    < annovar          : path to the ANNOVAR variant annotation software 
    < reference_version: version of the referenge genome (hg19/hg38)

    (linked)
    < sampleId         : patient ID value
    < method           : method used to annotate the .bam file 
                         (Varscan, Mutect2, Hapotypecaller (GATK) or Pindel)
    < vcf              : path to the vcf corresponding to the sampleId and
                         method used

    < python_vcf_to_csv: path to the python script to simplify the vcf to a csv
    < python_annot     : path to the python script filtering the csv table

OUTPUT
    (linked)
    > sampleId                                           : patient ID value
    > method                                             : method used to annotate the .bam file 
                                                           (Varscan, Mutect2, Hapotypecaller (GATK) or Pindel)
    > Filter_simple_annotation_${sampleId}_${method}.csv : comma separated values (.csv) file containing variants
                                                           annotated and filtered

INPUT FROM
    <- annovar              : annovar value channel
    <- reference_version    : reference_version value
    <- sampleId, method, vcf: VARSCAN, MUTECT2, HAPLOTYPECALLER, PINDEL processes
        (This module is called 4 times as ANNOVAR_VSC, ANNOVAR_MT2, ANNOVAR_HTC and ANNOVAR_PDL)
    <- python_vcf_to_csv    : python_vcf_to_csv value channel
    <- python_annot         : python_annot value channel

OUTPUT INTO
    -> MERGE_ANNOVAR_FILES
*/
    tag "${sampleId}"
    publishDir "${params.results}/${sampleId}/Variation", mode: 'copy', pattern: "Filter_simple_annotation*"
    errorStrategy 'ignore'

    input:
        path(annovar)
        val(reference_version)
        tuple val(sampleId), val(method), path(vcf)
        path(humandb_annovar)
        path(python_vcf_to_csv)
        path(python_annot)

    output:
        tuple val(sampleId), val(method), path("Filter_simple_annotation_${sampleId}_${method}.csv")

    script:
    """
    ${annovar}/table_annovar.pl ${vcf} ${humandb_annovar} -buildver hg${reference_version} -out ${sampleId}_${method}_annotated -remove -protocol refGene,cytoBand,cosmic92,cosmic89,avsnp138,gnomad211_exome,clinvar_20200316,dbnsfp35a,IARC,icgc21 -operation gx,r,f,f,f,f,f,f,f,f -nastring . -polish -vcfinput -xref ${humandb_annovar}/hg19_refGene.txt

    python ${python_vcf_to_csv} SimplifyVCF -inVCF ${sampleId}_${method}_annotated.hg${reference_version}_multianno.vcf -toType table -out simple_annotation_${sampleId}_${method}.csv

    python ${python_annot} -d . -f simple_annotation_${sampleId}_${method}.csv -o simple_annotation_${sampleId}_${method}.csv -m ${method}
    """
}


process VEP {
/*
INPUT
    (linked)
< sampleId                      : patient ID value
    < method                    : method used to annotate the .bam file
                                                           (Varscan, Mutect2, Hapotypecaller (GATK) or Pindel)
    < ${sampleId}_${method}.vcf : path to the vcf corresponding to the sampleId and
                                  method used

OUTPUT
    (linked)
    > sampleId                  : patient ID value
    > method                    : method used to annotate the .bam file
    > ${sampleId}_${method}.csv : comma separated values (.csv) file containing annotated variants

INPUT FROM
<- 

OUTPUT INTO
-> MERGE_VEP_FILES
*/
    cpus 4
    tag "${sampleId}"
    publishDir "${params.results}/${sampleId}/Variation", mode: 'copy'

    input:
    tuple val(sampleId), path("${sampleId}_merged_dedup.vcf")
    path(reference_genome)
    path(CADD_vep_plugin_indel) 
    path(CADD_vep_tabix_indel) 
    path(CADD_vep_plugin_snv) 
    path(CADD_vep_tabix_snv) 
    path(dbNSFP_vep_plugin) 
    path(dbNSFP_vep_tabix)
    path(python_vep_process)

    output:
    tuple val(sampleId), path("${sampleId}_annotated.*")

    script:
    """
    vep -i ${sampleId}_merged_dedup.vcf \
        --fasta ${reference_genome} \
        --offline \
        --dir_cache ${params.vep_cache} \
        --dir_plugins ${params.vep_plugins} \
        --fork $task.cpus \
        --pick_allele \
        --force_overwrite \
        --no_stats \
        --exclude_predicted \
        --sift b \
        --polyphen b \
        --ccds \
        --symbol \
        --numbers \
        --domains \
        --regulatory \
        --canonical \
        --protein \
        --biotype \
        --af \
        --max_af \
        --pubmed \
        --uniprot \
        --mane \
        --variant_class \
        --show_ref_allele \
        --keep_csq \
        --plugin CADD,${CADD_vep_plugin_indel},${CADD_vep_plugin_snv} \
        --plugin dbNSFP,${dbNSFP_vep_plugin},"CADD_phred,ClinPred_pred,ClinPred_rankscore,ExAC_AF,FATHMM_pred,HGVSc_VEP,HGVSp_VEP,MutPred_Top5features,MutPred_rankscore,REVEL_rankscore,UK10K_AF,clinvar_MedGen_id,clinvar_OMIM_id" \
        --fields "IMPACT,CLIN_SIG,Existing_variation,SYMBOL,MANE_SELECT,CDS_position,Codons,Protein_position,Amino_acids,VARIANT_CLASS,BIOTYPE,Feature,CANONICAL,EXON,INTRON,DOMAINS,miRNA,AF,ExAC_AF,MAX_AF,UK10K_AF,SIFT,PolyPhen,CADD_phred,ClinPred_pred,ClinPred_rankscore,FATHMM_pred,MutPred_Top5features,MutPred_rankscore,REVEL_rankscore,clinvar_MedGen_id,clinvar_OMIM_id" \
        --vcf \
        -o ${sampleId}_annotated.vcf
    
    python3.6 ${python_vep_process} ${sampleId}_annotated.vcf
    """
// TODO : expand with plugins
}


process MERGE_ANNOVAR_FILES {
/*
Gets all 4 annotated and filtered variation files and merges them to get the final variation file

INPUT
    < python_annot                                       : path to the python script merging variants csv 
                                                           into a single .csv file

    (linked)
    < sampleId                                           : patient ID value
    < method                                             : method used to annotate the .bam file 
                                                           (Varscan, Mutect2, Hapotypecaller (GATK) or Pindel)
    < Filter_simple_annotation_${sampleId}_${method}.csv : comma separated values (.csv) file containing variants
                                                           annotated and filtered

OUTPUT
    (linked)
    > sampleId                : patient ID value
    > variants_${sampleId}.csv: comma separated values (.csv) file containing variants
                                annotated and filtered, with merged values from all variant callers

INPUT FROM
    <- python_annot                                                        : python_annot value channel
    <- sampleId, method, Filter_simple_annotation_${sampleId}_${method}.csv: ANNOVAR_VSC, ANNOVAR_MT2, 
                                                                             ANNOVAR_HTC, ANNOVAR_PDL processes
    <- annotation_dict                                                     : python_dict value channel
                                                                             OR
                                                                             ANNOTATION_DICTIONNARY_SETUP process

OUTPUT INTO
    -> MERGE_ANNOVAR_FILES

// Possible improvement: annotation of a variable number of variant callers
*/
    tag "${sampleId}"
    publishDir "${params.results}/${sampleId}/Variation", mode: 'copy', pattern: 'Annotation_*.csv'
    publishDir "${params.results}/${sampleId}", mode: 'copy', pattern: 'Annotation_*.xlsx'

    input:
        path(python_annot)
        tuple val(sampleId), val("varscan"), path("Filter_simple_annotation_${sampleId}_varscan.csv"), val("mutect2"), path("Filter_simple_annotation_${sampleId}_mutect2.csv"), val("gatk"), path("Filter_simple_annotation_${sampleId}_gatk.csv"), val("pindel"), path("Filter_simple_annotation_${sampleId}_pindel.csv")
        path(annotation_dict)
        path(stats_dict)

    output:
        path("*")

    script:
    """
    python ${python_annot} -d . -f Filter_simple_annotation_${sampleId}_varscan.csv,Filter_simple_annotation_${sampleId}_mutect2.csv,Filter_simple_annotation_${sampleId}_gatk.csv,Filter_simple_annotation_${sampleId}_pindel.csv -o variants_${sampleId} -i ${annotation_dict} -r ${params.run_id} -m merge
    python ${python_annot} -d . -f Final_variants_${sampleId}.csv -o Annotation_${sampleId} -i ${annotation_dict} -m statistics -s ${stats_dict}
    mv Annotation_${sampleId} ./Annotation_${sampleId}.csv
    """
}
