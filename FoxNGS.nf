/// ----- F O X   N G S   M A I N -----

ANSI_GREEN = "\033[1;32m"
ANSI_RED = "\033[1;31m"
ANSI_ORANGE = "\033[48;2;255;165;0m"
ANSI_RESET = "\033[0m"

log.info """\
${ANSI_ORANGE}F O X - N G S   V $params.pipeline_version
===============================================================================
${ANSI_RESET}

run number                 : $params.run_id
make index                 : $params.make_index
make fai                   : $params.make_fai
make reference dictionary  : $params.make_ref_dict
make annotation dictionary : $params.make_annot_dict
reference genome           : $params.genome
reference_version          : $params.reference_version
index files location       : $params.index
bed_bait                   : $params.bed_bait
bed_exon                   : $params.bed_exon
bed_pindel                 : $params.bed_pindel
reads location             : $params.reads
picard                     : $params.picard
gatk                       : $params.gatk
varscan                    : $params.varscan
results location           : $params.results
help                       : $params.help
"""


def helpMessage() {
    log.info"""${ANSI_GREEN}
        ${ANSI_ORANGE}Fox NGS - Version ${params.pipeline_version}${ANSI_GREEN}
        
        Use: ${ANSI_RED}nextflow run FoxNFS.nf [-- parameters] -profile [profile] [meta]${ANSI_GREEN}
        
        Parameters
            Can be accessed using the ${ANSI_RED}--[parameter_name]${ANSI_GREEN} in the command line
            See nextflow.config for parameter list and default values
        
        Profiles (On local executor)
            ${ANSI_RED}normal${ANSI_GREEN} : 8 CPU
            ${ANSI_RED}fast${ANSI_GREEN} : 16 CPU
            ${ANSI_RED}urgent${ANSI_GREEN} : 46 CPU
        
        pipeline insight
            ${ANSI_RED}-with-report${ANSI_GREEN} for a complete report on ressource attribution/usage
            ${ANSI_RED}-with-timeline${ANSI_GREEN} for a complete timeline of process execution
            ${ANSI_RED}-with-dag [file.name]${ANSI_GREEN} for a workflow diagram

        pipeline restart
            If for some reason the pipeline has failes and you need to restart it
            Use ${ANSI_RED}-resume${ANSI_GREEN} in the command to call cached processes
    ${ANSI_RESET}""".stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// PIPELINE MODULES CALLING
include {
    FAI_SETUP;
    REFERENCE_DICT_SETUP;
    REFERENCE_INDEX_SETUP;
    VARIATION_INDEX_SETUP;
    BED_INTERVAL_SETUP as BAIT_BED_INTERVAL_SET;
    BED_INTERVAL_SETUP as EXON_BED_INTERVAL_SET;
    ANNOTATION_DICTIONNARY_SETUP;
    QUALITY_CONTROL;
    BAM_SETUP;
    BAM_MAPPING;
    ON_TARGET_MAPPING;
    DUPMARK_BAM_SETUP;
    COLLECT_HS_METRICS;
    COVERAGE_ANALYSIS;
    PATIENT_REPORT;
    VARSCAN;
    MUTECT2;
    HAPLOTYPECALLER;
} from './modules.nf'

// WORKFLOW FUNCTIONS AND PROCESSES CALLING
workflow {
    // FILES CHANNEL

    Channel.fromFilePairs("${params.reads}/*_S*_R{1,2}_001.fastq.gz")
        .map{ left, right ->  
            def sampleId = left.split('_')[0]
            def fastq_r1 = right[0]
            def fastq_r2 = right[1]
            tuple(sampleId, fastq_r1, fastq_r2)
        }
        .set { fastq_input }

    // VARIABLES SETUP
    reference_version = params.reference_version

    // FILES SETUP
    // Reference genome and derivatives (indexes, fasta index and dictionary)
    reference_genome      = file("${params.genome}/hg${params.reference_version}/*.fna")
    reference_index       = file("${params.index}/*.{ann,amb,bwt,pac,sa}")
    reference_fai         = file("${params.genome}/hg${params.reference_version}/*.fna.fai")
    reference_dict        = file("${params.genome}/hg${params.reference_version}/*.dict")
    // bed files
    bed_bait              = file(params.bed_bait)
    bed_exon              = file(params.bed_exon)
    bed_hotspots          = file(params.bed_hotspots)
    //Coverage report templates
    exon_template         = file(params.exon_template)
    hotspot_template      = file(params.hotspot_template)
    patient_report_config = file(params.patient_report_config)
    // Software files
    picard                = file(params.picard)
    gatk                  = file(params.gatk)
    varscan               = file(params.varscan)

    // DATA PROCESSING
    // Quality control
    QUALITY_CONTROL(fastq_input)

    // Indexes setup
    REFERENCE_DICT_SETUP(gatk, reference_genome)
    FAI_SETUP(reference_genome)
    BAIT_BED_INTERVAL_SET(picard, bed_bait, "bait", REFERENCE_DICT_SETUP.out)
    EXON_BED_INTERVAL_SET(picard, bed_exon, "exon", REFERENCE_DICT_SETUP.out)

    // SAM/BAM setup + Processing
    if ( params.make_index ) {
        REFERENCE_INDEX_SETUP(reference_genome)
        BAM_SETUP(
            fastq_input,
            reference_genome,
            REFERENCE_INDEX_SETUP.out)
    }

    else {
        BAM_SETUP(
            fastq_input,
            reference_genome,
            reference_index)
    }

    BAM_MAPPING(BAM_SETUP.out)
    ON_TARGET_MAPPING(BAM_MAPPING.out, bed_bait)
    DUPMARK_BAM_SETUP(ON_TARGET_MAPPING.out, picard)

    // Coverage analysis
    COLLECT_HS_METRICS(picard, BAM_MAPPING.out, reference_genome, reference_fai, BAIT_BED_INTERVAL_SET.out, EXON_BED_INTERVAL_SET.out)

    coverage_channel = BAM_MAPPING.out.combine(ON_TARGET_MAPPING.out, by: 0)

    COVERAGE_ANALYSIS(coverage_channel, bed_bait, bed_exon)

    patient_report_channel = HAPLOTYPECALLER.out.base_recalibration.combine(
        COVERAGE_ANALYSIS.out.mosdepth_stats.combine(
            COLLECT_HS_METRICS.out,
            by: 0),
        by: 0)

    PATIENT_REPORT(patient_report_channel, bed_hotspots, bed_exon, exon_template, hotspot_template, patient_report_config)

    // Variant calling
    VARSCAN(DUPMARK_BAM_SETUP.out.bam, reference_genome, varscan)
    MUTECT2(gatk, DUPMARK_BAM_SETUP.out.bam, reference_genome, FAI_SETUP.out, REFERENCE_DICT_SETUP.out)
    HAPLOTYPECALLER(gatk, DUPMARK_BAM_SETUP.out.bam, reference_genome, FAI_SETUP.out, REFERENCE_DICT_SETUP.out, dbsnp, VARIATION_INDEX_SETUP.out)

    // Variant annotation
}