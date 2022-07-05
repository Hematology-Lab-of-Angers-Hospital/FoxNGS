/// ----- F O X   N G S   M A I N -----

ANSI_GREEN = "\033[1;32m"
ANSI_RED = "\033[1;31m"
ANSI_RESET = "\033[0m"

log.info """\
${ANSI_GREEN}F O X - N G S   V 1.0a
===============================================================================
${ANSI_RESET}

bcl_runfolder           : $params.bcl
samplesheet             : $params.samplesheet
bcl2fastq_image         : $params.bcl2fastq_image
genome                  : $params.genome
reads                   : $params.reads
make index              : $params.make_index
make fai                : $params.make_fai
make dictionnary        : $params.make_dictionnary
reference_version       : $params.reference_version
index files location    : $params.index
bed_hemato              : $params.bed_hemato
bed_exon                : $params.bed_exon
bed_pindel              : $params.bed_pindel
picard                  : $params.picard
gatk                    : $params.gatk
varscan                 : $params.varscan     
annotation python script: $params.python_annot
dict python script      : $params.python_dict   
iarc python script      : $params.python_iarc 
results                 : $params.results
"""


def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info"""
        JeeJ.
    """.stripIndent()
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
    INDEX_SETUP;
    SAM_SETUP;
    BAM_SETUP;
    BAM_MAPPING;
    ON_TARGET_MAPPING;
    DUPMARK_BAM_SETUP;
    COVERAGE_ANALYSIS;
//    VARSCAN;
    MUTECT2;
//    HAPLOTYPECALLER;
} from './modules.nf'

// WORKFLOW FUNCTIONS AND PROCESSES CALLING
workflow {
    // Channel setup
    bcl_runfolder    = Channel.fromPath("${params.bcl}")
    samplesheet      = file("${params.samplesheet}")
    bcl2fastq_image  = file("${params.bcl2fastq_image}")

    reference_genome = file("${params.genome}/hg${params.reference_version}/*.fna")

    reference_index  = file("${params.index}/*.{ann,amb,bwt,pac,sa}")

    reference_fai    = file("${params.genome}/hg${params.reference_version}/*.fna.fai")
    reference_dict   = file("${params.genome}/hg${params.reference_version}/*.dict")
//    dbsnp            = Channel.fromPath(params.dbsnp)
    bed_hemato       = file(params.bed_hemato)
    bed_exon         = file(params.bed_exon)
//    report_rscript   = Channel.fromPath(params.report_rscript)
    picard           = file(params.picard)
    gatk             = file(params.gatk)
    varscan          = file(params.varscan)

    Channel.fromFilePairs("${params.reads}/*_R{1,2}_*.fastq.gz")
        .map{ left, right ->  
            def sampleId = left
            def fastq_r1 = right[0]
            def fastq_r2 = right[1]
            tuple(sampleId, fastq_r1, fastq_r2)
        }
        .set { fastq_input }

    // SAM/BAM setup + Processing
    if ( params.make_index ) {
        INDEX_SETUP(reference_genome)
        SAM_SETUP(
            fastq_input,
            reference_genome,
            INDEX_SETUP.out)
    }

    else {
        SAM_SETUP(
            fastq_input,
            reference_genome,
            reference_index)
    }

    REFERENCE_DICT_SETUP(gatk, reference_genome)
    FAI_SETUP(reference_genome)
    BAM_SETUP(SAM_SETUP.out)
    BAM_MAPPING(BAM_SETUP.out.sorted_bam)
    ON_TARGET_MAPPING(BAM_MAPPING.out, bed_hemato)
    DUPMARK_BAM_SETUP(ON_TARGET_MAPPING.out, picard)

    // Coverage analysis
    COVERAGE_ANALYSIS(BAM_MAPPING.out, ON_TARGET_MAPPING.out, DUPMARK_BAM_SETUP.out, BAM_SETUP.out.mapping_stats, bed_hemato, bed_exon) // report_rscript
    
    // Variant calling
//    VARSCAN(DUPMARK_BAM_SETUP.out, varscan, reference_genome)
    MUTECT2(DUPMARK_BAM_SETUP.out, reference_genome, FAI_SETUP.out, REFERENCE_DICT_SETUP.out, gatk)
//    HAPLOTYPECALLER(DUPMARK_BAM_SETUP.out, reference_genome, dbsnp, gatk)

    // Variant annotation'
}
