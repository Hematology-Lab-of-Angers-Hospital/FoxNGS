/// ----- F O X   N G S   M A I N -----

ANSI_GREEN = "\033[1;32m"
ANSI_RED = "\033[1;31m"
ANSI_ORANGE = "\033[48;2;255;165;0m"
ANSI_RESET = "\033[0m"

log.info """\
${ANSI_ORANGE}F O X - N G S   V $params.pipeline_version
===============================================================================
${ANSI_RESET}

run number                 : $params.run_number
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
entrypoint                 : $params.entrypoint
bcl_runfolder              : $params.bcl
samplesheet                : $params.samplesheet
bcl2fastq_image            : $params.bcl2fastq_image
reads location             : $params.reads
picard                     : $params.picard
gatk                       : $params.gatk
report R script            : $params.report_rscript
varscan                    : $params.varscan
pindel                     : $params.pindel
annovar                    : $params.annovar
human annovar database     : $params.humandb_annovar
vcf to csv python script   : $params.python_vcf_to_csv
annotation python script   : $params.python_annot
dict python script         : $params.python_dict   
iarc python script         : $params.python_iarc
transcripts base           : $params.base_transcript
artifacts base             : $params.base_artifact
variants_dictionary        : $params.annotation_dict
results location           : $params.results
help                       : $params.help
"""


def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
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
//    ANNOTATION_DICTIONNARY_SETUP;
    QUALITY_CONTROL;
    BAM_SETUP;
    BAM_MAPPING;
    ON_TARGET_MAPPING;
    DUPMARK_BAM_SETUP;
    COLLECT_HS_METRICS;
    COVERAGE_ANALYSIS;
    VARSCAN;
    MUTECT2;
    HAPLOTYPECALLER;
//    PINDEL;
    ANNOVAR;
    ANNOVAR as ANNOVAR_VSC;
    ANNOVAR as ANNOVAR_MT2;
    ANNOVAR as ANNOVAR_HTC;
//    ANNOVAR as ANNOVAR_PDL;
//    MERGE_ANNOTATION_FILES;
} from './modules.nf'

// WORKFLOW FUNCTIONS AND PROCESSES CALLING
workflow {
    // FILES CHANNEL
    // Raw files
    bcl_runfolder     = Channel.fromPath("${params.bcl}")

    if ( params.entrypoint == "fastq") {
        Channel.fromFilePairs("${params.reads}/*_R{1,2}_*.fastq.gz")
            .map{ left, right ->  
                def sampleId = left.split('_')[0]
                def fastq_r1 = right[0]
                def fastq_r2 = right[1]
                tuple(sampleId, fastq_r1, fastq_r2)
            }
            .set { fastq_input }
    }
    // VARIABLES SETUP
    reference_version = params.reference_version
    make_annot_dict = params.make_annot_dict

    // FILES SETUP
    // Sequence run samplesheet
    samplesheet       = file("${params.samplesheet}")
    // Reference genome and derivatives (indexes, fasta index and dictionary)
    reference_genome  = file("${params.genome}/hg${params.reference_version}/*.fna")
    reference_index   = file("${params.index}/*.{ann,amb,bwt,pac,sa}")
    reference_fai     = file("${params.genome}/hg${params.reference_version}/*.fna.fai")
    reference_dict    = file("${params.genome}/hg${params.reference_version}/*.dict")
    // annotation dictionary
    annotation_dict   = file("${params.annotation_dict}")
    // bed files
    bed_bait          = file(params.bed_bait)
    bed_exon          = file(params.bed_exon)
    bed_pindel        = file(params.bed_pindel)
    // Variation files
    dbsnp             = file(params.dbsnp)
    humandb_annovar   = file(params.humandb_annovar)
//    report_rscript    = Channel.fromPath(params.report_rscript)
    // Singularity images 
    bcl2fastq_image   = file("${params.bcl2fastq_image}")
    // Software files
    picard            = file(params.picard)
    gatk              = file(params.gatk)
    varscan           = file(params.varscan)
    pindel            = file(params.pindel)
    annovar           = file(params.annovar)
    // Python scripts file
    vcf_to_csv_python = file(params.python_vcf_to_csv)
    python_annot      = file(params.python_annot)
    iarc_python       = file(params.python_iarc)


    // DATA PROCESSING
    // Quality control
    QUALITY_CONTROL(fastq_input)

    // Indexes setup
    REFERENCE_DICT_SETUP(gatk, reference_genome)
    VARIATION_INDEX_SETUP(gatk, dbsnp)
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

    BAM_MAPPING(BAM_SETUP.out.sorted_bam)
    ON_TARGET_MAPPING(BAM_MAPPING.out, bed_bait)
    DUPMARK_BAM_SETUP(ON_TARGET_MAPPING.out, picard)

    // Coverage analysis
    COLLECT_HS_METRICS(picard, BAM_MAPPING.out, reference_genome, reference_fai, BAIT_BED_INTERVAL_SET.out, EXON_BED_INTERVAL_SET.out)
    COVERAGE_ANALYSIS(BAM_MAPPING.out, ON_TARGET_MAPPING.out, BAM_SETUP.out.mapping_stats, bed_bait, bed_exon) // report_rscript
    
    // Variant calling
    VARSCAN(DUPMARK_BAM_SETUP.out.bam, reference_genome, varscan)
    MUTECT2(gatk, DUPMARK_BAM_SETUP.out.bam, reference_genome, FAI_SETUP.out, REFERENCE_DICT_SETUP.out)
    HAPLOTYPECALLER(gatk, DUPMARK_BAM_SETUP.out.bam, reference_genome, FAI_SETUP.out, REFERENCE_DICT_SETUP.out, dbsnp, VARIATION_INDEX_SETUP.out)
//    PINDEL(DUPMARK_BAM_SETUP.out, BAM_MAPPING.out, reference_genome, pindel, bed_pindel)

    // Variant annotation
    ANNOVAR_VSC(annovar, reference_version, VARSCAN.out.varscan_variation, humandb_annovar, vcf_to_csv_python, python_annot)  
    ANNOVAR_MT2(annovar, reference_version, MUTECT2.out, humandb_annovar, vcf_to_csv_python, python_annot)
    ANNOVAR_HTC(annovar, reference_version, HAPLOTYPECALLER.out.gatk_variation, humandb_annovar, vcf_to_csv_python, python_annot)
//    ANNOVAR_PDL(annovar, reference_version, PINDEL.out, humandb_annovar, vcf_to_csv_python, python_annot)

/*
    // Variant formatting and filtering
    if ( make_annot_dict ) {
        python_dict     = file(params.python_dict)
        transcript_base = file(params.transcript_base)
        artifact_base   = file(params.artifact_base)
        
        MAKE_ANNOTATION_DICTIONNARY(python_dict, artifact_base, transcript_base, pipeline_version)
        MERGE_ANNOTATION_FILES(python_annot, ANNOVAR_VSC.out, ANNOVAR_MT2.out, HAPLOTYPECALLER.out, ANNOVAR_PDL.out, MAKE_ANNOTATION_DICTIONNARY.out)
    }

    else {
        MERGE_ANNOTATION_FILES(python_annot, ANNOVAR_VSC.out, ANNOVAR_MT2.out, HAPLOTYPECALLER.out, ANNOVAR_PDL.out)
    }
*/
}