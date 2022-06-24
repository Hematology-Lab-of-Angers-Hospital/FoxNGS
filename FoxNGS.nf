/* 
Define the default parameters
*/ 
params.genome         = "$baseDir/data/reference/*.fna"
params.makeIndex      = false
params.index          = "$baseDir/data/reference/*.{amb,ann,bwt,pac,sa}"
params.bed_hemato     = "$baseDir/data/bed/SureSelect-HEMATO-v7_UBA1.bed"
params.bed_exon       = "$baseDir/data/bed/DESIGN-FH-EXONS-gene_panel.bed"
params.bed_pindel     = "$baseDir/data/bed/Pindel_search_CALR-9_FLT3-14-15.bed"
params.picard         = "$baseDir/tools/picard.jar"
params.gatk           = "$baseDir/tools/gatk-package-4.2.6.1-local.jar"
params.varscan        = "$baseDir/tools/VarScan.v2.3.9.jar"
params.python_annot   = "$baseDir/tools/database_dictionnary.py"
params.python_dict    = "$baseDir/tools/annotation_process.py"  
params.python_iarc    = "$baseDir/tools/dvt_database_IARC.py"
params.reads          = "$baseDir/inputs/reads"
params.results        = "$baseDir/results"
params.report_rscript = "$baseDir/tools/R_Quality_SMPHD.Rmd"

log.info """\
F O X - N G S   V 1.0a
================================
genome     : $params.genome
make index : $params.makeIndex
index      : $params.index
bed_hemato : $params.bed_hemato
bed_exon   : $params.bed_exon
picard     : $params.picard
gatk       : $params.gatk
reads      : $params.reads
results    : $params.results
"""

include {
    SAM_SETUP;
    BAM_SETUP;
    BAM_MAPPING;
    ON_TARGET_MAPPING;
    DUPMARK_BAM_SETUP;
    COVERAGE_ANALYSIS;
    VARSCAN;
    MUTECT2;
    setup_paired_end_reads;
    setup_index;
} from './modules.nf'

workflow {
    fastq_input = Channel.fromFilePairs("${params.reads}/*_R{1,2}_*.fastq.gz")
    reference_genome = Channel.fromPath(params.genome)
    reference_index = Channel.fromPath(params.index).buffer(size:5)
    bed_hemato = Channel.fromPath(params.bed_hemato)
    bed_exon = Channel.fromPath(params.bed_exon)
    report_rscript = Channel.fromPath(params.report_rscript)
    picard = Channel.fromPath(params.picard)
    gatk = Channel.fromPath(params.gatk)
    varscan = Channel.fromPath(params.varscan)

    if ( params.makeIndex ) {
        INDEX_SETUP(reference_genome)
        SAM_SETUP(
            setup_paired_end_reads(fastq_input),
            reference_genome,
            INDEX_SETUP.out )
    }

    else {
        SAM_SETUP(
            setup_paired_end_reads(fastq_input),
            reference_genome,
            setup_index(reference_index))
    }

    BAM_SETUP(SAM_SETUP.out)
    BAM_MAPPING(BAM_SETUP.out.sorted_bam)
    ON_TARGET_MAPPING(BAM_MAPPING.out, bed_hemato)
    DUPMARK_BAM_SETUP(ON_TARGET_MAPPING.out, picard)
    COVERAGE_ANALYSIS(BAM_MAPPING.out, ON_TARGET_MAPPING.out, DUPMARK_BAM_SETUP.out, bed_hemato, bed_exon) // report_rscript
    VARSCAN(DUPMARK_BAM_SETUP.out, varscan, reference_genome)
    MUTECT2(DUPMARK_BAM_SETUP.out, reference_genome, gatk)
//    HAPLOTYPECALLER(DUPMARK_BAM_SETUP.out, reference_genome, gatk)
}