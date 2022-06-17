/* 
Define the default parameters
*/ 
params.genome     = "$baseDir/data/reference/*.fna"
params.makeIndex  = false
params.index      = "$baseDir/data/reference/*.{amb,ann,bwt,pac,sa}"
params.bed_hemato = "$baseDir/data/reference/"
params.bed_exon   = "$baseDir/data/reference/"
params.reads      = "$baseDir/inputs/reads/*_R{1,2}_*.fastq.gz"
params.results    = "$baseDir/results"

log.info """\
F O X - N G S   V 1.0a
================================
genome     : $params.genome
make index : $params.makeIndex
index      : $params.index
reads      : $params.reads
results    : $params.results
"""

include {
    SAM_SETUP;
    BAM_SETUP;
    BAM_MAPPING;
    setup_reads_channel;
    setup_index_channel;
} from './modules.nf'

workflow {
    fastq_input = Channel.fromFilePairs(params.reads)
    reference_genome = Channel.fromPath(params.genome)
    reference_index = Channel.fromPath(params.index).buffer(size:5)
    bed_hemato = Channel.fromPath(params.bed_hemato)
    bed_exon = Channel.fromPath(params.bed_exon)

    if ( params.makeIndex ) {
        INDEX_SETUP(reference_genome)
        SAM_SETUP(
            setup_reads_channel(fastq_input),
            reference_genome,
            INDEX_SETUP.out )
    }

    else {
        SAM_SETUP(
            setup_reads_channel(fastq_input),
            reference_genome,
            setup_index_channel(reference_index))
    }

    BAM_SETUP(SAM_SETUP.out)
    BAM_MAPPING(BAM_SETUP.out.sorted_bam)
    ON_TARGET_MAPPING(BAM_MAPPING.out, bed_hemato)
    
}