mkdir $1/fastq
mkdir $1/work
mkdir $1/results $1/results/reports
nextflow run /home/bioinfo/FoxNGS/FoxNGS.nf -profile normal -w $1/work --reads $1/fastq --results $1/results -with-report -with-timeline
mv timeline.html $1/results/reports
mv report.html $1/results/reports
multiqc $1/results -o $1/results/reports
nextflow clean -f;
