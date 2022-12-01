fox_run () {
    mkdir ../fastq
    mkdir ../results ../results/reports
    singularity exec -B /media/tfli-0070/DATASAVE1 /home/tfli-0070/BioNGSTools/Singularity_database/bcl2fastq_docker.sif bcl2fastq --runfolder . --output-dir ../fastq --sample-sheet $1 --use-bases-mask Y150,I8,I8,Y150 --no-lane-splitting -r 32 -p 32 -w 32
    nextflow run /home/tfli-0070/Bureau/Data/Recherche/Pipeline/SMPHD_Routine/FoxNGS.nf -profile fast -w /media/tfli-0070/DATASAVE2/work --reads ../fastq --results ../results -with-report -with-timeline
    mv timeline.html ../results/reports
    mv report.html ../results/reports
    multiqc ../results -o ../results/reports
    nextflow clean -f;
    awk -i inplace '{gsub(/\\n/,"\n")}1' ../results/Statistic_couverture.csv; sed -i /^$/d ../results/Statistic_couverture.csv
}

fox_run $1