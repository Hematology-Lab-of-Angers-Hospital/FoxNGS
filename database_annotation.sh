#!/bin/bash
# Creation of new database of annotation for pipeline_SMPHD_v2.7.sh
# Note all commande for creation of database for pipeline_SMPHD_v2.6.sh in annovar.sh dans Archive
# Faire pareil pour toutes les databases
# Date 03/07/20

# ~/Bureau/Recherche/Pipeline/./database_annotation.sh
# Header
ANNOVAR=~/BioNGSTools/annovar
PATH=/media/t-chu-027/DATAPART1/Database
cd $PATH

# Download Database
echo $ANNOVAR
# Rappel Commande Annotations
#$ANNOVAR/table_annovar.pl variant/${method}/${name}.${method}.vcf $ANNOVAR_DB -buildver hg19 -out $NAME_REP_ANNOVAR/${method}/annotation_simple_${name}.${method} -remove -protocol refGene,cytoBand,cosmic90,gnomad211_exome,clinvar_20190305,dbnsfp35a,IARC,icgc21 -operation gx,r,f,f,f,f,f,f -nastring . -thread 16 -polish -vcfinput -xref $ANNOVAR_DB/hg19_refGene.txt 
	
#************************
# Cytoband
echo -e "perl $ANNOVAR/annotate_variation.pl -buildver hg19 -downdb cytoBand humandb_annovar/"
perl $ANNOVAR/annotate_variation.pl -buildver hg19 -downdb cytoBand humandb_annovar/

#************************
# Cosmic
#************************
# ***********************
# Remplace by Cosmicv89 by v90 d
# For coding
# Rappel Erreur
#NOTICE: Finished reading 10068058 mutation ID from the VCF file CosmicCodingMuts.vcf
#Error: COSMIC MutantExport format error: column 17 or 12 should be 'Mutation ID' or 'ID_NCV'
#*****************************
# Solution
# Avec COSV parent
awk -v OFS="\t" -F "\t" '
gsub("GENOMIC_MUTATION_ID","Mutation ID");1' CosmicMutantExport_error.tsv > CosmicMutantExport.tsv

perl /home/t-chu-027/BioNGSTools/annovar/prepare_annovar_user.pl -dbtype cosmic Cosmic_v90/CosmicMutantExport.tsv -vcf Cosmic_v90/CosmicCodingMuts.vcf > humandb_annovar/hg19_cosmic90.txt

# Amelioration pour SMPHD 3.1
# cosmic 92
mv CosmicMutantExport.tsv CosmicMutantExport_error.tsv
# Transformation du nom de colonne nécessaire
awk -v OFS="\t" -F "\t" 'gsub("GENOMIC_MUTATION_ID","Mutation ID");1' CosmicMutantExport_error.tsv > CosmicMutantExport_2.tsv
sed '1d' mon_fichier.txt > CosmicMutantExport.tsv
# Création
perl /home/t-chu-027/BioNGSTools/annovar/prepare_annovar_user.pl -dbtype cosmic Cosmic_v92/CosmicMutantExport.tsv -vcf Cosmic_v92/CosmicCodingMuts.vcf > humandb_annovar/hg19_cosmic92.txt
# CNV 
# **************************************************

# CNV Search Contre Cosmic codingMut
awk -v OFS="\t" -F "\t" '
gsub("CNV_ID","Mutation ID");1' CosmicCompleteCNA_error.tsv > CosmicCompleteCNA.tsv

mv CosmicCompleteCNA.tsv CosmicCompleteCNA_e.tsv
awk -v OFS="\t" -F "\t" '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$1,$18,$19,$20}' CosmicCompleteCNA_e.tsv > CosmicCompleteCNA.tsv
perl /home/t-chu-027/BioNGSTools/annovar/prepare_annovar_user.pl -dbtype cosmic Cosmic_v90/CosmicCompleteCNA.tsv -vcf Cosmic_v90/CosmicCodingMuts.vcf > humandb_annovar/hg19_cosmic90_CNV.txt

# Error
#NOTICE: Finished reading 0 COSMIC records in DB file Cosmic_v90/CosmicCompleteCNA.tsv
#WARNING: 1008000 COSMIC ID from MutantExport file cannot be found in VCF file (this may be normal if the VCF file only contains coding or noncoding variants
# COSV whereis?

# For non coding 
# Preparation of Non coding
awk -v OFS="\t" -F "\t" '
gsub("GENOMIC_MUTATION_ID","Mutation ID");1' CosmicNCV_error.tsv > CosmicNCV.tsv

perl /home/t-chu-027/BioNGSTools/annovar/prepare_annovar_user.pl -dbtype cosmic Cosmic_v90/CosmicNCV.tsv -vcf Cosmic_v90/CosmicNonCodingVariants.vcf > humandb_annovar/hg19_cosmic90_Non_Coding.txt
# Same error 
# Error in cosmic Cosmic_v90/CosmicCompleteCNA.tsv 

#************************
# gnomad211_exome
perl /home/t-chu-027/BioNGSTools/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_exome humandb_annovar/
# clinvar_20190305
https://github.com/macarthur-lab/clinvar
perl /home/t-chu-027/BioNGSTools/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20200316 humandb_annovar/

# dbnsfp35a
perl /home/t-chu-027/BioNGSTools/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp35a humandb_annovar/
# IARC Database
# Position chromosomique in hg19 hg 38 or hg 37 
# Parcours IARC
python3 ~/Bureau/Recherche/Pipeline/SMPHD_v2.8/dvt_database_IARC.py

#icgc21
perl /home/t-chu-027/BioNGSTools/annovar/prepare_annovar_user.pl -buildver hg19 -downdb -webfrom annovar icgc21 humandb_annovar/
# ***************************************
# Other
# mcap for non synonimous
perl /home/t-chu-027/BioNGSTools/annovar/prepare_annovar_user.pl -buildver hg19 -downdb -webfrom annovar mcap13 humandb_annovar/


# Clinvar 20190305 to 20190722
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20190722.vcf.gz non fourni
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20190722.vcf.gz.tbi non fourni
# https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/#download
~BioNGSTools/vt/vt decompose clinvar_20190722.vcf.gz -o temp.split.vcf
perl /home/t-chu-027/BioNGSTools/annovar/prepare_annovar_user.pl   -dbtype clinvar_preprocess2 temp.split.vcf -out temp.split2.vcf
#critique
~/BioNGSTools/vt/vt normalize temp.split2.vcf -r /home/t-chu-027/BioNGSTools/Genome/hg19/hg19.fa -o temp.norm.vcf -w 2000000
# [E::faidx_adjust_position] The sequence "1" was not found
#[variant_manip.cpp:72 is_not_ref_consistent] failure to extract base from fasta file: 1:949695-949695
#FAQ: http://genome.sph.umich.edu/wiki/Vt#1._vt_cannot_retrieve_sequences_from_my_reference_sequence_file
#It is common to use reference files based on the UCSC browsers database and from the Genome Reference Consortium.
#For example, HG19 vs Grch37.  The key difference is that chromosome 1 is represented as chr1 and 1 respectively in the 
#FASTA files from these 2 sources.  Just use the appropriate FASTA file that was used to generate your VCF file originally.
# reference hg37
~BioNGSTools/vt/vt normalize temp.split2.vcf -r /home/t-chu-027/BioNGSTools/Genome/hg37-p13/reference/GCF_000001405.25_GRCh37.p13_genomic.fna -o temp.norm.vcf -w 2000000
#[E::faidx_adjust_position] The sequence "1" was not found
#[variant_manip.cpp:72 is_not_ref_consistent] failure to extract base from fasta file: 1:949695-949695
#FAQ: http://genome.sph.umich.edu/wiki/Vt#1._vt_cannot_retrieve_sequences_from_my_reference_sequence_file
# Solution utile aussi pour hg37.p13
sed -i -e "s/>chr/>/g" hg19_without_chr.fa 
~BioNGSTools/vt/vt normalize temp.split2.vcf -r /home/t-chu-027/BioNGSTools/Genome/hg19/hg19_without_chr.fa -o temp.norm.vcf -w 2000000
#[variant_manip.cpp:96 is_not_ref_consistent] reference bases not consistent: 1:984429-984512  TGCAGCTCAGGTGGGCGGGGAGGGGACGGGGCCGGGGCAGCTCAGGTGGGCGGGGAGGGGACGGGCGGGGGAGGGGGGGCCGGG(REF) vs TGCAGCTCAGGTgggcggggaggggacggggccggggcagctcaggtgggcggggaggggacgggcgggggagggggggccggg(FASTA)
#[normalize.cpp:209 normalize] Normalization not performed due to inconsistent reference sequences. (use -n or -m option to relax this)
# UPPER lettre hg19.fa
awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' hg19_without_chr.fa > hg19_clinvar.fa
# Abandon
