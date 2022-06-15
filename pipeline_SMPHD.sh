#!/bin/bash
# Original author : Maxime Renard
# Current author : Joris Argentin
# Bioinformatics
# Laboratoire d'hématologie
# Variant research pipeline for MPN on a 70 genes panel

# Last Modified 26/01/21

# A prevoir: Include dbsnp annotation, improve dbsnfp41a


echo "#############################################################"
echo "#-----------------------------------------------------------#"
echo "#-                                                         -#"
echo "#-                	 FoxNGS Pipeline                   -#"
echo "#-                      SURESELECQTX Data                  -#"
echo "#-                     Version 4.0 Routine                 -#"
echo "#-                                                         -#"
echo "#-----------------------------------------------------------#"
echo "#############################################################"

# ******************************************************************************************
# FUNCTION
# ******************************************************************************************
# ************************************** PROGRAM CALL **************************************

# help
function HELP {
	echo -e "FoxNGS Pipeline"
	echo -e "Lancement du pipeline HD"
	echo -e "./pipeline_SMPHD.sh"
	echo -e "A subsequent menu will ask for parameters"
	echo -e "Version 4.0 Routine"
#	echo -e "Cette nouvelle version contient la création d'un dictionnaire"
#	echo -e "Remplacement de cosmic 92, dbSNP138  at ajout qualite"
#	echo -e " Si une modification de database ou d'appel de software a été modifié."
#	echo -e "Il est nécessaire de créer un autre dictionnaire"
	echo -e "Question qualité:"
	echo -e "Dans le répertoire /home/$USER/BioNGSTools/ version_Tools or Version_genome_REFERENCE"
	echo -e "FASTQC check (run quality assessment)"
	echo -e "QUALITY coverage check : Look at a html file"
	exit 1
}

# User-input manual pipeline termination
function STOP_PROGRAM () {
	echo "Do you want to stop FoxNGS ?"
	read reponse
	case $reponse in
		[nN]*) echo "FoxNGS will keep going.";;
		[yYoO]*) echo "$0 FoxNGS has been manually terminated."
			exit 0;;
		*) echo "Input error"
			exit 1;;
	esac
}

# Agent info logging
function LOG_AGENT () {  # classic logger
	local prefix="[Pipeline $0 $(date +%Y/%m/%d\ %H:%M:%S)] by $agent"
	echo "${prefix} $@" >&2
	echo "${prefix} $@" >> $LOG
}

# Reference genome selection
function GENOME_SELECTION () {
	echo -e "1 : hg19 (GRCh37.p13 - latest)\n"
	echo -e "2 : hg38 (GRCh37.p12)\n"
	echo -e "3 : DEFAULT - hg19.fa (legacy version)"
	echo -e "Select your reference genome : "
	echo "********"
}
# MAJ : Passer à GRCh38.p13

# OUTPUT FILES CREATION
function OUTPUT_FILE_CREATION () {
	DIRECTORY=${RUN}/${OUTPUT}
	mkdir $DIRECTORY

	LOG=$DIRECTORY/log_pipeline_SMPHD_Routine_v3.3.txt
	WORKFLOW=$DIRECTORY/Workfow.xml
	CONFIGURATION=$DIRECTORY/Configuration_file
	DESIGN=$DIRECTORY/Experimental_Design
	# Localisation file
	# Exit report
	QUALITY=/media/$USER/Elements/Result_NGS/$RUN_ID
	COUV=/media/$USER/Elements/Result_NGS/$RUN_ID
	mkdir $QUALITY
	mkdir $COUV
	# Recherche automatique des fichiers brutes
	# Recuperation des data ngs brutes
	foldernamebcl=$(ls | grep -vE "^[a-z]|^[A-Z]")
	# Chemin vers le fichier brute de séquençage bcl
	BCL=$RUN/$foldernamebcl
	# Récupération du fichier du RUN contenant les échantillons
	# Note : Ce fichier doit contenir le nom du run
	TPL=$BCL/SureSelect_QXT_$RUN_ID.csv
}
# Selecion des fichiers du génome de référence choisi par l'utilisateur
function REFERENCE_GENOME () {

    case $REFERENCE in
        hg19)
			echo "Vous travaillez avec le chromosome de référence GCF_000001405.25_GRCh37.p13 dernère modification le 19/04/17 " >> $LOG
			BWA_REF=/media/$USER/DATAPART2/Database/Genome/hg37-p13/reference/GCF_000001405.25_GRCh37.p13_ref
			BWA_FASTA=/media/$USER/DATAPART2/Database/Genome/hg38-p12/reference/GCF_000001405.25_GRCh37.p13_genomic.fna
			;;
		hg38 )
			echo -e "Vous travaillez avec le chromosome de référence GCF_000001405.38_GRCh38.p13 dernière modification le 25/06/2019 " >> $LOG
			BWA_REF=/media/$USER/DATAPART2/Database/Genome/hg38-p12/reference/GCF_000001405.39_GRCh38.p13_ref
			BWA_FASTA=/media/$USER/DATAPART2/Database/Genome/hg38-p12/reference/GCF_000001405.39_GRCh38.p13_genomic.fna
			;;
  		DEFAULT | *)
			echo -e "Vous travaillez avec le chromosome de référence hg19 datant de 2013. " >> $LOG
			BWA_REF=/media/$USER/DATAPART2/Database/Genome/hg19/hg19bwaidx
			BWA_FASTA=/media/$USER/DATAPART2/Database/Genome/hg19/hg19.fa
			;;
	esac
	echo -e "*********************************"
	echo -e "Note : Genome de référence"
	echo -e "Pour en savoir plus sur la préparation d'un génome regarder le script /home/$USER/Bureau/Recherche/Pipeline/Preparation_genome/Preparation_genome.sh"
	echo -e $BWA_REF >> $LOG
	echo -e $BWA_FASTA >> $LOG
	echo -e "See ~/BioNGSTools/Version_genome_reference pour plus d'information sur la version des génomes de référence"  >> $LOG
	echo -e "En ligne de commande : gedit ~/BioNGSTools/Version_genome_reference" >> $LOG
}

# Validation d'un fichier d'analyse
function VALIDATION () {
	echo -e "Rentrez OK ou NO"
	read validation
	case $validation in
		[yYoO]*) echo -e "Quality OK" >> $LOG;;
		[nN]*) echo -e "ERROR : Quality files are incorrect. " >> $LOG
			echo "$0 ERROR - Invalid analysis files"
			exit 0;;
		DEFAULT | *) echo "Input error (Default: OK)" >> $LOG
	esac
}

# patient_id reminder: changing current directory to patient directory
function PATIENT_ID_REMINDER () {
	name=$1
	echo -e "***********************************************\n"
	echo "Patient name: ${name}:"
	echo -e "***************************************\n" >> $LOG
	echo "Patient name: ${name}:" >> $LOG
	cd $DIRECTORY/$name
}



# *********************************************
#  Logiciels et Outils utilisés
# *********************************************
function DATABASE () {
	# See BioNGSTools/version_Tools pour plus d'information sur la version des outils
	echo -e "Fichier qualité des outils dans /home/$USER/Bureau/Recherche/Diagnostique_routine/Qualite_kalilab" >> $LOG
	# Tools
	BEDTOOLS=~/BioNGSTools/bedtools-2.27.0/bedtools2/bin/bedtools
	PICARD=~/BioNGSTools/picard/picard.jar

	# Analyse qualité - Script R
	RSCRIPT=~/Bureau/Recherche/Pipeline/SMPHD_Routine/R_Quality_SMPHD.Rmd
	# Caller
	VARSCAN=~/BioNGSTools/varscan/VarScan.v2.4.3.jar
	# Somatic Variant
	GATK=~/BioNGSTools/gatk-4.1.5.0/gatk-package-4.1.5.0-local.jar
	# Deletion
	PINDEL=~/BioNGSTools/pindel/./pindel
	PINDEL_VCF=~/BioNGSTools/pindel/./pindel2vcf

	# Préparation de l'Annotation
	VCFTOCSV=~/BioNGSTools/VCF-Simplify/VCF-Simplify.py
	# Annotation
	TRAIT_ANNOT=~/Bureau/Recherche/Pipeline/SMPHD_Routine/Treatment_of_Annotation.py
	EXIST_FILE=~/Bureau/Recherche/Script/exist_file.sh
	DICTIONNARY=~/Bureau/Recherche/Pipeline/SMPHD_Routine/database_dictionnary.py
	# ********************************************
	# Essai samtools 1.10
	SAMTOOLS_1_10=/home/$USER/BioNGSTools/samtools-1.10/./samtools

	# ********************************************
	# Database
	# Fichier
	BED=/media/$USER/DATAPART2/Database/Fichier_intersection_bed/Sure_Select_design/SureSelect-HEMATO-v7_UBA1.bed
	# Bed pour couverture et l'analyse qualité R
	BEDEXON=/media/$USER/DATAPART2/Database/Fichier_intersection_bed/Analyse_coverage/DESIGN-FH-EXONS-gene_panel.bed
	# Variant
	# Fichier Bed Pindel Ajout de NPM1
	BED_PINDEL=/media/$USER/DATAPART2/Database/Variant/Pindel_search_CALR-9_FLT3-14-15.bed
	DBSNP=/media/$USER/DATAPART2/Database/DB/dbsnp_138.hg19.vcf

	# Annotation
	ANNOVAR=~/BioNGSTools/annovar/
	ANNOVAR_DB=/media/$USER/DATAPART2/Database/humandb_annovar

	ANNOTATION_REP=/media/$USER/DATAPART2/Database/Annotation
	BASE_TRANSCRIT=$ANNOTATION_REP/Transcript_reference/Liste_genes_transcript_26_01_21.csv
	BASE_ARTEFACT=$ANNOTATION_REP/Artefact/Base_artefact_120220.csv
	# Dictionnaire annotation
	DICT_ANNOTATION=$ANNOTATION_REP/Database_annotation_26_01_21_v3.3.json
	# Exit Database Result
	STAT_DICT=/media/$USER/Elements/Result_NGS/Dictionnary_database_Routine/Statistic_Database_dictionnary_routine_v3.3.csv
	DICT=/media/$USER/Elements/Result_NGS/Dictionnary_database_Routine/Dictionnary_Database_variant_routine_v3.3.json

	# Methode d'appel de variant
	METHODE1="GATK"
	METHODE2="Mutect2"
	METHODE3="Varscan"
	METHODE4="Pindel"
}

# Menu des analyses des fichiers de NGS
function CHOIX_ANALYSE () {
	echo "***********************************"
	echo -e "1: Demultiplexing (fastq file setup)\n"
	echo -e "2: Quality_bam (fastq.gz -> bam) \n"
	echo -e "3: Variant_calling (4 variant calling) \n"
	echo -e "4: Annotation (Annotation des variants, Filtres et Ecriture du Fichier d'annotation)\n"
	echo -e "5: Analyse_bam (Lancement des étapes de variant calling et d'annotation (Pour les Tests))\n"
	echo -e "6: All (2ème étape : Lancement de l'analyse par patient \n du fichier brute à la création du rapport (Combinaison du choix 2 à 4))\n"
	echo "***********************************"
	echo "Input a choice (Demultiplexing)"
	echo "Which analysis ould you like to perform? "
	echo "********"
	read ANALYSE
	echo "********"
}
# Test d'entrée des paramètres si erreur
function TEST_PARA () {
	while [ -z $SELECTION ]
	do
		read ANALYSE
	done
}

# ****************
# Menu
function INTERFACE () {
	echo "****************************************************"
	echo "SURESELECT data Analysis Pipeline - Version 4.0 Routine"
	echo "****************************************************"
	echo "Input RUN data"
	echo "********"
	echo -e "Input user name"
	echo "********"
	read agent
	cd /media/$USER/DATAPART1/RUN_NGS/SURESELECTQXT/
	exp=$(ls)
	echo -e "Nom des RUN Sureselect :\n${exp}"
	echo "***********************************"
	echo "Saisir le nom RUN :"
	read RUN_ID
	# Choix du RUN à analyser
	echo "********"
	RUN=/media/$USER/DATAPART1/RUN_NGS/SURESELECTQXT/$RUN_ID
	cd $RUN
	fichier=$(ls | grep -E "^[a-z]|^[A-Z]")
	echo "***********************************"
	echo "Nommer le répertoire de sortie des résultats :"
	echo "Note tous les répertoires d'analyse de données doivent commencer obligatoirement par une lettre Minuscule ou Majuscule.\n
	Puisque le fichier de sortie du séquenceur commence par un chiffre (Moyen de le discerner)"
	echo "Attention : Le nommer avec la date de lancement d'analyse et le nom du lanceur est préférable"
	echo "Ne pas Mettre d'espace dans le nom du fichier mais séparer les charactères par des underscores"
	echo "Exemple Analyse_NGS_MR_Date"
	echo -e "Répertoire d'analyse déja existant: \n(Le répertoire ne s'écrasera pas si vous voulez continuer l'analyse du RUN dans le même répertoire)"
	echo -e $fichier
	echo "********"
	echo -e "Rentrer le nom de votre répertoire d'analyse :"
	echo "********"
	read OUTPUT
	# Recuperation de la localisation des fichiers d'interet
	OUTPUT_FILE_CREATION
	echo "********"
	echo "***********************************"
	# Reference genome selection
	# DEFAULT parameter  pour accelerer la procedure de lancement
	# Amélioration 4.0 Choix des genomes de référence (hg19,hg37 et hg38)
	REFERENCE="DEFAULT"
	# Information sur la localisation des génomes de référence
	REFERENCE_GENOME
	# Intermediary step check
	# None by default
	qualite="NO"
	CHOIX_ANALYSE
	# Analysis range selection
	if [ "$ANALYSIS" != "Demultiplexing" ]; then
		echo "***********************************"
		echo "Voulez vous lancer l'analyse sur tous les patients ou sur un patient?"
		echo -e "Choix 1 - all (Analyse demandé Sur tous les patients du RUN sélectioné\n"*
		echo -e "Choix 2 - selection (Une analyse sur 1 patient \nou une selection de plusieurs patients délimité par des,)\n Choix des patients dessous\n" ;
		echo -e "Rentrez le nom du choix (Exemple : all) "
		echo "********"
		read SELECTION

	fi
}

# ***********
# Rappel des paramètres
function PARAMETER_REMINDER () {
	echo -e "**********************************************************" >> $LOG
	LOG_AGENT
	echo -e "Input parameters recap" >> $LOG
	echo -e "Raw sequence file: ${BCL}" >> $LOG
	echo -e "NGS patient file: ${TPL}" >> $LOG
	echo -e "Run ID: ${RUN_ID}" >> $LOG
	echo -e "Reference genome: ${REFERENCE}" >> $LOG
	echo -e "***********************************">> $LOG
	echo -e "Output directory:">> $LOG
	echo -e $DIRECTORY >> $LOG
	echo -e "Planned analyses $ANALYSIS " >> $LOG
	if [ "$ANALYSIS" != "Demultiplexing" ]; then
		echo -e "Running on ${SELECTION}"
	fi
}

# Error handling
# checking if there's a file at input path
function VERIFY_FILE () {
	path_file=$1
	# Fonction qui verifie si le fichier existe
	file_exist=$($EXIST_FILE ${path_file})

	# Si le fichier n'exite pa ou il est vide : Creation du dictionnaire
	if [ $file_exist == "File_exist" ]
		then
		echo -e "Verification: File: ${file} exist"
	# S'il n'existe pas arret du programme
	else
		echo -e "Error File doesn't exist stop of analysis of patent"
		exit 1
	fi

}



# ******************************* ANALYSIS PIPELINE *****************************************
# *******************************************************************
# Data demultiplexing for each patient
# INPUT: bcl2fastq
# OUTPUT: fastq
# *******************************************************************
function FASTQ_SETUP () {
	cd $DIRECTORY
	echo -e "**********************************************************************\n" >> $LOG
	echo "Génération des fichiers FASTQ à partir des BCL" >> $LOG
	date >> $LOG
	echo "start de bcl2fastq" >> $LOG
	echo -e "bcl2fastq --runfolder $BCL --output-dir $DIRECTORY/ --sample-sheet $TPL --use-bases-mask Y150,I8,I8,Y150 --no-lane-splitting -r 16 -p 16 -w 16" >> $LOG

	echo "fin de bcl2fastq" >> $LOG
	date >> $LOG
	echo -e "**********************************************************************\n" >> $LOG
}

# ************************
# BAM file setup
# ************************

function QC () {
	# Getting name parameter
	name=$1

	# Fichier log de sortie
	PREPARATION_BAM_FILE=$DIRECTORY/$name/log_bam_test_Routine.txt

	# Création d'un répertoire pour chaque patient
	mkdir $DIRECTORY/$name
	# Création d'un répertoire de sortie de résultat dans le bureau
	mkdir $QUALITY/$name

	echo -e "**********************************************************************\n" > $PREPARATION_BAM_FILE
	date > $PREPARATION_BAM_FILE
	echo -e "Génération des Fichiers d'alignement BAM pour ${name}\n" >> $PREPARATION_BAM_FILE

	# Si relancement d'un patient redéplacement dans le dossier des résultats pour l'analyse
	mv $name/*.fastq.gz .

	# Récupération des fichiers fastq sens R1 et R2 correspondant à un identifiant
	R1=$(ls | grep $name | grep _R1_)
	R2=$(ls | grep $name | grep _R2_)

	echo "R1" >> $PREPARATION_BAM_FILE
	echo $R1 >> $PREPARATION_BAM_FILE
	echo "R2" >> $PREPARATION_BAM_FILE
	echo $R2 >> $PREPARATION_BAM_FILE

	# Déplacement des FASTQ dans le fichier du patient
	mv $R1 $DIRECTORY/$name
	mv $R2 $DIRECTORY/$name

	#On rentre dans le fichier du  patient
	echo -e "Lancement Quality bam" >> $LOG
	PATIENT_ID_REMINDER $name
	echo "Nom du patient analysé ${name}:" >> $PREPARATION_BAM_FILE
	echo $name >> $PREPARATION_BAM_FILE
	VERIFY_FILE $DIRECTORY/$name/$R1
	VERIFY_FILE $DIRECTORY/$name/$R2
	# *************************************************
	# Elaboration du FASTQC pour la patient :
	# *************************************************
	# Génération du fichier qualité
	echo -e "fastqc -o . $R1 $R2 -t 16" >> $PREPARATION_BAM_FILE
	fastqc -o . $R1 $R2 -t 16
	# Extension
	html="_fastqc.html"

	# copie des fichiers d'analyse fastqc vers le repertoire qualite version 3.1 Simplification cp
	cp $DIRECTORY/$name/*$html $QUALITY/$name/

	# Suppression des fichiers brutes en attente fichier intermediaire
	rm -dr $DIRECTORY/$name/*fastqc.zip
	#Creation d'un fichier temporaire: Stockage Fichier SAM à BAM de préparation
	mkdir $DIRECTORY/$name/tmp
	# ****************************************
	# Alignement contre le génome de référence
	# *****************************************
	echo -e "**********************************************************************" >> $PREPARATION_BAM_FILE
	echo -e "Alignement pour l'échantillon identifié ${name}" >> $PREPARATION_BAM_FILE
	#Construction du fichier SAM via BAW-MEM
	echo -e "Alignement contre le génome de référence:\n" >>$PREPARATION_BAM_FILE
	echo -e "Construction du fichier SAM via BAW-MEM" >>$PREPARATION_BAM_FILE
	date >> $PREPARATION_BAM_FILE
	echo -e "bwa mem -t 16 -R '@RG\tID:C5-${name}\tPL:illumina\tPU:HXXX\tLB:Solexa\tSM:C5-${name}' $BWA_REF $R1 $R2 -o tmp/${name}.sam" >> $PREPARATION_BAM_FILE
	bwa mem -t 16 -R "@RG\tID:C5-${name}\tPL:illumina\tPU:HXXX\tLB:Solexa\tSM:C5-${name}" $BWA_REF $R1 $R2 -o tmp/${name}.sam
	echo -e "Alignement effectué" >> $PREPARATION_BAM_FILE
	date >> $PREPARATION_BAM_FILE

	# Génération du fichier bam
	echo -e "Génération du fichier bam:\n" >> $PREPARATION_BAM_FILE
	echo -e "samtools view -@ 16 -Sh tmp/${name}.sam -bo tmp/${name}.bam" >> $PREPARATION_BAM_FILE
	samtools view -@ 16 -Sh tmp/${name}.sam -bo tmp/${name}.bam
	# Tri du fichier
	echo -e "samtools sort -@ 16 tmp/${name}.bam -o tmp/${name}.sort.bam">> $PREPARATION_BAM_FILE
	samtools sort -@ 16 tmp/${name}.bam -o tmp/${name}.sort.bam >> $PREPARATION_BAM_FILE
	echo -e "samtools index -@ 16 -b tmp/${name}.sort.bam" >> $PREPARATION_BAM_FILE
	samtools index -@ 16 -b tmp/${name}.sort.bam

	# Vérification des reads, ils sont bien mappés?
	# Total of read
	total=$(samtools view -h -c tmp/${name}.sort.bam )
	echo -e "*******************************" >> $PREPARATION_BAM_FILE
	echo -e "Total of read : $total" >> $PREPARATION_BAM_FILE
	echo -e "*******************************" >> $PREPARATION_BAM_FILE
	echo -e " Keep only mapped\n" >> $PREPARATION_BAM_FILE
	echo -e "samtools view -F 0x4 -@ 16 -h -b tmp/${name}.sort.bam >tmp/${name}.sort_mapped.bam" >> $PREPARATION_BAM_FILE
	samtools view -F 0x4 -h -@ 16 -b tmp/${name}.sort.bam >tmp/${name}.sort_mapped.bam


	echo -e "*******************************" >> $PREPARATION_BAM_FILE
	# Read mapped
	mapped=$(samtools view -h -c tmp/${name}.sort_mapped.bam)
	echo -e "Only of map : $mapped" >> $PREPARATION_BAM_FILE
	echo -e "*******************************" >> $PREPARATION_BAM_FILE
	# Unmapped and chimeric
	unmapped=$(samtools view -f 0x4 -h -@ 16 -c -b tmp/${name}.sort.bam)
	chimeric=$(samtools view  -f 0x800 -h -@ 16 -c -b tmp/${name}.sort.bam)
	echo -e "Only unmapped : $unmapped" >> $PREPARATION_BAM_FILE
	echo -e "Only chimeric : $chimeric" >> $PREPARATION_BAM_FILE

	echo -e "*******************************" >> $PREPARATION_BAM_FILE
	date >> $PREPARATION_BAM_FILE


	# Analyse qualité coverage All read mapped
	samtools index -@ 16 -b tmp/${name}.sort_mapped.bam
	echo -e "$SAMTOOLS_1_10 coverage tmp/${name}.sort_mapped.bam -m > tmp/${name}.sort_mapped_analyse_coverage.bed" >> $PREPARATION_BAM_FILE
	$SAMTOOLS_1_10 coverage tmp/${name}.sort_mapped.bam -m > tmp/${name}.sort_mapped_analyse_coverage.bed


	# Intersection avec le fichier du design
	# avec les région intronique aussi  pour déceler les mutations de type épissage en plus
	echo -e "$BEDTOOLS intersect -a tmp/${name}.sort_mapped.bam -b $BED  > tmp/${name}-on-target.bam" >> $PREPARATION_BAM_FILE
	$BEDTOOLS intersect -a tmp/${name}.sort_mapped.bam -b $BED  > tmp/${name}-on-target.bam

	# Analyse qualité coverage on target
	echo -e "samtools index  -@ 16 -b tmp/${name}-on-target.bam" >> $PREPARATION_BAM_FILE
	samtools index -@ 16 -b tmp/${name}-on-target.bam
	echo -e "$SAMTOOLS_1_10 coverage tmp/${name}-on-target.bam -m > tmp/${name}-on-target_analyse_coverage.bed" >> $PREPARATION_BAM_FILE
	$SAMTOOLS_1_10 coverage tmp/${name}-on-target.bam -m > tmp/${name}-on-target_analyse_coverage.bed

	echo -e "$BEDTOOLS intersect -a tmp/${name}.sort_mapped.bam -b $BED -v  > tmp/${name}-off-target.bam" >> $PREPARATION_BAM_FILE
	$BEDTOOLS intersect -a tmp/${name}.sort_mapped.bam -b $BED -v > tmp/${name}-off-target.bam

	# Analyse qualité coverage on target
	echo -e "samtools index  -@ 16 -b tmp/${name}-off-target.bam" >> $PREPARATION_BAM_FILE
	samtools index -@ 16 -b tmp/${name}-off-target.bam
	echo -e "$SAMTOOLS_1_10 coverage tmp/${name}-off-target.bam -m > tmp/${name}-off-target_analyse_coverage.bed" >> $PREPARATION_BAM_FILE
	$SAMTOOLS_1_10 coverage tmp/${name}-off-target.bam -m > tmp/${name}-off-target_analyse_coverage.bed



	# Verification off target
	totalmapnoffintersect=$(samtools view -h -@ 16 -c tmp/${name}-off-target.bam )
	echo -e "*******************************" >> $PREPARATION_BAM_FILE
	echo -e "Total of read mapped in panel gene: $totalmapnoffintersect" >> $PREPARATION_BAM_FILE

	# Mapped and intersected in region
	totalmapintersect=$(samtools view -h -@ 16 -c tmp/${name}-on-target.bam )
	echo -e "*******************************" >> $PREPARATION_BAM_FILE
	echo -e "Total of read mapped in panel gene: $totalmapintersect" >> $PREPARATION_BAM_FILE
	# Calcul du ratio on target
	ratio=$(echo "($totalmapintersect/$mapped)*100" | bc -l )
	# Calcul du ratio off target
	ratio_off=$(echo "($totalmapnoffintersect/$mapped)*100" | bc -l )
	# Conversion float to int to comparison
	int_ratio=${ratio%.*}
	int_ratio_off=${ratio_off%.*}
	echo -e "Ratio of read mapped in panel gene: $int_ratio" >> $PREPARATION_BAM_FILE
	echo -e "Ratio of read mapped out panel gene: $int_ratio_off" >> $PREPARATION_BAM_FILE
	limite_ratio=60
	echo -e "Seuil ratio of read mapped in panel gene: $limite_ratio" >> $PREPARATION_BAM_FILE

	# Condition
	# Librairie non correcte
	if [ "$int_ratio" -lt "$limite_ratio" ]
		then
		echo -e "Error of library of patent - Stop program : ${name}. ratio of map intersect is only ${ratio}." >> $PREPARATION_BAM_FILE
		echo -e "Error of library of patent - Stop program : ${name}."
	# Sinon librairie correcte
	else
		echo -e "Librairie correcte pour le patient : ${name}." >> $PREPARATION_BAM_FILE
	fi

	VERIFY_FILE	tmp/${name}-on-target.bam
	# Génération des statistiques: compte rendu qualité
	echo -e "Compte rendu qualité"

	echo -e "samtools flagstat -@ 16 tmp/${name}-on-target.bam > tmp/${name}.bam.sort.stat" >> $PREPARATION_BAM_FILE
	samtools flagstat -@ 16 tmp/${name}-on-target.bam > tmp/${name}.bam.sort.stat
	echo -e "Reindexation" >> $PREPARATION_BAM_FILE

	# Marquages des duplicates sans les enlever
	echo -e "Marquage des duplicates" >> $PREPARATION_BAM_FILE
	echo -e "java -jar $PICARD MarkDuplicates I= tmp/${name}-on-target.bam O=$name.sort.dupmark.bam M=$name.marked_dup.metrics.txt" >> $PREPARATION_BAM_FILE
	java -jar $PICARD MarkDuplicates I=tmp/${name}-on-target.bam O=${name}.sort.dupmark.bam M=$name.marked_dup.metrics.txt >> $PREPARATION_BAM_FILE
	# Reindexation
	samtools index -@ 16 -b ${name}.sort.dupmark.bam

	echo "************************************" >> $PREPARATION_BAM_FILE
	echo "Alignement pour l'échantillon $name terminé" >> $PREPARATION_BAM_FILE
	# ****************************************************
	# Couverture
	# ****************************************************
	# Pour analyser la couverture nous nous intéressons qu'aux exons pour la génération du fichier qualité : les introns sont enlevés
	echo -e "Analyse couverture" >> $PREPARATION_BAM_FILE
	echo -e "$BEDTOOLS coverage -a $BEDEXON -b $DIRECTORY/$name/${name}.sort.dupmark.bam -d >$DIRECTORY/$name/${name}_couverture_analyse.bed" >> $PREPARATION_BAM_FILE
	$BEDTOOLS coverage -a $BEDEXON -b $DIRECTORY/$name/${name}.sort.dupmark.bam -d >$DIRECTORY/$name/${name}_couverture_analyse.bed

	# ****************************************************
	# Analyse de la qualité de coverage Rmarkdown
	echo -e "R -e rmarkdown::render('${RSCRIPT}',
	params = list( directory='$(pwd)',file='${name}_couverture_analyse.bed',agent='${agent}',pipeline='$0',output='${COUV}/Statistic_couverture.csv',output_gene='/media/$USER/Elements/Result_NGS/Stat_gene/Statistic_couverture_gene.csv',ratio_library='${int_ratio}'),
		output_file='$(pwd)/${name}_couverture_analyse.bed.html')" >> $PREPARATION_BAM_FILE
	R -e "rmarkdown::render('${RSCRIPT}', params = list(directory='$(pwd)',file='${name}_couverture_analyse.bed',agent='${agent}',pipeline='${0}',output='${COUV}/Statistic_couverture.csv',output_gene='/media/$USER/Elements/Result_NGS/Stat_gene/Statistic_couverture_gene.csv',ratio_library='${int_ratio}'),output_file='$(pwd)/${name}_couverture_analyse.bed.html')"

	echo -e "**********************************************************************\n" >> $PREPARATION_BAM_FILE
	# Copie vers le répertoire qualité
	# Mise en commentaire car le code R ne génère pas fichier pour quelques patients
	#VERIFY_FILE $(pwd)/${name}_couverture_analyse.bed.html
	cp $(pwd)/${name}_couverture_analyse.bed.html $QUALITY/$name/
	# Si analyse du fichier qualité
	if [ $qualite = "OK" ]
		then
			echo -e "Test analyse qualité"
			echo -e "Afficher OK:  si le fichier de qualité est correcte\nNO: si incorrecte\n ********\n"
			firefox $(pwd)/${name}_couverture_analyse.bed.html
			VALIDATION
	fi
	date >> $PREPARATION_BAM_FILE

	# Suppresion des fichiers tmp
	# Suppression du fichier temporaire .sam, du bam ininial et du bam sort_mapped et du bam target
	rm -dr tmp/*sam tmp/${name}.sort_mapped.bam tmp/${name}.bam tmp/*2048*

	# Copie des fichiers BAM
	cp ${name}.sort.dupmark.bam ${name}.sort.dupmark.bam.bai tmp/*analyse_coverage.bed $QUALITY/$name/

}

# ***********************
# Variant calling
# ***********************

# ***
# Main function for variant calling
# ***
function VARIANT_CALLING () {

	# Récupération paramètre
	name=$1
	mkdir $QUALITY/$name
	# Fichier de sortie
	VARIANT_FILE=$DIRECTORY/$name/logvariantcalling_Routine.txt
	# function RAPPEL patient
	echo -e "Lancement Variant_calling" >> $LOG
	PATIENT_ID_REMINDER $name

	echo -e "**********************************************************************"  > $VARIANT_FILE
	echo -e "Variant calling" >> $VARIANT_FILE
	date >> $VARIANT_FILE
	VERIFY_FILE ${name}.sort.dupmark.bam
	# ***********************************************
	# Création du repertoire qui stocke les fichiers d'analyse variants
	mkdir variant
	# ******************************************
	# Detection par Varscan
	# ******************************************
	mkdir variant/${METHODE3}
	echo -e " ***************************" >> $VARIANT_FILE
	echo "Variant calling (Varscan 2)" >> $VARIANT_FILE
	echo "Start à :" >> $VARIANT_FILE
	date >> $VARIANT_FILE
	# Preparation varscan
	#  the following command lines call SNPs and short INDEL
	echo -e "samtools mpileup -Q 13 -q 0 -A -B -d 100000 -f $BWA_FASTA ${name}.sort.dupmark.bam > variant/${METHODE3}/${name}.${METHODE3}.mpileup" >> $VARIANT_FILE
	samtools mpileup -Q 13 -q 0 -A -B -d 100000 -f $BWA_FASTA ${name}.sort.dupmark.bam > variant/${METHODE3}/${name}.${METHODE3}.mpileup
	echo "java -jar $VARSCAN mpileup2cns variant/${METHODE3}/${name}.${METHODE3}.mpileup --min-coverage 50 --min-reads2 8 --min-avg-qual 30 --min-var-freq 0.02 --p-value 0.1 --strand-filter 0 --output-vcf --variants > variant/${METHODE3}/${name}.${METHODE3}.vcf" >> $VARIANT_FILE
	java -jar $VARSCAN mpileup2cns variant/${METHODE3}/${name}.${METHODE3}.mpileup --min-coverage 50 --min-reads2 8 --min-avg-qual 30 --min-var-freq 0.02 --p-value 0.1 --strand-filter 0 --output-vcf --variants > variant/${METHODE3}/${name}.${METHODE3}.vcf
	#VERIFY_FILE variant/${METHODE3}/${name}.${METHODE3}.vcf
	echo -e "Variant calling ${METHODE3} terminé\n" >> $VARIANT_FILE
	date >> $VARIANT_FILE
	# Supprimer File
	rm variant/${METHODE3}/${name}.${METHODE3}.mpileup
	# ******************************************
	# Détection par GATK
	# ******************************************
	mkdir variant/${METHODE1}
	echo -e " ***************************" >> $VARIANT_FILE
	echo "Variant calling (GATK)"
	echo "Variant calling (GATK)" >> $VARIANT_FILE
	echo "Start à :" >> $VARIANT_FILE
	date >> $VARIANT_FILE

	# Recalibration
	echo -e "Recalibration\n"  >> $VARIANT_FILE
	echo -e "java -jar $GATK BaseRecalibrator
			-I ${name}.sort.dupmark.bam
			-R $BWA_FASTA
			--known-sites $DBSNP
			-O variant/${METHODE1}/${name}_recal_data.table" >> $VARIANT_FILE
	java -jar $GATK BaseRecalibrator \
		-I ${name}.sort.dupmark.bam \
		-R $BWA_FASTA \
		--known-sites $DBSNP \
		-O variant/${METHODE1}/${name}_recal_data.table
	#Note for hg19 we use dbSNP138 because last release dbSNP 152 contain more error!!!
	echo -e "Apply BQSR\n"  >> $VARIANT_FILE
	echo -e "java -jar $GATK ApplyBQSR
				-I ${name}.sort.dupmark.bam
				-R $BWA_FASTA
				-bqsr variant/${METHODE1}/${name}_recal_data.table
				-O variant/${METHODE1}/${name}.bqsr.bam"	>> $VARIANT_FILE
	java -jar $GATK ApplyBQSR \
		-I ${name}.sort.dupmark.bam \
		-R $BWA_FASTA \
		-bqsr variant/${METHODE1}/${name}_recal_data.table \
		-O variant/${METHODE1}/${name}.bqsr.bam

	echo -e "Haplotype caller\n"  >> $VARIANT_FILE
	echo -e "java -jar $GATK HaplotypeCaller  \
			-R $BWA_FASTA \
			-I variant/${METHODE1}/$name.bqsr.bam \
			--native-pair-hmm-threads 16 \
			-O variant/${METHODE1}/${name}.${METHODE1}.vcf \
			--min-base-quality-score 30 \
			--minimum-mapping-quality 20 \
			--dont-use-soft-clipped-bases true" >> $VARIANT_FILE
	java -jar $GATK HaplotypeCaller  \
		-R $BWA_FASTA \
		-I variant/${METHODE1}/$name.bqsr.bam \
		--native-pair-hmm-threads 16 \
		-O variant/${METHODE1}/${name}.${METHODE1}.vcf \
		--min-base-quality-score 30 \
		--minimum-mapping-quality 20 \
		--dont-use-soft-clipped-bases true
	#remove date
	awk '{gsub(",Date=[^>]+>",">");}1' variant/${METHODE1}/$name.${METHODE1}.vcf
	VERIFY_FILE variant/${METHODE1}/${name}.${METHODE1}.vcf
	date >> $VARIANT_FILE
	rm variant/${METHODE1}/$name.bqsr.bam variant/${METHODE1}/$name.bqsr.bai
	# *************************************
	# MUTECT2
	# *************************************

	mkdir variant/${METHODE2}
	date >> $VARIANT_FILE
	# Suppression des anciens fichiers
	rm variant/${METHODE2}/${name}.${METHODE2}.vcf.gz variant/${METHODE2}/${name}.${METHODE2}.vcf
	echo -e " ***************************" >> $VARIANT_FILE
	echo -e "MUTECT2" >> $VARIANT_FILE
	echo -e "java -jar $GATK Mutect2 -R $BWA_FASTA -I ${name}.sort.dupmark.bam  --min-base-quality-score 30 --native-pair-hmm-threads 16 --dont-use-soft-clipped-bases true -O variant/${METHODE2}/${name}.${METHODE2}.vcf.gz" >> $VARIANT_FILE
	java -jar $GATK Mutect2 -R $BWA_FASTA -I ${name}.sort.dupmark.bam --min-base-quality-score 30 --dont-use-soft-clipped-bases true --native-pair-hmm-threads 16 -O variant/${METHODE2}/${name}.${METHODE2}.vcf.gz
	echo -e "gunzip variant/${METHODE2}/${name}.${METHODE2}.vcf.gz" >> $VARIANT_FILE
	gunzip variant/${METHODE2}/${name}.${METHODE2}.vcf.gz
	VERIFY_FILE variant/${METHODE2}/${name}.${METHODE2}.vcf
	# Test Filter Mutect Call
	#java -jar $GATK Mutect2 FilterMutectCalls -V variant/${METHODE2}/${name}.${METHODE2}.vcf -R $BWA_FASTA -O variant/${METHODE2}/${name}.${METHODE2}_filter.vcf

	# *************************************
	# PINDEL
	# *************************************
	mkdir variant/${METHODE4}
	# create fichier configfile_pindel
	# Par defaut
	# 150 + 150 + 300 ins/del + marge 200
	date >> $VARIANT_FILE
	echo -e "Pindel" >> $VARIANT_FILE
	echo -e "${name}.sort.dupmark.bam 800 Duplicate_mark\ntmp/${name}.sort.bam 800 All_read" >> $VARIANT_FILE
	echo -e "${name}.sort.dupmark.bam 800 Duplicate_mark\ntmp/${name}.sort.bam 800 All_read" > variant/${METHODE4}/config_file_pindel.txt
	# ITD300
	echo -e "$PINDEL -f $BWA_FASTA -i $DIRECTORY/$name/variant/${METHODE4}/config_file_pindel.txt -j $BED_PINDEL -T 14 -o variant/${METHODE4}/ITD_${name}" >> $VARIANT_FILE
	$PINDEL -f $BWA_FASTA -i $DIRECTORY/$name/variant/${METHODE4}/config_file_pindel.txt -j $BED_PINDEL -T 14 -o variant/${METHODE4}/ITD_${name}

	#D: Deletion + TD:Tandem Duplication + INV: Inversion + INS short and long insertion
	echo -e "$PINDEL_VCF -p variant/${METHODE4}/ITD_${name}_D -r $BWA_FASTA  -R x -d 00000000 -G  -v variant/${METHODE4}/${name}.${METHODE4}.vcf" >> $VARIANT_FILE
	$PINDEL_VCF -P variant/${METHODE4}/ITD_${name} -r $BWA_FASTA  -R x -d 00000000 -G -v variant/${METHODE4}/${name}.${METHODE4}.vcf
	VERIFY_FILE variant/${METHODE4}/${name}.${METHODE4}.vcf
	date >> $VARIANT_FILE
	echo -e "**********************************************************************"  >> $VARIANT_FILE
}


# ******************
# Annotation
# ******************

# VCF to csv conversion
# Input: VCF file
# Output:
function VCFTOTABLE (){
	name=$1
	method=$2
	mkdir $NAME_REP_ANNOVAR/Table/

	echo -e "*********\nVCF to Table\n*********"
	date >> $ANNOTATION_FILE
	echo -e "python3 $VCFTOCSV SimplifyVCF -inVCF $NAME_REP_ANNOVAR/${method}/annotation_simplified_${name}.${method}.hg19_multianno.vcf -toType table -out  $NAME_REP_ANNOVAR/Table/annotation_simplified_${name}.${method}.csv" >> $ANNOTATION_FILE
	python3 $VCFTOCSV SimplifyVCF -inVCF $NAME_REP_ANNOVAR/${method}/annotation_simplified_${name}.${method}.hg19_multianno.vcf -toType table -out  $NAME_REP_ANNOVAR/Table/annotation_simplified_${name}.${method}.csv
	date >> $ANNOTATION_FILE
}

# Analyse du fichier d'annotation brute (application des filtres)
function ANNOTATION_ANALYSE (){
	name=$1
	method=$2
	date >> $ANNOTATION_FILE
	VERIFY_FILE $DIRECTORY/$name/$NAME_REP_ANNOVAR/Table/annotation_simplified_${name}.${method}.csv
	echo -e "*********\nANNOTATION_ANALYSE\n*********"
	echo -e "python3 $TRAIT_ANNOT -d $DIRECTORY/$name/$NAME_REP_ANNOVAR/Table/ -f annotation_simplified_${name}.${method}.csv  -o processed_annotation_simplified_${name}.${method}.csv  -m ${method}" >> $ANNOTATION_FILE
	python3 $TRAIT_ANNOT -d $DIRECTORY/$name/$NAME_REP_ANNOVAR/Table/ -f annotation_simplified_${name}.${method}.csv  -o processed_annotation_simplified_${name}.${method}.csv  -m ${method}
	date >> $ANNOTATION_FILE
}

# Launch ANNOVAR annotation
function ANNOVAR (){
	name=$1
	method=$2
	# Creation repertory
	mkdir $NAME_REP_ANNOVAR/${method}
	echo -e "*********\n${method} - ANNOVAR vcf input\n*********\n" >> $ANNOTATION_FILE

	# Annotation
	date >> $ANNOTATION_FILE
	echo -e "$ANNOVAR/table_annovar.pl variant/${method}/${name}.${method}.vcf $ANNOVAR_DB -buildver hg19 -out $NAME_REP_ANNOVAR/${method}/annotation_simplified_${name}.${method} -remove -protocol refGene,cytoBand,cosmic92,cosmic89,avsnp138,gnomad211_exome,clinvar_20200316,dbnsfp35a,IARC,icgc21 -operation gx,r,f,f,f,f,f,f,f,f -nastring . -thread 16 -polish -vcfinput -xref $ANNOVAR_DB/hg19_refGene.txt" >> $ANNOTATION_FILE
	$ANNOVAR/table_annovar.pl variant/${method}/${name}.${method}.vcf $ANNOVAR_DB -buildver hg19 -out $NAME_REP_ANNOVAR/${method}/annotation_simplified_${name}.${method} -remove -protocol refGene,cytoBand,cosmic92,cosmic89,avsnp138,gnomad211_exome,clinvar_20200316,dbnsfp35a,IARC,icgc21 -operation gx,r,f,f,f,f,f,f,f,f -nastring . -thread 16 -polish -vcfinput -xref $ANNOVAR_DB/hg19_refGene.txt
	date >> $ANNOTATION_FILE
	# VCT to CSV
	VCFTOTABLE $name $method
	# Filter
	ANNOTATION_ANALYSE $name $method

	VERIFY_FILE $DIRECTORY/$name/$NAME_REP_ANNOVAR/Table/filtered_processed_annotation_simplified_${name}.${method}.csv
}

# ***
# Main annnotation function
# ***
function LAUNCH_ANNOTATION () {
	name=$1
	# Creation du repertoire si ce n'est déja fait
	mkdir $QUALITY/$name
	# Ecriture terminal et fichier log principale
	echo -e "Lancement Annotation" >> $LOG
	PATIENT_ID_REMINDER $name
	# Fichier de sortie
	ANNOTATION_FILE=$DIRECTORY/$name/log_annotation_Routine.txt
	echo -e "************************************************\n" > $ANNOTATION_FILE
	echo -e "Lancement annotation ANNOVAR:" >> $ANNOTATION_FILE
	date >> $ANNOTATION_FILE
	# **************************
	# ANNOVAR
	# **************************

	NAME_REP_ANNOVAR=Annotation_ANNOVAR
	mkdir $NAME_REP_ANNOVAR
	# GATK
	# ****************
	ANNOVAR $name $METHODE1
	# Mutect
	# ****************
	ANNOVAR $name $METHODE2
	# Varscan
	# ****************
	ANNOVAR $name $METHODE3
	# Pindel
	# ****************
	# Only in CALR FLT3 NPM1 KMT2A
	ANNOVAR $name $METHODE4

	date >> $ANNOTATION_FILE
	# Dictionnary of information by biologist artefact
	echo -e "**********************************" >> $ANNOTATION_FILE
	echo -e "Fusion ANNOTATION" >> $ANNOTATION_FILE
	date >> $ANNOTATION_FILE
	dico_annotation_exists=$($EXIST_FILE ${DICT_ANNOTATION} )

	echo -e $dico_annotation_exists
	if [ $dico_annotation_exists == "No_file" ]
		then
		# Creation Dictionnary of annotation if don't exist
		echo -e "python3 $TRAIT_ANNOT -a $BASE_ARTEFACT -t ${BASE_TRANSCRIT} -c ${DICT_ANNOTATION}" >>$ANNOTATION_FILE
		python3 $TRAIT_ANNOT -a $BASE_ARTEFACT -t $BASE_TRANSCRIT -c $DICT_ANNOTATION
	fi
	# ***************
	# All fusion annotation
	# Annotation simple
	echo -e "python3 $TRAIT_ANNOT -d  $DIRECTORY/$name/$NAME_REP_ANNOVAR/Table/ -f filtered_processed_annotation_simplified_${name}.${METHODE1}.csv,filtered_processed_annotation_simplified_${name}.${METHODE2}.csv,filtered_processed_annotation_simplified_${name}.${METHODE3}.csv,filtered_processed_annotation_simplified_${name}.${METHODE4}.csv  -o Fusion_annotation_simplified_${name} -i ${DICT_ANNOTATION}  -r ${RUN_ID} -m All" >> $ANNOTATION_FILE
	python3 $TRAIT_ANNOT -d  $DIRECTORY/$name/$NAME_REP_ANNOVAR/Table/ -f filtered_processed_annotation_simplified_${name}.${METHODE1}.csv,Filter_Fichier_annotation_simplified_${name}.${METHODE2}.csv,Filter_Fichier_annotation_simplified_${name}.${METHODE3}.csv,Filter_Fichier_annotation_simplified_${name}.${METHODE4}.csv  -o Fusion_annotation_simplified_${name}  -i $DICT_ANNOTATION -r ${RUN_ID} -m All

	# ***************
	# Insert Dictionnary
	# Insertion des données en ne supprimant pas les variants avec une AF_raw > 0.01
	echo -e "Ecriture Dictionnaire :" >> $ANNOTATION_FILE
	# Verification si le dictionnaire existe ou non à ce chemin
	dico_exist=$($EXIST_FILE ${DICT})

	# Si le fichier n'existe pas ou il est vide : Creation du dictionnaire
	echo "**************************************"
	if [ $dico_exist == "No_file" ]
		then
		echo -e "python3 $DICTIONNARY -d $DIRECTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Fusion_annotation_simplified_${name}.csv -o ${DICT} -c True" >> $ANNOTATION_FILE
		python3 $DICTIONNARY -d $DIRECTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Fusion_annotation_simplified_${name}.csv -o $DICT -c True
	# S'il existe deja implementation du dictionnaire
	else
		echo -e "python3 $DICTIONNARY -d $DIRECTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Fusion_annotation_simplified_${name}.csv -o ${DICT}" >> $ANNOTATION_FILE
		python3 $DICTIONNARY -d $DIRECTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Fusion_annotation_simplified_${name}.csv -o $DICT
	fi
	# ***************
	# Creation du fichier d'annotation file
	echo -e "python3 $TRAIT_ANNOT -d  $DIRECTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Final_Fusion_annotation_simplified_${name}.csv  -o $QUALITY/$name/Annotation_patient_${name}.csv  -m Statistic -s $STAT_DICT -i $DICT_ANNOTATION" >> $ANNOTATION_FILE
	# Retourne un fichier d'annotation finale avec la colonne de Freq database
	python3 $TRAIT_ANNOT -d  $DIRECTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Final_Fusion_annotation_simplified_${name}.csv  -o $QUALITY/$name/Annotation_patient_${name}.csv  -m Statistic -s $STAT_DICT -i $DICT_ANNOTATION
	# Verification de l'écriture du fichier
	VERIFY_FILE $QUALITY/$name/Annotation_patient_${name}.csv
	date >> $ANNOTATION_FILE
}



# Choix des patients à analyser pour le cas SELECTION=UNIQUE
function ID_SELECTION (){
	listepatient=$(cut -d "," -f1 $TPL | grep -e "[-]")
	echo -e "Patient list:\n${listepatient}"
	echo -e "************"
	echo -e "Rentrez Identifiant patient ou la selection de patient (STOP pour arreter le programme)"
	echo -e "Input one or more (space-separated) patient ids"
	read patient
	echo -e "************"
}

# PATIENT_ANALYSIS user input menu
function PATIENT_ANALYSIS_MENU () {
	# getting name parameter
	name=$1
	# Going back to output directory
	cd $DIRECTORY
	# Menu
	case $ANALYSIS in
		Quality_bam )
			echo -e "***\nLaunching $ANALYSIS\n***   " >> $LOG
			QC $name
			;;
		Variant_calling )
			echo -e "***\nLaunching $ANALYSIS\n***   " >> $LOG
			VARIANT_CALLING $name
			;;

		Annotation )
			echo -e "***\nLaunching $ANALYSIS\n***   " >> $LOG
			LAUNCH_ANNOTATION $name
			;;

		Analyse_bam )
			echo -e "***\nLaunching $ANALYSIS\n***   " >> $LOG
			VARIANT_CALLING $name
			LAUNCH_ANNOTATION $name
			;;

		All )
			echo -e "***\nLaunching $ANALYSIS\n***   " >> $LOG
			QC $name
			VARIANT_CALLING $name
			LAUNCH_ANNOTATION $name
			;;
		* )
			echo "ERREUR de saisie lors du programme : $0 Arret du programme " >> $LOG
			exit 0
			;;
	esac
}

# Once the analysis mode is selected in PATIENT_ANALYSIS_MENU, $
# LAUNCH_PATIENT_ANALYSIS is used to select on which files
# the analysis will be performed
function LAUNCH_PATIENT_ANALYSIS () {

	# Demultiplexing only
	if [ "$ANALYSIS" = "Demultiplexing" ]; then
		echo "Launching bcl2fastq on $ANALYSIS" >> $LOG
		FASTQ_SETUP

		# Program exit
		date >> $LOG
		exit 0
	# Getting the last generated dictionnary before the run
# IDEA: partir d'une même pour les patients lancés en même temps
	elif [ "$ANALYSIS" = "Annotation" ] || [ "$ANALYSIS" = "All" ] ; then
		# Statistic one time in dictionnary
		echo -e "python3 $DICTIONNARY -o ${DICT} -s True -outstat ${STAT_DICT}" >> $LOG
		python3 $DICTIONNARY -o $DICT -s True -outstat $STAT_DICT
		VERIFY_FILE $STAT_DICT
	else
		echo "Aucune preanalyse a effectué pour : $ANALYSIS" >> $LOG
	fi
	# *****************************
	# Lancement du code
	# *****************************
	# Si l'utilisateur travaille sur tous les patients du RUN
	if [ "$SELECTION" = "all" ]
		then
		# Setting up a counter
		i=0
		for name in $(cut -d "," -f1 $TPL | grep -e "-")
		do
			echo -e "All :${name}" >> $LOG
			PATIENT_ANALYSIS_MENU $name
			echo -e "Patient ${i} terminé ${name} pour l'analyse ${ANALYSIS}"

			date >> $LOG
			i=`expr $i + 1`
		done
		echo "Vous avez analysé $i patient pour l'analyse ${ANALYSIS}:" >> $LOG
		echo -e "**********************************************************************\n" >> $LOG
	# *****************************
	# if the user works on one or a few patients in the set
	elif [ "$SELECTION" = "selection" ]; then
		# Select patient id
		ID_SELECTION
		while test  "$patient" != "STOP" ; do
			listepatient=$(echo $patient | tr "," "\n")

			for namer in $listepatient

			do
				name=$(cut -d "," -f1 $TPL | grep -e "^${namer}")
				echo "patient ${name}:"
				echo "patient ${name}:" >> $LOG

				PATIENT_ANALYSIS_MENU $name
				echo -e "************\n Fin de l'analyse pour le patient ${name}\n************"
				date >> $LOG
			done
			# Select patient ID again
			ID_SELECTION
		done
	# ************************
	# Variable selection error
	# ************************
	else
		echo "$0 FoxNGS ERROR: Input error ${SELECTION} not found" >> $LOG
		date >> $LOG
		exit 1
	fi
	echo -e "***********\nEnd of program\n***********"
}

# ******************************************************************************************
# ***************************************  MAIN  *******************************************
# ******************************************************************************************
# Lancement interface utilisateur
INTERFACE
# Launching pipeline
echo -e "**********************************************************" >> $LOG
echo -e "Lauching ${ANALYSIS}:" >> $LOG
date >> $LOG
PARAMETER_REMINDER
# Tool update
DATABASE
# Lauching analysis
LAUNCH_PATIENT_ANALYSIS
echo "Analysis - ${ANALYSIS} complete" >> $LOG
date >> $LOG
exit 0