 #!/bin/bash
# Auteur: Maxime Renard
# Bioinformatique
# Laboratoire d'hématologie
# Pipeline Recherche de variant sur le panel de 70 gènes pour le SMP
# Dernière Modification: 04/08/20
# But :
# Verification et notation de la prise en compte de NM 007313 à NM05157 de ABL1. 
# Une réféencxe plus connue dans la littérature
# suppresion variable site qui ne sont plus utilisé (exac, 1000G gardé en achivvage)
echo "#############################################################"
echo "#-----------------------------------------------------------#"
echo "#-                                                         -#"
echo "#-                 Pipeline bioinformatique SMPHD          -#"
echo "#-                   Données SURESELECQTX                  -#"
echo "#-                       Version 3.0                       -#"
echo "#-                                                         -#"
echo "#-----------------------------------------------------------#"
echo "#############################################################"

# ****************************************************************************************** 
# FUNCTION
# ******************************************************************************************
# ******************************** APPEL PROGRAMME *****************************************

# help
function HELP {
	echo -e "Pipeline bioinformatique SMPHD"
	echo -e "Lancement du pipeline HD"
	echo -e "./pipeline_SMPHD_v3.0.sh"
	echo -e "Un menu demandera les paramètres à rentrer"
	echo -e "Version 3.0"
	echo -e "Cette nouvelle version contient la création d'un dictionnaire"
	echo -e "La mise à jour de Cosmic 91 à la place de Cosmic 90 et de clinvar"
	echo -e " Si une modification de database ou d'appel de software a été modifié."
	echo -e "Il est nécessaire de créer un autre dictionnaire"
	echo -e "Question qualité:"
	echo -e "Dans le répertoire /home/t-chu-027/BioNGSTools/ version_Tools or Version_genome_REFERENCE"
	echo -e "Vérification FASTQC: Regarder qualité du Séquencage"
	echo -e "Vérification QUALITY coverage: Look a file html"
	exit 1
}

# STOP PROGRAM
function STOP_PROGRAM () {
	echo "Voulez vous continuer le programme ?"
	read reponse
	case $reponse in
	  [yYoO]*) echo "Suite du programme";;
	  [nN]*) echo "$0 arret du programe par l'utilisateur ;-)"
			 exit 0;;
	  *) echo "ERREUR de saisie"
		 exit 1;;
	esac
}

# Information sur l'utilisateur
function LOG_USER () {  # classic logger
	local prefix="[Pipeline $0 $(date +%Y/%m/%d\ %H:%M:%S)]: par $USER " 
	echo "${prefix} $@" >&2
	echo "${prefix} $@" >>$LOG
}

# Selection du génome de référence
function CHOIX_GENOME () {
	echo -e "Choix 1 : hg19 (GRCh37.p13)\n"
	echo -e "Choix 2 : hg38 (GRCh37.p12)\n"
	echo -e "Choix 3 : DEFAULT - hg19.fa (Version ancienne v1 pipeline)"
	echo -e "Rentrer votre choix:"
	echo "********"	
}
# CREATION DES FICHIERS DE SORTIE
function CREATION_RECHERCHE_FILE () {
	REPERTORY=${RUN}/${SORTIE}
	mkdir $REPERTORY 

	LOG=$REPERTORY/log_pipeline_SMPHD_v3.0.txt
	WORKFLOW=$REPERTORY/Workfow.xml
	CONFIGURATION=$REPERTORY/Configuration_file
	DESIGN=$REPERTORY/Experimental_Design
	# Localisation file
	# Exit report
	QUALITY=/media/t-chu-027/Elements/Result_NGS/$NAME_RUN
	COUV=/media/t-chu-027/Elements/Result_NGS/$NAME_RUN
	mkdir $QUALITY
	mkdir $COUV
	# Recherche automatique des fichiers brutes
	# Recuperation des data ngs brutes
	foldernamebcl=$(ls | grep -vE "^[a-z]|^[A-Z]")
	# Chemin vers le fichier brute de séquençage bcl
	BCL=$RUN/$foldernamebcl	
	# Récupération du fichier du RUN contenant les échantillons
	# Note : Ce fichier doit contenir le nom du run
	TPL=$BCL/SureSelect_QXT_$NAME_RUN.csv
}
# Selecion des fichiers du génome de référence choisi par l'utilisateur
function GENOME_REFERENCE () {

 case $REFERENCE in
  hg19) 
		echo "Vous travaillez avec le chromosome de référence GCF_000001405.25_GRCh37.p13 dernère modification le 19/04/17 " >> $LOG
		BWA_REF=/media/t-chu-027/DATAPART2/Database/Genome/hg37-p13/reference/GCF_000001405.25_GRCh37.p13_ref
		BWA_FASTA=/media/t-chu-027/DATAPART2/Database/Genome/hg38-p12/reference/GCF_000001405.25_GRCh37.p13_genomic.fna
		;;
  hg38 ) 
		echo -e "Vous travaillez avec le chromosome de référence GCF_000001405.38_GRCh38.p13 dernière modification le 25/06/2019 " >> $LOG
		BWA_REF=/media/t-chu-027/DATAPART2/Database/Genome/hg38-p12/reference/GCF_000001405.39_GRCh38.p13_ref
		BWA_FASTA=/media/t-chu-027/DATAPART2/Database/Genome/hg38-p12/reference/GCF_000001405.39_GRCh38.p13_genomic.fna  	
		;;
  DEFAULT | *) 
		echo -e "Vous travaillez avec le chromosome de référence hg19 datant de 2013. " >> $LOG
		BWA_REF=/media/t-chu-027/DATAPART2/Database/Genome/hg19/hg19bwaidx
		BWA_FASTA=/media/t-chu-027/DATAPART2/Database/Genome/hg19/hg19.fa   
		;;
	esac
	echo -e "*********************************"
	echo -e "Note : Genome de référence "
	echo -e "Pour en savoir plus sur la préparation d'un génome regarder le script /home/t-chu-027/Bureau/Recherche/Pipeline/Preparation_genome/Preparation_genome.sh"
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
		[yYoO]*) echo -e "Les fichiers de qualité sont correcte " >> $LOG;;
		[nN]*) echo -e "Les fichiers de qualité ne sont pas correctes : Arret programme. " >> $LOG
				echo "$0 arret du programe - Fichier d'analyse non valide ;-)"
			 exit 0;;
	  DEFAULT | *) echo "ERREUR de saisie Par DEFAULT OK" >> $LOG
	esac
}

# Rappel du patient et entrée dans le répertoire de sortie du patient
function RAPPEL_PATIENT () {
	name=$1
	echo -e "***********************************************"
	echo "Nom du patient analysé ${name}:" 
	echo -e "***************************************\n" >> $LOG
	echo "Nom du patient analysé ${name}:" >> $LOG
	cd $REPERTORY/$name
}
# *********************************************
#  Logiciels et Outils utilisés
# *********************************************
function DATABASE () {
	# See BioNGSTools/version_Tools pour plus d'information sur la version des outils
	echo -e "Fichier qualité des outils dans /home/t-chu-027/Bureau/Recherche/Diagnostique_routine/Qualite_kalilab" >> $LOG
	# Tools
	BEDTOOLS=~/BioNGSTools/bedtools-2.27.0/bedtools2/bin/bedtools
	PICARD=~/BioNGSTools/picard/picard.jar
	
	# Analyse qualité - Script R
	RSCRIPT=~/Bureau/Recherche/Pipeline/SMPHD_v3.0/R_Quality_SMPHD_v3.0.Rmd
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
	TRAIT_ANNOT=~/Bureau/Recherche/Pipeline/SMPHD_v3.0/Treatment_of_Annotation_v3.0.py
	EXIST_FILE=~/Bureau/Recherche/Script/exist_file.sh
	DICTIONNARY=~/Bureau/Recherche/Pipeline/SMPHD_v3.0/database_dictionnary_v3.0.py
	# ********************************************
	# Database 
	# Fichier
	BED=/media/t-chu-027/DATAPART2/Database/Fichier_intersection_bed/Sure_Select_design/SureSelect-HEMATO-v5.sorted.bed
	# Bed pour couverture et l'analyse qualité R
	BEDEXON=/media/t-chu-027/DATAPART2/Database/Fichier_intersection_bed/Analyse_coverage/DESIGN-FH-EXONS-gene_panel.bed
	# Variant
	# Fichier Bed Pindel
	BED_PINDEL=/media/t-chu-027/DATAPART2/Database/Variant/Pindel_search_CALR-9_FLT3-14-15.bed
	DBSNP=/media/t-chu-027/DATAPART2/Database/DB/dbsnp_138.hg19.vcf

	# Annotation
	ANNOVAR=~/BioNGSTools/annovar/
	ANNOVAR_DB=/media/t-chu-027/DATAPART2/Database/humandb_annovar

	ANNOTATION_REP=/media/t-chu-027/DATAPART2/Database/Annotation
	BASE_TRANSCRIT=$ANNOTATION_REP/Transcript_reference/Liste_genes_transcript_27-07-20.csv
	BASE_ARTEFACT=$ANNOTATION_REP/Artefact/Base_artefact_120220.csv
	# Dictionnaire annotation 
	# Changement de nom 
	DICT_ANNOTATION=$ANNOTATION_REP/Database_annotation_27_07_20-v3.0.json
	# Exit Database Result
	STAT_DICT=/media/t-chu-027/Elements/Result_NGS/Dictionnary_database/Project_FIM_statistic_dictionnary.csv
	DICT=/media/t-chu-027/Elements/Result_NGS/Dictionnary_database/Project_FIM_Database_variant.json
	# Methode d'appel de variant
	METHODE1="GATK"
	METHODE2="Mutect2"
	METHODE3="Varscan"
	METHODE4="Pindel"
}

# Menu des analyses des fichiers de NGS
function CHOIX_ANALYSE () {
	echo "***********************************"
	echo -e "Choix 1 : Demultiplexage (1ère étape : Preparation des fichiers FASTQ par patient)\n"
	echo -e "Choix 2 : Quality_bam (fastq.gz to bam) \n"
	echo -e "Choix 3 : Variant_calling (4 variant calling) \n" 
	echo -e "Choix 4 : Annotation (Annotation des variants, Filtres et Ecriture du Fichier d'annotation)\n" 
	echo -e "Choix 5 : Analyse_bam (Lancement des étapes de variant calling et d'annotation (Pour les Tests))\n"
	echo -e "Choix 6 : All (2ème étape : Lancement de l'analyse par patient \n du fichier brute à la création du rapport (Combinaison du choix 2 à 4))\n"
	echo "***********************************"
	echo "Entrez le nom du choix. (Exemple All)  "
	echo "Quel analyse souhaitez vous effectuer ? "
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
	echo "Pipeline d'analyse des Données SURESELECT SMPHD version 3.0 "
	echo "****************************************************"
	echo "Saisie des Données de lancement du RUN :"
	echo "********"
	echo -e "Rentrez votre nom d'utilisateur :" 
	echo "********" 
	read USER
	cd /media/t-chu-027/DATAPART1/RUN_NGS/SURESELECTQXT/
	exp=$(ls)
	echo -e "Nom des RUN Sureselect :\n${exp}"
	echo "***********************************"
	echo "Saisir le nom RUN :"
	read NAME_RUN
	# Choix du RUN à analyser
	echo "********"
	RUN=/media/t-chu-027/DATAPART1/RUN_NGS/SURESELECTQXT/$NAME_RUN
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
	read SORTIE
	# Recuperation de la localisation des fichiers d'interet
	CREATION_RECHERCHE_FILE
	echo "********"
	echo "***********************************"
	# Choix du génome de référence
	# Par default choix DEFAULT pour accelerer la procedure de lancement
	# Amélioration 4.0 Choix des genomes de référence (hg19,hg37 et hg38)
	REFERENCE="DEFAULT"
	# Information sur la localisation des génomes de référence
	GENOME_REFERENCE
	# Verification des étapes intermediaires
	# Par default il n'y en a pas 
	qualite="NO"
	# Choix de l'Analyse du Pipeline
	CHOIX_ANALYSE
	# Choix de la portée de l'analyse
	if [ "$ANALYSE" != "Demultiplexage" ]; then
		echo "***********************************"
		echo "Voulez vous lancer l'analyse sur tous les patients ou sur un patient?"
		echo -e "Choix 1 - all (Analyse demandé Sur tous les patients du RUN sélectioné\n" 
		echo -e "Choix 2 - unique (Une analyse sur 1 patient \nou une selection de plusieurs patients délimité par des,)\n Choix des patients dessous\n"
		echo -e "Rentrez le nom du choix (Exemple : all) "
		echo "********"
		read SELECTION

	fi
}

# ***********
# Rappel des paramètres
function RAPPEL_PARAMETRE () {
	echo -e "**********************************************************" >> $LOG 
	LOG_USER
	echo -e "Récapitulatif des arguments rentrés en paramètre" >> $LOG
	echo -e "Dossier sequence brute: ${BCL}" >> $LOG
	echo -e "Fiche NGS patient: ${TPL}" >> $LOG
	echo -e "RUN: ${NAME_RUN}" >> $LOG
	echo -e "Génome de référence: ${REFERENCE}" >> $LOG
	echo -e "***********************************">> $LOG
	echo -e "Répertoire de sortie des résultats:">> $LOG
	echo -e $REPERTORY >> $LOG
	echo -e "Analyse Effectué $ANALYSE " >> $LOG
	if [ "$ANALYSE" != "Demultiplexage" ]; then
		echo -e "Lancement sur patient ${SELECTION}" 		
	fi
}

# Fonction Gestion des erreurs
# Verification de l'existence des fichiers
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

# ******************************* ANALYSE PIPELINE *****************************************
# *******************************************************************
# Demultiplexage des données pour chaque patient
# Entrée Fichier bcl2fastq 	SORTIE Fichier fastq
# *******************************************************************
function PREPARATION_FASTQ () {
	cd $REPERTORY
	echo -e "**********************************************************************\n" >> $LOG
	echo "Génération des fichiers FASTQ à partir des BCL" >> $LOG
	date >> $LOG                                            
	echo "start de bcl2fastq" >> $LOG
	echo -e "bcl2fastq --runfolder $BCL --output-dir $REPERTORY/ --sample-sheet $TPL --use-bases-mask Y150,I8,I8,Y150 --no-lane-splitting -r 16 -p 16 -w 16" >> $LOG
	bcl2fastq --runfolder $BCL --output-dir $REPERTORY/ --sample-sheet $TPL --use-bases-mask Y150,I8,I8,Y150 --no-lane-splitting -r 16 -p 16 -w 16 
	echo "fin de bcl2fastq" >> $LOG
	date >> $LOG
	echo -e "**********************************************************************\n" >> $LOG
}

# **********************************************************
# Préparation fichier BAM
function LANCEMENT_QUALITY_BAM () {
	
	# Récupération paramètre 
	name=$1
	
	# Fichier log de sortie
	PREPARATION_BAM_FILE=$REPERTORY/$name/log_bam.txt
	
	# Création d'un répertoire pour chaque patient
	mkdir $REPERTORY/$name
	# Création d'un répertoire de sortie de résultat dans le bureau 
	mkdir $QUALITY/$name

	echo -e "**********************************************************************\n" > $PREPARATION_BAM_FILE
	date > $PREPARATION_BAM_FILE
	echo -e "Génération des Fichiers d'alignement BAM pour ${name}\n" >> $PREPARATION_BAM_FILE
	mv $name/*.fastq.gz .

	# Récupération des fichiers fastq sens R1 et R2 correspondant à un identifiant
	R1=$(ls | grep $name | grep _R1_)
	R2=$(ls | grep $name | grep _R2_)

	echo "R1" >> $PREPARATION_BAM_FILE
	echo $R1 >> $PREPARATION_BAM_FILE
	echo "R2" >> $PREPARATION_BAM_FILE
	echo $R2 >> $PREPARATION_BAM_FILE

	# Déplacement des FASTQ dans le fichier du patient
	mv $R1 $REPERTORY/$name
	mv $R2 $REPERTORY/$name

	#On rentre dans le fichier du  patient
	echo -e "Lancement Quality bam" >> $LOG
	RAPPEL_PATIENT $name
	echo "Nom du patient analysé ${name}:" >> $PREPARATION_BAM_FILE
	echo $name >> $PREPARATION_BAM_FILE
	VERIFY_FILE $REPERTORY/$name/$R1
	VERIFY_FILE $REPERTORY/$name/$R2
	# *************************************************
	# Elaboration du FASTQC pour la patient :
	# *************************************************
	# Génération du fichier qualité
	echo -e "fastqc -o . $R1 $R2 -t 16" >> $PREPARATION_BAM_FILE
	fastqc -o . $R1 $R2 -t 16
	R1name=$(echo $R1 |cut -f1 -d.)
	R2name=$(echo $R2 |cut -f1 -d.)
	echo $R1name
	echo $R2name
	# Extension
	html="_fastqc.html"
	# copie des fichiers d'analyse fastqc vers le repertoire qualite
	cp  $REPERTORY/$name/$R1name$html $REPERTORY/$name/$R2name$html $QUALITY/$name/
	# Verification fichier
	if [ $qualite = "OK" ]
		then
			# Commande pour validation
			echo -e "Attente de la validation du fichier générer par fastqc (Quality Control):"
			echo -e "Afficher OK:  si le fichier de qualité est correcte\nNO: si incorrecte\n ********\n"
			firefox $REPERTORY/$name/$R1name$html
			firefox $REPERTORY/$name/$R2name$html
			VALIDATION
			echo -e "*******"
	fi

	# Suppression des fichiers brutes en attente fichier intermediaire
	rm -dr $REPERTORY/$name/*fastqc.zip 
	#Creation d'un fichier temporaire: Stockage Fichier SAM à BAM de préparation
	mkdir $REPERTORY/$name/tmp
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
	bwa mem -t 16 -R '@RG\tID:C5-${name}\tPL:illumina\tPU:HXXX\tLB:Solexa\tSM:C5-${name}' $BWA_REF $R1 $R2 -o tmp/${name}.sam 
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
	total=$(samtools view -h -c tmp/${name}.sort.bam )
	echo -e "*******************************" >> $PREPARATION_BAM_FILE
	echo -e "Total mapped: $total" >> $PREPARATION_BAM_FILE
	echo -e "*******************************" >> $PREPARATION_BAM_FILE
	echo -e " Keep only mapped\n" >> $PREPARATION_BAM_FILE
	echo -e "samtools view -F 0x4 -@ 16 -h -b tmp/${name}.sort.bam >tmp/${name}.sort_mapped.bam" >> $PREPARATION_BAM_FILE
	samtools view -F 0x4 -h -@ 16 -b tmp/${name}.sort.bam >tmp/${name}.sort_mapped.bam
	echo -e "samtools view -f 0x800 -@ 16 -h -b tmp/${name}.sort.bam >tmp/${name}.sort_2048.bam" >> $PREPARATION_BAM_FILE
	samtools view -f 0x800 -@ 16 -h -b tmp/${name}.sort.bam >tmp/${name}.sort_2048.bam
	echo -e "*******************************" >> $PREPARATION_BAM_FILE
	mapped=$(samtools view -h -@ 14 -c tmp/${name}.sort_mapped.bam)
	echo -e "Only mapped $mapped" >> $PREPARATION_BAM_FILE
	echo -e "*******************************" >> $PREPARATION_BAM_FILE
	unmapped=$(samtools view -f 0x4 -h -@ 10 -c -b tmp/${name}.sort.bam)
	chimeric=$(samtools view  -f 0x800 -h -@ 10 -c -b tmp/${name}.sort.bam)
	echo -e "Only unmapped: $unmapped" >> $PREPARATION_BAM_FILE
	echo -e "Only chimeric: $chimeric" >> $PREPARATION_BAM_FILE
	echo -e "*******************************" >> $PREPARATION_BAM_FILE
	date >> $PREPARATION_BAM_FILE
	# Intersection avec le fichier du design 
	# avec les région intronique aussi  pour déceler les mutations de type épissage en plus
	echo -e "$BEDTOOLS intersect -a tmp/${name}.sort_mapped.bam -b $BED  > tmp/${name}-on-target.bam" >> $PREPARATION_BAM_FILE
	$BEDTOOLS intersect -a tmp/${name}.sort_mapped.bam -b $BED  > tmp/${name}-on-target.bam
	VERIFY_FILE	tmp/${name}-on-target.bam
	# Génération des statistiques: compte rendu qualité
	echo -e "Compte rendu qualité"
	echo -e "samtools flagstat -@ 8 tmp/${name}-on-target.bam > tmp/${name}.bam.sort.stat" >> $PREPARATION_BAM_FILE
	samtools flagstat -@ 10 tmp/${name}-on-target.bam > tmp/${name}.bam.sort.stat 		
	echo -e "Reindexation" >> $PREPARATION_BAM_FILE
	echo -e "samtools index  -@ 8 -b tmp/${name}-on-target.bam" >> $PREPARATION_BAM_FILE
	samtools index -@ 14 -b tmp/${name}-on-target.bam

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
	echo -e "$BEDTOOLS coverage -sorted -a $BEDEXON -b $REPERTORY/$name/${name}.sort.dupmark.bam -d >$REPERTORY/$name/${name}_couverture_analyse.bed" >> $PREPARATION_BAM_FILE
	$BEDTOOLS coverage -sorted -a $BEDEXON -b $REPERTORY/$name/${name}.sort.dupmark.bam -d >$REPERTORY/$name/${name}_couverture_analyse.bed	
			
	# ****************************************************
	# Analyse de la qualité de coverage Rmarkdown 
	echo -e "R -e rmarkdown::render('${RSCRIPT}', 
		params = list( directory='$(pwd)',file='${name}_couverture_analyse.bed',user='${USER}',pipeline='$0',output='${COUV}/Statistic_couverture.csv',output_gene='/media/t-chu-027/Elements/Result_NGS/Stat_gene/Statistic_couverture_gene.csv'),
		output_file='$(pwd)/${name}_couverture_analyse.bed.html')" >> $PREPARATION_BAM_FILE
	R -e "rmarkdown::render('${RSCRIPT}', params = list(directory='$(pwd)',file='${name}_couverture_analyse.bed',user='${USER}',pipeline='${0}',output='${COUV}/Statistic_couverture.csv',output_gene='/media/t-chu-027/Elements/Result_NGS/Stat_gene/Statistic_couverture_gene.csv'),output_file='$(pwd)/${name}_couverture_analyse.bed.html')"

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
	rm -dr tmp/*sam tmp/${name}.sort_mapped.bam tmp/${name}.bam tmp/*on-target.bam* tmp/*2048* 
	
	# Copie des fichiers BAM
	cp ${name}.sort.dupmark.bam ${name}.sort.dupmark.bam.bai $QUALITY/$name/
	
}

# *****************************************
# Partie Variant calling
#

# Main function for variant calling
function LANCEMENT_VARIANT_CALLING () {

	# Récupération paramètre 
	name=$1
	mkdir $QUALITY/$name
	# Fichier de sortie 
	VARIANT_FILE=$REPERTORY/$name/logvariantcalling.txt
	# function RAPPEL patient
	echo -e "Lancement Variant_calling" >> $LOG
	RAPPEL_PATIENT $name
	
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
	VERIFY_FILE variant/${METHODE3}/${name}.${METHODE3}.vcf
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
	echo -e "$PINDEL -f $BWA_FASTA -i $REPERTORY/$name/variant/${METHODE4}/config_file_pindel.txt -j $BED_PINDEL -T 14 -o variant/${METHODE4}/ITD_${name}" >> $VARIANT_FILE
	$PINDEL -f $BWA_FASTA -i $REPERTORY/$name/variant/${METHODE4}/config_file_pindel.txt -j $BED_PINDEL -T 14 -o variant/${METHODE4}/ITD_${name}

	#D: Deletion + TD:Tandem Duplication + INV: Inversion + INS short and long insertion
	echo -e "$PINDEL_VCF -p variant/${METHODE4}/ITD_${name}_D -r $BWA_FASTA  -R x -d 00000000 -G  -v variant/${METHODE4}/${name}.${METHODE4}.vcf" >> $VARIANT_FILE
	$PINDEL_VCF -P variant/${METHODE4}/ITD_${name} -r $BWA_FASTA  -R x -d 00000000 -G -v variant/${METHODE4}/${name}.${METHODE4}.vcf
	VERIFY_FILE variant/${METHODE4}/${name}.${METHODE4}.vcf
	date >> $VARIANT_FILE
	echo -e "**********************************************************************"  >> $VARIANT_FILE
}

# *****************************************
# Partie annotation 

# Conversion des vcf en table csv pour l'analyse 
function VCFToTable (){
	name=$1
	method=$2
	mkdir $NAME_REP_ANNOVAR/Table/

	echo -e "*********\nVCF to Table\n*********"
	date >> $ANNOTATION_FILE
	echo -e "python3 $VCFTOCSV SimplifyVCF -inVCF $NAME_REP_ANNOVAR/${method}/annotation_simple_${name}.${method}.hg19_multianno.vcf -toType table -out  $NAME_REP_ANNOVAR/Table/annotation_simple_${name}.${method}.csv" >> $ANNOTATION_FILE
	python3 $VCFTOCSV SimplifyVCF -inVCF $NAME_REP_ANNOVAR/${method}/annotation_simple_${name}.${method}.hg19_multianno.vcf -toType table -out  $NAME_REP_ANNOVAR/Table/annotation_simple_${name}.${method}.csv
	date >> $ANNOTATION_FILE
}

# Analyse du fichier d'annotation brute (application des filtres) 
function Annotation_analyse (){
	name=$1
	method=$2
	date >> $ANNOTATION_FILE
	VERIFY_FILE $REPERTORY/$name/$NAME_REP_ANNOVAR/Table/annotation_simple_${name}.${method}.csv
	echo -e "*********\nAnnotation_analyse\n*********"
	echo -e "python3 $TRAIT_ANNOT -d $REPERTORY/$name/$NAME_REP_ANNOVAR/Table/ -f annotation_simple_${name}.${method}.csv  -o Fichier_annotation_simple_${name}.${method}.csv  -m ${method}" >> $ANNOTATION_FILE
	python3 $TRAIT_ANNOT -d $REPERTORY/$name/$NAME_REP_ANNOVAR/Table/ -f annotation_simple_${name}.${method}.csv  -o Fichier_annotation_simple_${name}.${method}.csv  -m ${method}
	date >> $ANNOTATION_FILE
}
# Launch Annovar annotation
function Annovar (){
	name=$1
	method=$2
	# Creation repertory
	mkdir $NAME_REP_ANNOVAR/${method}
	echo -e "*********\n${method} - Annovar vcf input\n*********\n" >> $ANNOTATION_FILE 
	
	# Annotation 
	date >> $ANNOTATION_FILE
	echo -e "$ANNOVAR/table_annovar.pl variant/${method}/${name}.${method}.vcf $ANNOVAR_DB -buildver hg19 -out $NAME_REP_ANNOVAR/${method}/annotation_simple_${name}.${method} -remove -protocol refGene,cytoBand,cosmic91,cosmic89,gnomad211_exome,clinvar_20200316,dbnsfp35a,IARC,icgc21 -operation gx,r,f,f,f,f,f,f,f -nastring . -thread 16 -polish -vcfinput -xref $ANNOVAR_DB/hg19_refGene.txt" >> $ANNOTATION_FILE
	$ANNOVAR/table_annovar.pl variant/${method}/${name}.${method}.vcf $ANNOVAR_DB -buildver hg19 -out $NAME_REP_ANNOVAR/${method}/annotation_simple_${name}.${method} -remove -protocol refGene,cytoBand,cosmic91,cosmic89,gnomad211_exome,clinvar_20200316,dbnsfp35a,IARC,icgc21 -operation gx,r,f,f,f,f,f,f,f -nastring . -thread 16 -polish -vcfinput -xref $ANNOVAR_DB/hg19_refGene.txt 
	date >> $ANNOTATION_FILE
	# VCT to CSV
	VCFToTable $name $method
	# Filter
	Annotation_analyse $name $method

	VERIFY_FILE $REPERTORY/$name/$NAME_REP_ANNOVAR/Table/Filter_Fichier_annotation_simple_${name}.${method}.csv
}


# Function principal
function LANCEMENT_ANNOTATION () {
	name=$1
	# Creation du repertoire si ce n'est déja fait
	mkdir $QUALITY/$name
	# Ecriture terminal et fichier log principale
	echo -e "Lancement Annotation" >> $LOG
	RAPPEL_PATIENT $name
	# Fichier de sortie
	ANNOTATION_FILE=$REPERTORY/$name/log_annotation.txt
	echo -e "************************************************\n" > $ANNOTATION_FILE
	echo -e "Lancement annotation ANNOVAR:" >> $ANNOTATION_FILE
	date >> $ANNOTATION_FILE
	# **************************
	# ANNOVAR
	# **************************	

	NAME_REP_ANNOVAR=Annotation_Annovar
	echo $NAME_REP_ANNOVAR
	mkdir $NAME_REP_ANNOVAR
	# GATK
	# ****************
	Annovar $name $METHODE1
	# Mutect
	# ****************
	Annovar $name $METHODE2
	# Varscan
	# ****************
	Annovar $name $METHODE3
	# Pindel 
	# ****************
	# Only in CALR  
	Annovar $name $METHODE4

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
	echo -e "python3 $TRAIT_ANNOT -d  $REPERTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Filter_Fichier_annotation_simple_${name}.${METHODE1}.csv,Filter_Fichier_annotation_simple_${name}.${METHODE2}.csv,Filter_Fichier_annotation_simple_${name}.${METHODE3}.csv,Filter_Fichier_annotation_simple_${name}.${METHODE4}.csv  -o Fusion_annotation_simple_${name} -i ${DICT_ANNOTATION} -m All" >> $ANNOTATION_FILE
	python3 $TRAIT_ANNOT -d  $REPERTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Filter_Fichier_annotation_simple_${name}.${METHODE1}.csv,Filter_Fichier_annotation_simple_${name}.${METHODE2}.csv,Filter_Fichier_annotation_simple_${name}.${METHODE3}.csv,Filter_Fichier_annotation_simple_${name}.${METHODE4}.csv  -o Fusion_annotation_simple_${name}  -i $DICT_ANNOTATION -m All

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
		echo -e "python3 $DICTIONNARY -d $REPERTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Fusion_annotation_simple_${name}.csv -o ${DICT} -c True" >> $ANNOTATION_FILE
		python3 $DICTIONNARY -d $REPERTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Fusion_annotation_simple_${name}.csv -o $DICT -c True
	# S'il existe deja implementation du dictionnaire
	else
		echo -e "python3 $DICTIONNARY -d $REPERTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Fusion_annotation_simple_${name}.csv -o ${DICT}" >> $ANNOTATION_FILE
		python3 $DICTIONNARY -d $REPERTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Fusion_annotation_simple_${name}.csv -o $DICT
	fi
	# ***************
	# Creation du fichier d'annotation file
	echo -e "python3 $TRAIT_ANNOT -d  $REPERTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Final_Fusion_annotation_simple_${name}.csv  -o $QUALITY/$name/Annotation_patient_${name}.csv  -m Statistic -s $STAT_DICT -i $DICT_ANNOTATION" >> $ANNOTATION_FILE
	# Retourne un fichier d'annotation finale avec la colonne de Freq database
	python3 $TRAIT_ANNOT -d  $REPERTORY/$name/$NAME_REP_ANNOVAR/Table/ -f Final_Fusion_annotation_simple_${name}.csv  -o $QUALITY/$name/Annotation_patient_${name}.csv  -m Statistic -s $STAT_DICT -i $DICT_ANNOTATION
	# Verification de l'écriture du fichier
	VERIFY_FILE $QUALITY/$name/Annotation_patient_${name}.csv
	date >> $ANNOTATION_FILE
}



# Choix des patients à analyser pour le cas SELECTION=UNIQUE
function CHOIX_IDENTIFIANT (){
	listepatient=$(cut -d "," -f1 $TPL | grep -e "[-]") 
	echo -e "Liste patient:\n${listepatient}" 
	echo -e "************"
	echo -e "Rentrez Identifiant patient ou la selection de patient (STOP pour arreter le programme)"
	read patient
	echo -e "************"
}

# MENU de lancement analyse_patient
function MENU_ANALYSE_PATIENT () {
	# Récupération paramètre 
	name=$1
	# Retour au répertoire de sortie
	cd $REPERTORY
	# Menu 
	case $ANALYSE in 
		Quality_bam ) 
			echo -e "***\nLancement Analyse $ANALYSE\n***   " >> $LOG
			LANCEMENT_QUALITY_BAM $name
			;;
		Variant_calling )
			echo -e "***\nLancement Analyse $ANALYSE\n***   " >> $LOG
			LANCEMENT_VARIANT_CALLING $name
			;;
		
		Annotation ) 
			echo -e "***\nLancement Analyse $ANALYSE\n***   " >> $LOG
			LANCEMENT_ANNOTATION $name
			;;

		Analyse_bam )
			echo -e "***\nLancement Analyse $ANALYSE\n***   " >> $LOG
			LANCEMENT_VARIANT_CALLING $name
			LANCEMENT_ANNOTATION $name
			;;

		All )
			echo -e "***\nLancement Analyse $ANALYSE\n***   " >> $LOG
			LANCEMENT_QUALITY_BAM $name
			LANCEMENT_VARIANT_CALLING $name 
			LANCEMENT_ANNOTATION $name
			;;
		* ) 
			echo "ERREUR de saisie lors du programme : $0 Arret du programme " >> $LOG
			exit 0
			;;
	esac
}

# Lancement de l'analyse du patient en prenant en compte le mode de séléction
function LANCEMENT_ANALYSE_PATIENT () {

	# Fonction Demultiplexage only
	if [ "$ANALYSE" = "Demultiplexage" ]; then
		echo "Lancement bcl2fastq pour $ANALYSE" >> $LOG
		PREPARATION_FASTQ
		# Sortie du programme
		date >> $LOG
		exit 0
	# Recuperation du dernier dictionnaire avant le lancement du RUN
	# Idée partir d'une même pour les patients lancés en même temps
	elif [ "$ANALYSE" = "Annotation" ] || [ "$ANALYSE" = "All" ] ; then
		# Statistic one time in dictionnary
		echo -e "python3 $DICTIONNARY -o ${DICT} -s True -outstat ${STAT_DICT}" 
		python3 $DICTIONNARY -o $DICT -s True -outstat $STAT_DICT
		VERIFY_FILE $STAT_DICT
	else
		echo "Aucune preanalyse a effectué pour : $ANALYSE" >> $LOG
	fi
	# *****************************
	# Lancement du code 
	# *****************************
	# Si l'utilisateur travaille sur tous les patients du RUN
	if [ "$SELECTION" = "all" ]
		then
		# Mise en place compteur
		i=0
		for name in $(cut -d "," -f1 $TPL | grep -e "-") 
		do
			echo -e "All :${name}" >> $LOG
			MENU_ANALYSE_PATIENT $name
			echo -e "Analyse du patient ${i} terminé ${name} pour l'analyse ${ANALYSE}"

			date >> $LOG
			i=`expr $i + 1`
		done
		echo "Vous avez analysé $i patient pour l'analyse ${ANALYSE}:" >> $LOG	
		echo -e "**********************************************************************\n" >> $LOG
	# *****************************
	# Si l'utilisateur travaille sur un patient ou une liste de patients en particulier
	elif [ "$SELECTION" = "unique" ]; then
		# Choix des identifiants du patients
		CHOIX_IDENTIFIANT
		while test  "$patient" != "STOP" ; do
			listepatient=$(echo $patient | tr "," "\n")

			for namer in $listepatient
			
			do
				name=$(cut -d "," -f1 $TPL | grep -e "^${namer}")
				echo "patient ${name}:"
				echo "patient ${name}:" >> $LOG
			
				MENU_ANALYSE_PATIENT $name
				echo -e "************\n Fin de l'analyse pour le patient ${name}\n************"
				date >> $LOG
			done
			# Relancement du choix des identifiants
			CHOIX_IDENTIFIANT
		done
	# *****************************
	# Error of variable selection		
	else
		echo "ERREUR de saisie lors du programme pour la variable ${SELECTION}: $0 Arret du programme" >> $LOG
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
# Lancement du pipeline
echo -e "**********************************************************" >> $LOG
echo -e "Lancement de l'analyse ${ANALYSE}:" >> $LOG
date >> $LOG
RAPPEL_PARAMETRE
# Mise à jour des outils
DATABASE
# Lancement des analyses
LANCEMENT_ANALYSE_PATIENT
echo "Analyse - ${ANALYSE} Réalisé avec succès:" >> $LOG
date >> $LOG
exit 0
