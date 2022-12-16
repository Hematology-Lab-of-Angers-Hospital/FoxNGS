function LOG_AGENT () {  # classic logger
	local prefix="[Pipeline $0 $(date +%Y/%m/%d\ %H:%M:%S)]: par $agent " 
	echo "${prefix} $@" >&2
	echo "${prefix} $@" >>$LOG
}

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

# Choix des patients à analyser pour le cas SELECTION=UNIQUE
function CHOIX_IDENTIFIANT (){
	listepatient=$(cut -d "," -f1 $TPL | grep -e "[-]") 
	echo -e "Liste patient:\n${listepatient}" 
	echo -e "************"
	echo -e "Rentrez Identifiant patient ou la selection de patient (STOP pour arreter le programme)"
	read patient
	echo -e "************"
}


function CHOIX_ANALYSE () {
	echo "***********************************"
	echo -e "Choix 1 : Demultiplexage (1ère étape : Preparation des fichiers FASTQ par patient)\n"
	echo -e "Choix 2 : Analyse (2ème étape : Lancement de l'analyse par patient \n du fichier brute à la création du rapport (Combinaison du choix 2 à 4))\n"
	echo "***********************************"
	echo "Entrez le nom du choix. (Exemple All)  "
	echo "Quel analyse souhaitez vous effectuer ? "
	echo "********"
	read ANALYSE
	echo "********"
}

function CREATION_RECHERCHE_FILE () {
	ANALYSIS_FOLDER=${RUN}/${SORTIE}
	mkdir $ANALYSIS_FOLDER
    mkdir $ANALYSIS_FOLDER/fastq

	LOG=$ANALYSIS_FOLDER/log_lancement_FoxNGS_Routine_v1_2_b.txt
	# Recherche automatique des fichiers brutes
	# Recuperation des data ngs brutes
	foldernamebcl=$(ls | grep -vE "^[a-z]|^[A-Z]")
	# Chemin vers le fichier brute de séquençage bcl
	BCL=$RUN/$foldernamebcl
	# Récupération du fichier du RUN contenant les échantillons
	# Note : Ce fichier doit contenir le nom du run
	TPL=$BCL/SureSelect_QXT_$NAME_RUN.csv
	CAT $TPL
}

function PREPARATION_FASTQ () {
	cd $ANALYSIS_FOLDER
	echo -e "**********************************************************************\n" >> $LOG
	date >> $LOG                                            
	echo "Lancement bcl2fastq" >> $LOG
    echo -e "singularity exec -B /media/tfli-0070/DATASAVE1 /home/tfli-0070/BioNGSTools/Singularity_database/bcl2fastq_docker.sif bcl2fastq --runfolder $BCL --output-dir $ANALYSIS_FOLDER/fastq --sample-sheet $TPL --use-bases-mask Y150,I8,I8,Y150 --no-lane-splitting -r 32 -p 32 -w 32" >> $LOG
    singularity exec -B /media/tfli-0070/DATASAVE1 /home/tfli-0070/BioNGSTools/Singularity_database/bcl2fastq_docker.sif bcl2fastq --runfolder $BCL --output-dir $ANALYSIS_FOLDER/fastq --sample-sheet $TPL --use-bases-mask Y150,I8,I8,Y150 --no-lane-splitting -r 32 -p 32 -w 32 
	echo "FIN bcl2fastq" >> $LOG
	date >> $LOG
	echo -e "**********************************************************************\n" >> $LOG
}

function INTERFACE () {
	echo "***************************************************************"
	echo "Pipeline d'analyse des Données SURESELECT FoxNGS version 1.2.b "
	echo "***************************************************************"
	echo "Saisie des Données de lancement du RUN :"
	echo "********"
	echo -e "Rentrez votre nom d'utilisateur :" 
	echo "********" 
	read USER
	cd /media/tfli-0070/DATASAVE1/Data/Data_Brute
	exp=$(ls)
	echo -e "Nom des RUN Sureselect :\n${exp}"
	echo "***********************************"
	echo "Saisir le nom RUN :"
	read NAME_RUN
	# Choix du RUN à analyser
	echo "********"
	RUN=/media/tfli-0070/DATASAVE1/Data/Data_Brute/$NAME_RUN
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
	# Choix de l'Analyse du Pipeline
	CHOIX_ANALYSE
	# Choix de la portée de l'analyse
    echo "***********************************"
    echo "Voulez vous lancer l'analyse sur tous les patients ou sur un patient?"
    echo -e "Choix 1 - All (Analyse demandé Sur tous les patients du RUN sélectioné\n" 
    echo -e "Choix 2 - Unique (Une analyse sur 1 patient \nou une selection de plusieurs patients délimité par des,)\n Choix des patients dessous\n"
    echo -e "Rentrez le nom du choix (Exemple : All)"
    echo "***********************************"
    read SELECTION
}

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
	elif [ "$ANALYSE" = "All" ] ; then
		python3 $DICTIONNARY -o $DICT -s True -outstat $STAT_DICT
		VERIFY_FILE $STAT_DICT
	else
		echo "Aucune preanalyse a effectué pour : $ANALYSE" >> $LOG
	fi
	# *****************************
	# Lancement du code 
	# *****************************
	# Si l'utilisateur travaille sur tous les patients du RUN
	if [ "$SELECTION" = "All" ]
		then

		for name in $(cut -d "," -f1 $TPL | grep -e "-") 
		do
            date >> $LOG
            echo -e "DEBUT ANALYSE NGS COMPLETE"
			echo -e "DEBUT ANALYSE NGS COMPLETE" >> $LOG
            #COMMANDE FOXNGS
            #nextflow run /home/tfli-0070/Bureau/Data/Recherche/Pipeline/SMPHD_Routine/FoxNGS.nf -profile fast -w /media/tfli-0070/DATASAVE2/work --reads ../fastq --results ../results -with-report -with-timeline
            # mv timeline.html ../results/reports
            # mv report.html ../results/reports
            # multiqc ../results -o ../results/reports
            echo -e "FIN ANALYSE NGS COMPLETE"
            echo -e "DEBUT ANALYSE NGS COMPLETE" >> $LOG
			date >> $LOG
		done
		echo "Vous avez analysé $i patient pour l'analyse ${ANALYSE}:" >> $LOG	
		echo -e "**********************************************************************\n" >> $LOG
	# *****************************
	# Si l'utilisateur travaille sur un patient ou une liste de patients en particulier
	elif [ "$SELECTION" = "Unique" ]; then
		# Choix des identifiants du patients
		CHOIX_IDENTIFIANT
		while test "$patient" != "STOP" ; do
			listepatient=$(echo $patient | tr "," "\n")

			for namer in $listepatient;
			
			do
				name=$(cut -d "," -f1 $TPL | grep -e "^${namer}");
				echo "ANALYSE ${name}:";
				echo "ANALYSE ${name}:" >> $LOG;

				echo -e "************\n FIN ANALYSE\n************"
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

function MAIN () {
#    mkdir ../fastq
#    mkdir ../results ../results/reports
    INTERFACE
    LANCEMENT_ANALYSE_PATIENT
#    nextflow clean -f;
#    awk -i inplace '{gsub(/\\n/,"\n")}1' ../results/Statistic_couverture.csv; sed -i /^$/d ../results/Statistic_couverture.csv
}

MAIN
exit 0