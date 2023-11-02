function LOG_AGENT () {  # classic logger
	local prefix="[Pipeline $0 $(date +%Y/%m/%d\ %H:%M:%S)]: par $USER " 
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
	IFS=', ' read -r -a patients_analyses <<< $patient
	echo -e "************"
}


function CHOIX_ANALYSE () {
	echo "***********************************"
	echo -e "Choix 1 : Demultiplexage (1ère étape : Preparation des fichiers FASTQ par patient)\n"
	echo -e "Choix 2 : Analyse (2ème étape : Lancement de l'analyse par patient \n du fichier brute à la création du rapport (Combinaison du choix 2 à 4))\n"
	echo -e "Choix 3 : All (Demultiplexage + Analyse)\n"
	echo -e "Choix 0 - Retour"
	echo "***********************************"
	echo "Entrez le nom du choix. (Exemple All)  "
	echo "Quelle analyse souhaitez vous effectuer ? "
	echo "********"
	read ANALYSE
	echo "********"
}


function CREATION_RECHERCHE_FILE () {
	ANALYSIS_FOLDER=${RUN}/${SORTIE}
	create_new_folders=$1

	LOG=$ANALYSIS_FOLDER/log_lancement_FoxNGS_Routine_v1_2_b.txt
	# Recherche automatique des fichiers brutes
	# Recuperation des data ngs brutes
	foldernamebcl=$(ls $RUN | grep ??????_*_????_*)
	# Chemin vers le fichier brute de séquençage bcl
	TPL=$( ls ${ANALYSIS_FOLDER}/*.csv | head -n 1) 
	BCL=$RUN/$foldernamebcl

	if [ "$create_new_folders" -eq "1" ];
	then
		mkdir $ANALYSIS_FOLDER
    	mkdir $ANALYSIS_FOLDER/fastq/
		mkdir $ANALYSIS_FOLDER/results/
		mkdir $ANALYSIS_FOLDER/results/reports/

		TPL=$BCL/*.csv

	elif [ "$create_new_folders" -eq "2" ];
	then
		rm -r $ANALYSIS_FOLDER/fastq/*
		rm -r $ANALYSIS_FOLDER/results/*

		TPL=$ANALYSIS_FOLDER/*.csv
	fi
	# Récupération du fichier du RUN contenant les échantillons
	# Note : Ce fichier doit contenir le nom du run
	# test -f $TPL
}


function VERIFY_FILE () {
	file=$1
	# Fonction qui verifie si le fichier existe

	# Si le fichier n'exite pa ou il est vide : Creation du dictionnaire
	if test -f $file ;
		then
		echo 1
	# S'il n'existe pas arret du programme
	else
		echo 0
	fi
}


function CHECK_FASTQ () {
	path=$1
	samplesheet=$2
	return=1

	for name in $(cut -d "," -f1 $samplesheet | grep -e "-");
		do
		for R in 1 2;
			do
				if test -f ${path}/${name}_S*_R${R}_001.fastq.gz;
				then
					return=0
					((R++))
				else
					return=1
					break 2
				fi
			done
		done
}


function NETTOYAGE () {
	echo -e "Le nottoyage des processus supprimerait $(nextflow clean -n | wc -l) files\n"
	echo -e "Souhaitez-vous effectuer le nettoyage des fichiers intermédiaires ?\n"
	echo -e "Choix 1 : Oui\n"
	echo -e "Choix 0 : Non\n" 
	echo "********"
	read nettoyage

	if [[ "$nettoyage" -eq "1" ]];
	then
		/home/tfli-0070/BioNGSTools/nextflow clean -f;
	fi
}


function PREPARATION_FASTQ () {
	cd $ANALYSIS_FOLDER
	echo -e "**********************************************************************\n" >> $LOG
	echo -e "$(date) : LAUNCHING BCL2FASTQ" >> $LOG
    echo -e "singularity exec -B /media/tfli-0070/DATASAVE1 /home/tfli-0070/BioNGSTools/Singularity_database/bcl2fastq_docker.sif bcl2fastq --runfolder $BCL --output-dir $ANALYSIS_FOLDER/fastq --sample-sheet $TPL --use-bases-mask Y$read,I$,I8,Y$read --no-lane-splitting -r 32 -p 32 -w 32" >> $LOG
	
	read="$(echo "cat /RunParameters/Setup/Read1/text()" | xmllint --nocdata --shell $BCL/RunParameters.xml | sed '1d;$d')";
	index="$(echo "cat /RunParameters/Setup/Index1Read/text()" | xmllint --nocdata --shell $BCL/RunParameters.xml | sed '1d;$d')";

	singularity exec -B /media/tfli-0070/DATASAVE1 /home/tfli-0070/BioNGSTools/Singularity_database/bcl2fastq_docker.sif bcl2fastq --runfolder $BCL/ --output-dir $ANALYSIS_FOLDER/fastq --sample-sheet $TPL --use-bases-mask Y$read,I$index,I$index,Y$read --stats-dir $ANALYSIS_FOLDER/results/reports/ --reports-dir $ANALYSIS_FOLDER/results/reports/ --no-lane-splitting -r 32 -p 32 -w 32

	if [ $? -eq 0 ];
		then
			echo -e "$(date) : BCL2FASTQ COMPLETE" >> $LOG
	else
		echo -e "$(date) : BCL2FASTQ INCOMPLETE\n" >> $LOG
		return
	fi

	cp $TPL $ANALYSIS_FOLDER
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
	if [ ! -d $SORTIE ]
	then
		echo "Creating new results folders"
		CREATION_RECHERCHE_FILE 1
	else
		echo "Results folder found"
		echo "********" 
		echo -e "\n Voulez vous recréer le dossier de résultats ?"
		echo -e "Choix 1 - Oui\n" 
		echo -e "Choix 2 - Non\n"
		echo -e "Choix 0 - Retour\n"
		echo -e "Rentrez le nom du choix (Exemple : 1)"
		echo "***********************************"
		read REMAKE_RESULTS

		if [ "$REMAKE_RESULTS" = "0" ];
			then
			INTERFACE
		elif [ "$REMAKE_RESULTS" = "1" ];
			then
			CREATION_RECHERCHE_FILE 2
		else
			CREATION_RECHERCHE_FILE 0
		fi
	fi

	if [[ $? -ne 0 ]];
	then
		echo -e "$(date) : UNABLE TO FIND SAMPLESHEET"
		echo -e "$(date) : UNABLE TO FIND SAMPLESHEET" >> $LOG
		exit 1
	fi

	echo "********"
	echo "***********************************"
	# Choix du génome de référence
	# Par default choix DEFAULT pour accelerer la procedure de lancement
	# Amélioration 4.0 Choix des genomes de référence (hg19,hg37 et hg38)
	REFERENCE="DEFAULT"
	# Choix de l'Analyse du Pipeline
	CHOIX_ANALYSE

	if [ "$ANALYSE" = "0" ];
		then
		INTERFACE
	fi

	# Choix de la portée de l'analyse
    echo "***********************************"
	listepatient=$(cut -d "," -f1 $TPL | grep -e "[-]") 
	echo -e "Liste des patients dans le $RUN:\n${listepatient}" 
    echo -e "\n Voulez vous lancer l'analyse sur tous les patients ou sur un patient?"
    echo -e "Choix 1 - All (Analyse demandé Sur tous les patients du run\n" 
    echo -e "Choix 2 - Unique (Une analyse sur 1 patient \nou une selection de plusieurs patients délimité par des virgules)\n Choix des patients dessous\n"
    echo -e "Choix 0 - Retour\n"
	echo -e "Rentrez le nom du choix (Exemple : 1)"
    echo "***********************************"
    read SELECTION

	if [ "$SELECTION" = "0" ];
		then
		CHOIX_ANALYSE
	fi

	if [ "$ANALYSE" = "0" ];
		then
		INTERFACE
	fi
}


function LANCEMENT_ANALYSE_PATIENT () {
	# Fonction Demultiplexage only
	if [ "$ANALYSE" = "1" ] || [[ $ANALYSE = "3" && $demultiplexage -eq 0 ]]; then
		PREPARATION_FASTQ
		# Sortie du programme
		demultiplexage=1
	# Recuperation du dernier dictionnaire avant le lancement du RUN
	# Idée partir d'une même pour les patients lancés en même temps
	elif [ "$ANALYSE" = "2" ] ; then
		echo -e "Passage en mode analyse après sélection des patient.\n"
		echo -e "Passage en mode analyse après sélection des patient.\n" >> $LOG	
	else
		echo "Aucune preanalyse a effectué pour : $ANALYSE"
		echo "Aucune preanalyse a effectué pour : $ANALYSE" >> $LOG
	fi
	# *****************************
	# Lancement du code 
	# *****************************
	# Si l'utilisateur travaille sur tous les patients du RUN
	if [ "$SELECTION" = "1" ]
		then
			echo -e "**********************************************************************" >> $LOG
            echo -e "$(date) : STARTING NGS ANALYSIS ON ALL SAMPLES"
			echo -e "$(date) : STARTING NGS ANALYSIS ON ALL SAMPLES" >> $LOG
            #COMMANDE FOXNGS
            echo -e "/home/tfli-0070/BioNGSTools/nextflow run /home/tfli-0070/Bureau/Data/Recherche/Pipeline/SMPHD_Routine/FoxNGS.nf -profile fast -w /media/tfli-0070/DATASAVE2/work --reads $ANALYSIS_FOLDER/fastq --results $ANALYSIS_FOLDER/results -with-report -with-timeline"
			/home/tfli-0070/BioNGSTools/nextflow run /home/tfli-0070/Bureau/Data/Recherche/Pipeline/SMPHD_Routine/FoxNGS.nf -profile urgent -w /media/tfli-0070/DATASAVE2/work --reads $ANALYSIS_FOLDER/fastq --results $ANALYSIS_FOLDER/results -with-report -with-timeline --resume
			
			if [ $? -eq 0 ];
			then
				echo -e "$(date) : NGS ANALYSIS ON ALL SAMPLES COMPLETE" >> $LOG
			else
				echo -e "$(date) : NGS ANALYSIS INCOMPLETE\nCHECK NEXTFLOW REPORT FOR MORE INFORMATION" >> $LOG
				return
			fi

			echo -e "**********************************************************************\n" >> $LOG
	# *****************************
	# Si l'utilisateur travaille sur un patient ou une liste de patients en particulier
	elif [ "$SELECTION" = "2" ]; then
		# Choix des identifiants du patients
		mkdir $ANALYSIS_FOLDER/fastq/tmp
		stop_trigger=0
		patients_analyses=()
		unique_patient=${#patients_analyses[@]}

        listepatient=$(echo $patient | tr "," "\n")

		while [ $stop_trigger -eq 0 ];
		do
			CHOIX_IDENTIFIANT

			for unique_patient in ${patients_analyses[@]};
			do

				if [ "$unique_patient" == "STOP" ];
					then
					stop_trigger=1
					date >> $LOG
					echo -e "$(date) : STARTING NGS ANALYSIS ON A SUBSET OF SAMPLES" >> $LOG
					/home/tfli-0070/BioNGSTools/nextflow run /home/tfli-0070/Bureau/Data/Recherche/Pipeline/SMPHD_Routine/FoxNGS.nf -profile fast -w /media/tfli-0070/DATASAVE2/work --reads $ANALYSIS_FOLDER/fastq/tmp/ --results $ANALYSIS_FOLDER/results -with-report -with-timeline --resume
					rm -rf $ANALYSIS_FOLDER/fastq/tmp/

					if [ $? -eq 0] ;
						then
							echo -e "$(date) : NGS ANALYSIS ON A SUBSET OF SAMPLES COMPLETE" >> $LOG
					else
						echo -e "$(date) : NGS ANALYSIS INCOMPLETE\nCHECK NEXTFLOW REPORT FOR MORE INFORMATION" >> $LOG
						return
					fi

				else
					echo $unique_patient
					echo $unique_patient >> $LOG
					cp $ANALYSIS_FOLDER/fastq/${unique_patient}_S*_R{1,2}_001.fastq.gz $ANALYSIS_FOLDER/fastq/tmp/
				fi
			done
		done

		echo -e "***********\nFIN ANALYSE\n***********"
	# *****************************
	# Error of variable selection
	else
		echo "ERREUR de saisie lors du programme pour la variable ${SELECTION}: $0 Arret du programme" >> $LOG
		date >> $LOG
		exit 1
	fi
	echo -e "***********\nEnd of program\n***********"
}


function RAPPORTS_QUALITE () {
	mv timeline*.html $ANALYSIS_FOLDER/results/reports/
    mv report*.html $ANALYSIS_FOLDER/results/reports/
	multiqc $ANALYSIS_FOLDER/results -o $ANALYSIS_FOLDER/results/reports/
	awk -i inplace '{gsub(/\\n/,"\n")}1' $ANALYSIS_FOLDER/results/Statistic_couverture.csv; sed -i /^$/d $ANALYSIS_FOLDER/results/Statistic_couverture.csv
}


function MAIN () {
	demultiplexage=0
    INTERFACE
	LOG_AGENT
    LANCEMENT_ANALYSE_PATIENT
	RAPPORTS_QUALITE
	NETTOYAGE
}


MAIN