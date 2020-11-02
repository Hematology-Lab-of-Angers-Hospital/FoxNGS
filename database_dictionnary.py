# -*- coding: utf-8 -*-
#!/usr/bin/python3
# Auteur Maxime RENARD
# Date:28/08/19
# Laboratoire d'Hematologie
# ******************************************
# Creation of Database of variant and implementation of database 
# Write Database in result file
# ******************************************
# Import bibliothèque
import sys
import pandas as pd
import argparse
import os
import re
import sys
import json

# ******************************************
# Function

# Initialisation du dictionnaire
def initialize_dico():
	"Initialise database and rename si un changement notable dans le pipeline est réalisé"
	dico_database = {}
	dico_database["Echantillon_patient"] = []
	return(dico_database)

# Ecriture du dictionnaire
def write_dictionnary(write_dico,name):
	"""
		Write Dictionnary in jsonformat
	""" 
	with open(name, 'w') as outfile:
		json.dump(write_dico, outfile)

# Chargement du dictionnaire
def load_dictionnary(name):
	"""
	load dictionnary 
	"""
	with open(name,"r") as infile:
		dico_database = json.load(infile)
	return(dico_database)

# Insertion des valeurs dans le dictionnaire
def Insertion_dictionnaire(dico_database,name_file):
	"""
		IN : Database
		OUT: Dictionnaire with insertion of data
	"""

	# Preparation des elements pour bien ordonner le dictionnaire
	# Partie "Patient"
	liste_VAF = ["Mean_VAF","Method_Variant_GATK","Total_DP_GATK","VAF_GATK","Method_Variant_Mutect2","Total_DP_Mutect2","VAF_Mutect2","Method_Variant_Varscan","Total_DP_Varscan","VAF_Varscan","Method_Variant_Pindel","Total_DP_Pindel","VAF_Pindel"]
	echantillon = name_file.split(".")[0].split("_")[3]

	# Declaration du nom de l'echantillon
	# Non ecriture de TEM HORIZON dans le dictionnaire == Passage de l'insertion
	liste_exclure = ["TEM-HO","EEQ"]


	if echantillon.startswith(liste_exclure[0]) or echantillon.startswith(liste_exclure[1]):
		print('No match echantillon of control',echantillon)
		sys.exit(0)
	else:
		print("Nom de l echantillon en cours d insertion du dictionnaire : " , echantillon)


	liste_patient = ['0409-GJ-09112010','0422-RC-19072006',"0013-SJ-30092010","0025-CR-13012010","0034-GA-29102009",
					"0083-WP-07012011","0161-WF-07102011","0165-ND-30062011","0224-BC-10062009","0737-KJ-16112011",
					"0113-GM-17092013","0117-GP-25012010","0294-CR-14112008","0041-LB-21032011","0628-BJ-07062016",
					"0859-DG-16042015"]
	if echantillon in liste_patient:
		print("Patient à repasser")
		sys.exit(0)
	# Ecriture dans le dictionnaire
	if echantillon not in dico_database["Echantillon_patient"]:
		dico_database["Echantillon_patient"].append(echantillon)
	# Lecture du fichier d'entrée d'annotation 
	with open(name_file,"r+") as file_annotation:
		# Action sur les lignes
		for ligne in file_annotation:
			# Si la ligne ne correspond pas au header
			if ligne[0:2] != "ID":
				# Récupération de ID du variant
				info=ligne.split(";")
				ID= info[0]
				# Nouveau variant identifié qui n'est pas un témoin horizon
				if ID not in dico_database:
					# Création d'un sous dictionnaire
					dico_database[ID] = {}
					# Sous dictionnaire pour un patient
					dico_database[ID][echantillon] = {}
					# Stockage de la liste des échantillons ayant cet ID
					dico_database[ID]["Liste_patient_ID"] = []
					dico_database[ID]["Liste_patient_ID"].append(echantillon)
					# Stockage de l'annotation information générale
					dico_database[ID]["Annotation"] = {}
					# Ecriture dans le dictionnaire
					i = 1
					while i < len(header):	
						nom_attribut = header[i]
						value = info[i]
						# Element ID à la base de ce sous dictionnaire
						if nom_attribut == "ID":
							continue
						# Insertion des valeurs de VAF
						elif nom_attribut in liste_VAF:
							dico_database[ID][echantillon][nom_attribut] = value
						# Insertion des valeurs d'Annotation
						else:	
							dico_database[ID]["Annotation"][nom_attribut] = value
						i+=1
				
				# Implementation de valeur d'un patient pour un ID reconnu
				elif ID in dico_database and echantillon not in dico_database[ID]["Liste_patient_ID"] :
					dico_database[ID][echantillon]= {}
					dico_database[ID]["Liste_patient_ID"].append(echantillon)
					i = 1
					# Ecriture seulement des VAF pour chaque patient
					while i < len(header):
						nom_attribut = header[i]
						value = info[i]
						# Insertion des valeurs de VAF pour le patient
						if nom_attribut in liste_VAF:
							dico_database[ID][echantillon][nom_attribut] = value
						i+=1
				
				# ID déja connu pour ce patient: Signalement
				elif ID in dico_database and echantillon in dico_database[ID]["Liste_patient_ID"] :
					print("Ce patient {:20s} et cet ID {:20s} a déja été stocké.".format(echantillon,ID))

				# Other error to understand
				else :
					print("Other error")
					sys.exit(1)
			# Recuperation du header du tableau
			else:
				header = ligne.split(";")

	return(dico_database)

# Statistic sur le dictionnaire
def statistic(dico_database,out):
	# filout CSV of statistic
	# Mise en place Tabulation
	# liste patient to remove
	"""
		IN: Dictionnaire
		OUT: Ecriture des informations principales dans un fichier
	"""
	liste_patient = ['0409-GJ-09112010','0422-RC-19072006',"0013-SJ-30092010","0025-CR-13012010","0034-GA-29102009",
					"0083-WP-07012011","0161-WF-07102011","0165-ND-30062011","0224-BC-10062009","0737-KJ-16112011",
					"0113-GM-17092013","0117-GP-25012010","0294-CR-14112008","0041-LB-21032011","0628-BJ-07062016",
					"0859-DG-16042015"]
	tab = "\t" 
	with open(out,"w") as stat_out:
		
		print(" Patient",dico_database["Echantillon_patient"])

		count_total = len(dico_database["Echantillon_patient"])
		print("Number Patient",count_total) 
		# Mutation par gene 
		stat_out.write("ID\tGene\tIGV\tcytoBand\tFreq\tPatient\tRedondance\tExonicFunc\t Notation\n")

		# variable temporaire pour voir s'il y a aucun doublon d'echantillon
		tmp = 0
		# suppression de patient à ne pas prendre en compte
		dico_database["Echantillon_patient"] = [i for i in dico_database["Echantillon_patient"] if i not in liste_patient]
	 
		# Search in dictionnary
		for cle in sorted(dico_database.keys()):
			
			# Suppression des fichiers doublons
			if tmp == cle:
				print("Error of duplication ID")
				sys.exit(1)
			# Si c'est bien un variant	
			if cle.startswith("chr") == True:
				tmp = cle
				# Afficher l'information sur un variant
				#if cle == "chr17_74732935_CGGCGGCTGTGGTGTGAGTCCGGGG_C":
				#	print(dico_database[cle])
				
				# Si la cle n'est pas vide
				if len(dico_database[cle]["Liste_patient_ID"]) > 0:
					# Suppresion ID du patient
					dico_database[cle]["Liste_patient_ID"] = [i for i in dico_database[cle]["Liste_patient_ID"] if i not in liste_patient]
					for echantillon in liste_patient: 
						if echantillon in dico_database[cle]:
							# Suppression des echantillons
							del dico_database[cle][echantillon]	
						

					# Operation in dictionnary
					countID = len(dico_database[cle]["Liste_patient_ID"])
					freq_ID = round(float(len(dico_database[cle]["Liste_patient_ID"]) / count_total),4)
					if freq_ID > 0.90:
						Note_freq = "Always-Artefact"
					elif freq_ID > 0.70 and freq_ID <= 0.90:
						Note_freq = "Often"
					elif freq_ID <= 0.70 and freq_ID > 0.40:
						Note_freq = "Usually"
					elif freq_ID <= 0.40 and freq_ID >= 0.10:
						Note_freq = "Sometimes"
					else:
						Note_freq = "Rarely"
				# Détermination des patients ayant cet ID
					line_patient = ""
					for patient in dico_database[cle]["Liste_patient_ID"]:
						if line_patient != "":
							line_patient = line_patient + "," + patient
						# First line
						else:
							line_patient = line_patient + patient
					line = cle + tab + str(dico_database[cle]["Annotation"]["Gene.refGene"]) + tab + str(dico_database[cle]["Annotation"]["IGV"]) + tab \
					+ str(dico_database[cle]["Annotation"]["cytoBand"]) + tab  + str(freq_ID) + tab + line_patient + tab \
					+ Note_freq + tab + str(dico_database[cle]["Annotation"]["ExonicFunc.refGene"]) + tab + str(dico_database[cle]["Annotation"]["AAChange.refGene"]) + "\n" 
					stat_out.write(line)
					
					
	return(dico_database)

# *****************************************************
# ********************* Main **************************
# *****************************************************
if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog="database_dictionnary_v2")
	# Create sub_parsers (one for Treatment_Annotation, other for Combine)
	
	parser.add_argument("-d", "--directory",
		dest="directory",
		required=False,
		type=str,
		help="Directory of result annotation"
	)
	parser.add_argument("-f", "--fileresult",
		dest="file",
		required=False,
		type=str,
		help="File of result annotation for annovar fusion"
	)
	parser.add_argument("-o", "--outout",
		dest="out",
		required=True,
		type=str,
		help="Name of new dictionnary and his location"
	)
	parser.add_argument("-c", "--newdictionnary",
		dest="dictionnary",
		required=False,
		type=bool,
		help="Argument to indicate creation of new dictionnary -c TRUE sinon ne pas afficher cette option"
	)
	parser.add_argument("-s", "--statistic",
		dest="stat",
		required=False,
		type=bool,
		help="Argument to activate statistic option"
	)
	parser.add_argument("-outstat", "--outstatistic",
		dest="outstat",
		required=False,
		type=str,
		help="Argument to indicate fileout of statistic csv"
	)
	args = parser.parse_args()
	# Si args.dictionnary == TRUE : Initialisation du dictionnnaire
	if args.dictionnary:
		
		#initialisation dictionnaire si un changement notable dans le pipeline est réalisé
		dico_database = initialize_dico()
	# Sinon chargement du dictionnaire	
	else:
		# Load dictionnary
		dico_database = load_dictionnary(args.out)
		
	# if args.stat == True Statistic on dictionnary
	if args.stat:
		# Appel de la fonction statistic
		dico = statistic(dico_database,args.outstat)
		# Ecriture des changements suite à d'éventuel suppression
		write_dictionnary(dico,args.out)
	# Insertion value in dictionnary	
	else:
		# Deplacement dans le repertoire
		os.chdir(args.directory)
		dico = Insertion_dictionnaire(dico_database,args.file)
		# Save new dictionnary
		write_dictionnary(dico,args.out)
		print("L'insertion des données {:20s} dans la database s'est bien déroulée.".format(args.file))
	
	
	sys.exit(0)
#python3 /home/t-chu-027/Bureau/Recherche/Pipeline/SMPHD_v2.8/database_dictionnary_v2.8.py -o /media/t-chu-027/Elements/Result_NGS/Dictionnary_database/Project_FIM_Database_variant.json -s True -outstat /media/t-chu-027/Elements/Result_NGS/Dictionnary_database/Project_FIM_Database_variant.csv
