# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# Auteur Maxime RENARD
# Date:03/07/20
# Laboratoire d'Hematologie
# Ce code traite les fichiers d'annotations sorties par les différents logiciels d'appel de variant par le pipeline SMPHD:
# 1. Calcul du Variant Allele Frequency 
# 2- Mise en place des filtres sur les fichiers d'annotations
# 3- Merging de différents fichiers d'annotation provenant de variants différents (GATK, Mutect, Varscan, Pindel)
# 4- Fichier préparé pour l'insertion futur dans un dictionnaire exhaustif
# 5- Mise en place d"un dernier filtre sur GnomAD pour la lecture facilité par les biologistes
# 6 - Vu dans la database et mark SNP

# Import bibliothèque
import pandas as pd
import argparse
import os
import sys
import datetime
import json
import re
# Graphique
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles


# *******************************************
#  Function
# *******************************************

# Calcul des VAF et creation d'un identifiant de variation
def preparation_file(file,output,method):
	"""
		1- Creation d'un identifiant de ligne pour aider au merging
		2- Calcul du VAF
		
	"""
	data = pd.read_csv(file, sep='\t', index_col=False)

	# IARC
	data.loc[:,"IARC"] = data.loc[:,"IARC"].replace(to_replace=r"x3b", value=";", regex=True)
	data.loc[:,"IARC"] = data.loc[:,"IARC"].replace(to_replace=r"\\", value="", regex=True)

	# Cosmic 90
	data.loc[:,"cosmic91"] = data.loc[:,"cosmic91"].replace(to_replace=r"x3b", value=";", regex=True)
	data.loc[:,"cosmic91"] = data.loc[:,"cosmic91"].replace(to_replace=r"x3d", value="", regex=True)
	data.loc[:,"cosmic91"] = data.loc[:,"cosmic91"].replace(to_replace=r"\\", value="", regex=True)
	data.loc[:,"cosmic91"] = data.loc[:,"cosmic91"].replace(to_replace=":", value="", regex=True)
	name_Method_variant = "Method_Variant"
	# Cosmic 89
	data.loc[:,"cosmic89"] = data.loc[:,"cosmic89"].replace(to_replace=r"x3b", value=";", regex=True)
	data.loc[:,"cosmic89"] = data.loc[:,"cosmic89"].replace(to_replace=r"x3d", value="", regex=True)
	data.loc[:,"cosmic89"] = data.loc[:,"cosmic89"].replace(to_replace=r"\\", value="", regex=True)
	data.loc[:,"cosmic89"] = data.loc[:,"cosmic89"].replace(to_replace=":", value="", regex=True)
	name_VAF= "VAF"

	i=0
	if data.loc[i,'CHROM'] != '.':
		while i <= data.shape[0] -1 :
			pass
			# colonne ID
			data.loc[i,"ID"] = data.loc[i,'CHROM']+ "_" + str(data.loc[i,'POS']) + "_" + str(data.loc[i,'REF']) + "_" + str(data.loc[i,'ALT'])
			# colonne method
			
			data.loc[i,name_Method_variant] = method
			#Colonne IGV 
			data.loc[i,"IGV"] = data.loc[i,'CHROM']+ ":" + str(data.loc[i,'POS'])

			# Cosmic annotation
			# Cosmic 90
			# If exist
			if data.loc[i,"cosmic91"] != ".":
				
				ID=data.loc[i,"cosmic91"].split(';')[0]
				OCC=data.loc[i,"cosmic91"].split(';')[1]
				
				data.loc[i,"ID_cosmic91"] = ID
				data.loc[i,"OCC_cosmic91"] = OCC
			# elif column vide
			else:
				data.loc[i,"ID_cosmic91"] = "."
				data.loc[i,"OCC_cosmic91"] = "."
			# Cosmic 89

			if data.loc[i,"cosmic89"] != ".":
				
				ID=data.loc[i,"cosmic89"].split(';')[0]
				OCC=data.loc[i,"cosmic89"].split(';')[1]
				
				data.loc[i,"ID_cosmic89"] = ID
				data.loc[i,"OCC_cosmic89"] = OCC
			# elif column vide
			else:
				data.loc[i,"ID_cosmic89"] = "."
				data.loc[i,"OCC_cosmic89"] = "."
			# ******

			# IARC Annotation
			if data.loc[i,"IARC"] != ".":
				
				ID=data.loc[i,"IARC"].split(';')[0]
				OCC=data.loc[i,"IARC"].split(';')[1]
				
				data.loc[i,"Genomic_IARC"] = ID
				data.loc[i,"Transactivation_Structure_Function_IARC"] = OCC
			# elif column vide
			else:
				data.loc[i,"Genomic_IARC"] = "."
				data.loc[i,"Transactivation_Structure_Function_IARC"] = "."


			# Calcul VAF
			if method == "GATK" or method == "Mutect2":
				
				ALT = data.loc[i,"C5-${name}:AD"].split(',')[1]
				REF = data.loc[i,"C5-${name}:AD"].split(',')[0]
				data.loc[i,"GATK_ALT"] = float(ALT)
				data.loc[i,"GATK_REF"] = float(REF)
				data.loc[i,name_VAF] = round(float(data.loc[i,"GATK_ALT"])/data.loc[i,"C5-${name}:DP"],4)
					
			elif method == "Varscan":
				
				data.loc[i,name_VAF] = round(float(data.loc[i,"Sample1:AD"])/data.loc[i,"Sample1:DP"],4)
			
			elif method == "Pindel":
				# All read
				All_ALT = data.loc[i,"All_read:AD"].split(',')[1]
				All_REF = data.loc[i,"All_read:AD"].split(',')[0]
				data.loc[i,"All_ALT"] = float(All_ALT)
				data.loc[i,"All_REF"] = float(All_REF)
				DP_All = data.loc[i,"All_ALT"] + data.loc[i,"All_REF"]
				data.loc[i,name_VAF] = round(float(data.loc[i,"All_ALT"]/(DP_All)),4)

				# Dupmark
				Dup_ALT = data.loc[i,"All_read:AD"].split(',')[1]
				Dup_REF = data.loc[i,"All_read:AD"].split(',')[0]
				data.loc[i,"Dup_ALT"] = float(Dup_ALT)
				data.loc[i,"Dup_REF"] = float(Dup_REF)
				data.loc[i,"Total_DP"] = data.loc[i,"Dup_ALT"] + data.loc[i,"Dup_REF"]
				data.loc[i,name_VAF + "Dup"] = round(float(data.loc[i,"Dup_ALT"] / data.loc[i,"Total_DP"] ),4)

			else: 
				print("FLT3 or error")
				sys.exit(1)
			# Cosmic 90
			data.loc[:,"ID_cosmic91"] = data.loc[:,"ID_cosmic91"].replace(to_replace="ID", value="", regex=True)
			data.loc[:,"OCC_cosmic91"] = data.loc[:,"OCC_cosmic91"].replace(to_replace="OCCURENCE", value="", regex=True)
			# Cosmic 89
			data.loc[:,"ID_cosmic89"] = data.loc[:,"ID_cosmic89"].replace(to_replace="ID", value="", regex=True)
			data.loc[:,"OCC_cosmic89"] = data.loc[:,"OCC_cosmic89"].replace(to_replace="OCCURENCE", value="", regex=True)

			data.loc[:,"ANNOVAR_DATE"] = data.loc[:,"ANNOVAR_DATE"].replace(to_replace="\n", value="", regex=True)
			# Compteur
			i+=1
	# If file of annotation is empty
	else: 
		print("Fichier vide")
		filter_out ='Filter_' + output
		data.to_csv(filter_out, sep = ';')
		return("error")
	
	name = 'ALL_VAF_' + output
	# Write file
	data.to_csv(name, sep = ';')
	return(data)

# Function Delete columns and clean table First filter to suppress column
def action_column(data,method):
	"""
		Input : Table annotation with column redundance and reorder column
		Output: Table annotation without column redundance and reorder column
	"""
	# ****************
	# First suppresion column
	commun_info=['Xref.refGene','GeneDetail.refGene','ALLELE_END','cosmic91','IARC']
	data.drop(commun_info, axis='columns', inplace=True)
	name_DP = "Total_DP" 
	name_VAF = "VAF"
	# ,'DS'
	if method == "GATK":
		data.rename(columns={'C5-${name}:DP': name_DP}, inplace=True)
		info_delete=['AC','AF','AN','BaseQRankSum','ExcessHet','FS','InbreedingCoeff',
					'MLEAC','MLEAF','MQRankSum','QD','ReadPosRankSum','SOR','C5-${name}:AD',
					'C5-${name}:PL','GATK_REF','GATK_ALT','DP','C5-${name}:GT','QUAL','MQ','C5-${name}:GQ']
		
	elif method == "Mutect2":
		
		data = data.rename(columns={'C5-${name}:DP':name_DP})
		info_delete=['QUAL','ECNT','MBQ','MFRL','MMQ','MPOS','NCount','PON','POPAF','RPA','RU',
					'SEQQ','STR','STRANDQ','STRQ','TLOD','UNIQ_ALT_READ_COUNT','CONTQ','GERMQ','NALOD'
					,'NLOD','OCM','ROQ','C5-${name}:PID','C5-${name}:PL','C5-${name}:AD','C5-${name}:GT',
					'GATK_REF','GATK_ALT','DP','C5-${name}:AF','C5-${name}:F1R2','C5-${name}:F2R1','C5-${name}:GQ','C5-${name}:PGT','C5-${name}:PS','C5-${name}:SB','AF']
		

	elif method == "Varscan":
		

		data = data.rename(columns={'Sample1:DP':name_DP})
		info_delete=['QUAL','ADP','NC','HET','HOM','WT','Sample1:SDP','Sample1:AD','Sample1:RD',
					'Sample1:FREQ','Sample1:RBQ','Sample1:ABQ','Sample1:RDF','Sample1:RDR',
					'Sample1:ADF','Sample1:ADR','Sample1:GT','Sample1:GQ','Sample1:PVAL','AF']
	
	elif method == "Pindel":
		diff_Pindel = ["END","HOMLEN","HOMSEQ","SVLEN","SVTYPE","NTLEN"]

		info_delete=["PF","All_read:PL","All_read:GT","All_read:RD","All_read:AD","Duplicate_mark:PL","Duplicate_mark:GT","Duplicate_mark:RD","Duplicate_mark:AD"]

	else:
		print ("error of choice methode")
		sys.exit(1)
	# delete	
	data.drop(info_delete, axis='columns', inplace=True)

	# ****************
	#First Reorder column
	cols = list(data.columns.values)
	reorder=[]
	for c, value in enumerate(cols, 0):

		# Si ANNOVAR_DATE ajout des informations VAF et Method
		if value == "ANNOVAR_DATE":
			reorder.append("Method_Variant")
			reorder.append(name_DP)
			reorder.append(name_VAF)

		# Si value = cytoband reinsertion des colonnes Occurences
		elif value == "cytoBand":
			reorder.append(value)
			reorder.append("IGV")
			# Cosmic 90
			reorder.append("ID_cosmic91")
			reorder.append("OCC_cosmic91")

			# Cosmic 89
			reorder.append("ID_cosmic89")
			reorder.append("OCC_cosmic89")
			# IARC
			reorder.append("Genomic_IARC")
			reorder.append("Transactivation_Structure_Function_IARC")
			
		# Fin des colonnes
		elif value == "Method_Variant":

			reorder.append("ANNOVAR_DATE")
			break

		# Passage de la colonne name_DP pour ne pas la noter en double	
		elif value == name_DP:
			continue
		# Condition normal insertion des colonnes dans le meme ordre
		else:
			reorder.append(value)
	#reindexation des colonnes
	reindex_data = data.reindex(columns =reorder)

	return reindex_data

# Réalisation des filtres simples sur les variants avant l'entrée dans le dictionnaire
def filter_annotation_dico(annotation,out,method):
	"""
		Input: Fichier combiné
		Output: Fichier Filtré sur les variants d'interet 
	"""
	# 
	annotation.loc[:,"Func.refGene"] = annotation.loc[:,"Func.refGene"].replace(to_replace=r"x3b", value="", regex=True)
	annotation.loc[:,"Gene.refGene"] = annotation.loc[:,"Gene.refGene"].replace(to_replace=r"x3b", value="", regex=True)
	# Keep "ncRNA_exonic" and UTR aund upstream
	Function_remove = ["intronic","UTR5","UTR3","ncRNA_intronic","upstream","upstream\downstream","downstream","upstream\downstream","ncRNA_exonic"]
	# keep splicing in exonic fonction only and upstream\downstream
	Exonic_remove = ["upstream\downstream","synonymous_SNV"]

	# ****************
	#Pindel
	# To keep
	# FLT3 exon 14: chr13:28608219-28608351
	# FLT3 exon 15: chr13:28608024-28608128
	FLT3_type = ["INS","DUP:TANDEM"]
	# CALR exon 9 ch19:13054527-13055303
	CALR_type = ["DEL"]
	#initialization of list of line removed by filter
	remove = []
	i = 0
	name_VAF = "VAF"
	# Filtre des VAF <2% sont eliminées
	Filter_VAF = 0.02
	
	# Filter
	# Pour le dictionnaire ne pas enlever les polymorphismes
	while i <= annotation.shape[0] - 1 :
		pass

		# filter function	
		if annotation.loc[i,"Func.refGene"] in Function_remove:
			remove.append(i)
		# Remove only synonymous_SNV and keep frameshift_deletion,insertion, nonframeshift,nonsynonimous_SNV and stopgain
		elif annotation.loc[i,"ExonicFunc.refGene"] in Exonic_remove:
			remove.append(i)


		# Partie AF
		elif annotation.loc[i,name_VAF] < Filter_VAF  and method != "Pindel":
			remove.append(i)
		
		elif method == "Pindel":
			# Cas chr13:FLT3
			# premier filtre sur le type
			if annotation.loc[i,"Gene.refGene"] == "FLT3" and annotation.loc[i,"SVTYPE"] not in FLT3_type:
				remove.append(i)
			# Deuxieme filtre sur lA VAF
			if annotation.loc[i,"Gene.refGene"] == "FLT3"and annotation.loc[i,name_VAF] < 0.01:
				remove.append(i)
			# Cas chr19 CALR
			# premier filtre sur le type
			elif annotation.loc[i,"Gene.refGene"] == "CALR" and annotation.loc[i,"SVTYPE"] not in CALR_type:
				remove.append(i)	
			# Deuxieme filtre sur lA VAF
			elif annotation.loc[i,"Gene.refGene"] == "CALR" and annotation.loc[i,name_VAF] < Filter_VAF:
				remove.append(i)
			# Passge du filtre
			else:	
				annotation.loc[i,"FILTER"]="PASS"


		# Colonne verification du passage des filtres  
		else:	
			annotation.loc[i,"FILTER"]="PASS"		
		i+=1

	# Remove line
	final_list=sorted(list(set(remove))) 
	print("Sur {:4d} lignes au total, il y a  {:4d} lignes eliminé par les filtres donc {:4d} lignes sauvegardées".format(annotation.shape[0],len(final_list),annotation.shape[0] - len(final_list)))
	annotation.drop(final_list, inplace=True)
	# Creation du nouvelle index ID 
	annotation.set_index('ID',inplace=True)
	# Call function Clean columns
	annotation_filter = action_column(annotation,method)
	# enregistrement 
	filter_out ='Filter_' + out
	annotation_filter.to_csv(filter_out, sep = ';')
	



# **************************************************************
# Partie Concatenation des fichiers et figures

# Figure sur les fichiers d'annotation 
def Figure (all_file,output):

	# Concatenate
	concatenate = pd.concat(all_file, axis=0,join='outer',sort=False,keys=["GATK","Mutect2","Varscan","Pindel"])
	GATK = concatenate.loc["GATK"].index
	Mutect = concatenate.loc["Mutect2"].index
	Varscan = concatenate.loc["Varscan"].index
	repartition_variant = concatenate['Method_Variant'].value_counts()
	
	set_GATK = set(GATK) 
	
	set_Mutect = set(Mutect)
	
	set_Varscan = set(Varscan) 
	
	Liste = list(set_GATK) + list (set_Mutect) + list(set_Varscan) 

	unique_liste = sorted(list(set(Liste)))
	print("unique_list",len(unique_liste))
	v = venn3(subsets=[set_GATK, set_Mutect, set_Varscan],set_labels=('GATK', 'Mutect', 'Varscan'))
	c = venn3_circles(subsets=(set_GATK, set_Mutect, set_Varscan),linestyle='dashed', linewidth=0.5)
	name = output + ".png"
	plt.savefig(name)
	# GATK en rouge, Mutect en vert, varscan en bleu
	# Statistique sur les groupes d'ID 
	linebygroup = concatenate.groupby('ID').size()
	#print(linebygroup)

# Fusion des fichiers en un seul avec les differents appel de variant
def Fusion_file(data1,data2,data3,data4,all_file,output):
	"""
		Input: Fichiers d'annotation des 4 appels de variants
		Output: Fichier (Fusion) Fusion des Fichiers et calcul 
	"""
	# Preparation du fichier ne contenant que les annotations
	concatenate = pd.concat(all_file, axis=0,join='outer',sort=False)
	liste = ["Method_Variant","Total_DP","VAF"]
	# suppression des colonnes
	concatenate.drop(liste, axis='columns', inplace=True)
	# Suppression des doublons
	concatenate_one_line = concatenate.loc[~concatenate.index.duplicated(keep='first')]

	# Concatenation des fichiers VAF
	merge_VAF1 = pd.merge(left=data1[liste],right=data2[liste], left_on ='ID', right_on ='ID',copy=False, how="outer",suffixes=('_GATK','_Mutect2'))
	merge_VAF2  = pd.merge(left=data3[liste],right=data4[liste], left_on ='ID', right_on ='ID',copy=False,how="outer",suffixes=('_Varscan','_Pindel'))
	All_merge = pd.merge(left=merge_VAF1,right=merge_VAF2, left_on ='ID', right_on='ID',copy =False,how="outer")
	
	# Concatenation du fichier d'annotation avec le fichier VAF
	combine_file = pd.concat([All_merge,concatenate_one_line], axis=1,join='inner',sort=False)
	# ************************************

	# Creation New column
	# Column of validation for engineer
	combine_file["Validation_BIO_ING"] = "."
	# Column for class of Mutation
	combine_file["Classification_mutation"] = "."
	# Calcul de statistique
	Col_DP = ['Total_DP_GATK','Total_DP_Mutect2','Total_DP_Varscan','Total_DP_Pindel']
	col_VAF = ['VAF_GATK','VAF_Mutect2','VAF_Varscan','VAF_Pindel']
	# Calcul VAF et DP mean
	combine_file['Mean_VAF'] = combine_file[col_VAF].mean(skipna=True,axis=1)
	# Calcul Mean DP
	combine_file['Mean_DP'] = combine_file[Col_DP].mean(skipna=True,axis=1)
	combine_file['Mean_DP'] = combine_file['Mean_DP'].round(3)
	combine_file['Mean_VAF'] = combine_file['Mean_VAF'].round(3)	
	# call column transcript

	# Reorder column
	cols = list(combine_file.columns.values)
	reorder=[]
	liste_VAF = ["Mean_VAF","Mean_DP","Method_Variant_GATK","Total_DP_GATK","VAF_GATK","Method_Variant_Mutect2","Total_DP_Mutect2","VAF_Mutect2","Method_Variant_Varscan","Total_DP_Varscan","VAF_Varscan","Method_Variant_Pindel","Total_DP_Pindel","VAF_Pindel"]
	liste_annot = ["CHROM","Validation_BIO_ING","Classification_mutation","cytoBand","IGV","POS","REF","ALT","Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene"]
	liste_methode = ["Method_Variant_GATK","Method_Variant_Mutect2","Method_Variant_Varscan","Method_Variant_Pindel"]
	liste_annotation = ["ID_cosmic91","OCC_cosmic91","ID_cosmic89","OCC_cosmic89","Genomic_IARC","Transactivation_Structure_Function_IARC","AF_popmax","AF_male","AF_female","AF_raw","AF_afr","AF_sas","AF_amr","AF_eas","AF_nfe","AF_fin","AF_asj","AF_oth","non_topmed_AF_popmax","non_neuro_AF_popmax","non_cancer_AF_popmax","controls_AF_popmax","CLNALLELEID_2020","CLNDN_2020","CLNDISDB_2020","CLNREVSTAT_2020","CLNSIG_2020","SIFT_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","MutationTaster_pred","FATHMM_pred","M-CAP_pred","CADD_phred","Interpro_domain","ICGC_Id","ICGC_Occurrence"]
	for i in liste_annot:
		reorder.append(i)
	# Ajout des VAF seulement et pas de Methode variant
	for el in liste_VAF:
		if el not in liste_methode:
			reorder.append(el)
	#Ajout des annotations
	for annot in liste_annotation:
		reorder.append(annot)
	reorder.append("FILTER")
	reorder.append("ANNOVAR_DATE")
	#reindexation des colonnes
	combine_file_reorder = combine_file.reindex(columns =reorder)
	#transformation des Nan en .
	fusion_file = combine_file_reorder.fillna(".")
	# Ecriture du fichier
	fusion_file.to_csv(output + ".csv", sep = ';')
	return(fusion_file)

# Concatenation of file
def Combine_figure(files,output):
	"""
	Input 1: All file separate by ,
	Input 2: output file concatenate all files and make venn diagrame in png format
	"""
	file1 = files.split(",")[0]
	file2 = files.split(",")[1]
	file3 = files.split(",")[2]
	file4 = files.split(",")[3]
	# GATK
	data1 = pd.read_csv(file1, sep=';', index_col='ID')
	#Mutect
	data2 = pd.read_csv(file2, sep=';', index_col='ID')

	#Varscan
	data3 = pd.read_csv(file3, sep=';', index_col='ID')

	# Pindel
	data4 = pd.read_csv(file4, sep=';', index_col='ID')

	All_file = [data1, data2,data3,data4]
	# Creation of diagram de venn
	#Figure(All_file,output)
	# Creation of final combine file
	combine_annotation = Fusion_file(data1,data2,data3,data4,All_file,output)
	return(combine_annotation)

# Filter annotation for Analyse patient against gnomAD
def filter_annotation(annotation,out):
	"""
		Input: Fichier d'annotation avec les variants des 4 techniques
		Output: Fichier Filtré sur les variants d'interet sur GnomAD
		Si AF_raw > 0.01 ce variant est supprime
		warning AF topmax contient autre catégorie
	"""
	#initialization of list of line removed by filter
	remove = []
	# Filtre polymorphisme
	filter_raw = 0.01
	# Filter in GnomAD for polymorphism
	for ID in annotation.itertuples():

		if annotation.loc[ID.Index,"AF_raw"] != ".":
			if float(annotation.loc[ID.Index,"AF_raw"]) >= filter_raw:
				remove.append(ID.Index)

		else:	
			annotation.loc[ID.Index,"FILTER"]="PASS"

	# Remove line
	final_list = sorted(list(set(remove))) 
	print("Filter GnomAD:\nSur {:4d} lignes au total, il y a {:4d} lignes eliminées par les filtres donc {:4d} lignes sauvegardées.".format(annotation.shape[0],len(final_list),annotation.shape[0] - len(final_list)))
	# Suppression des lignes
	annotation.drop(final_list, inplace=True)
	return(annotation)

# Remove variant of file fusion
def remove_artefact(file_annotation,dico,out):
	"""
		Input: Fichier d'annotation filtré et fichier des artefacts
		Output: Fichier d'annotation sans les artefacts
	"""
	# Liste des Index du fichier d'annotation
	list_index = list(file_annotation.index)
	remove = []
	# Ouverture du fichier des artefact
	for artefact in dico["Artefact"].keys():
		artefact_f = artefact.strip(",")
		# S'il i exixte dans l'index suppression de la ligne
		if artefact_f.startswith("chr") == True and artefact_f in list_index:
			remove.append(artefact_f)

	final_list = sorted(list(set(remove)))
	
	print("Filter Artefact:\nSur {:4d} lignes au total, il y a  {:4d} lignes eliminées par les filtres donc {:4d} lignes sauvegardées".format(file_annotation.shape[0],len(final_list),file_annotation.shape[0] - len(final_list)))
	# Suppresseion des lignes artefacts
	file_annotation.drop(final_list, inplace=True)
	file_annotation.to_csv("Final_" + out + ".csv", sep = ';')

# Create Column Freq_database à un seuil de patient atteint at colonne avec un transcript
def final_statistic_database (name_patient,file_database,dico_annotation,out):
	# Creation Column Pourcentage de fois trouvé dans la database 
	# Ouverture du fichier de statistic de database
	"""
		IN: File of annotation,dictionnary of annotation
		OUT: Final file of annotation with column of frequence and transcript
	"""
	file_patient = pd.read_csv(name_patient, sep=';', index_col='ID')
	# Appel Fonction
	# Creation de la colonne FreqDatabase
	file_patient = freq_database(file_patient,file_database)
	# Creation de la colonne AA.chage.RegGene_transcrit
	file_patient = transcript_selection(file_patient,dico_annotation)
	# Suppression des populations Gnomad > 1%
	file_patient_popmax = mark_popmax(file_patient)
	# Warning si ERROR revoir la fonction
	#Replacer colonne
	cols = list(file_patient.columns.values)
	reorder=[]
	# All line VAF
	liste_VAF = ["Freq_Database","Mean_VAF","Mean_DP","Method_Variant_GATK","Total_DP_GATK","VAF_GATK","Method_Variant_Mutect2","Total_DP_Mutect2","VAF_Mutect2","Method_Variant_Varscan","Total_DP_Varscan","VAF_Varscan","Method_Variant_Pindel","Total_DP_Pindel","VAF_Pindel"]
	# Annotation Primaire
	liste_annot = ["CHROM","Validation_BIO_ING","Classification_mutation","cytoBand","IGV","POS","REF","ALT","Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","Favorite_transcript"]
	# Liste méthode (necessaire suppresion)
	liste_methode = ["Method_Variant_GATK","Method_Variant_Mutect2","Method_Variant_Varscan","Method_Variant_Pindel"]
	# Liste annotation secondaire
	liste_annotation = ["ID_cosmic91","OCC_cosmic91","ID_cosmic89","OCC_cosmic89","Genomic_IARC","Transactivation_Structure_Function_IARC","AF_popmax","AF_male","AF_female","AF_raw","AF_afr","AF_sas","AF_amr","AF_eas","AF_nfe","AF_fin","AF_asj","AF_oth","non_topmed_AF_popmax","non_neuro_AF_popmax","non_cancer_AF_popmax","controls_AF_popmax","CLNALLELEID_2020","CLNDN_2020","CLNDISDB_2020","CLNREVSTAT_2020","CLNSIG_2020","SIFT_pred","Polyphen2_HDIV_pred","Polyphen2_HVAR_pred","MutationTaster_pred","FATHMM_pred","M-CAP_pred","CADD_phred","Interpro_domain","ICGC_Id","ICGC_Occurrence"]
	for i in liste_annot:
		reorder.append(i)
	# Ajout des VAF seulement et pas de Methode variant
	for el in liste_VAF:
		if el not in liste_methode:
			reorder.append(el)
	#Ajout des annotations
	for annot in liste_annotation:
		reorder.append(annot)
	reorder.append("FILTER")
	reorder.append("ANNOVAR_DATE")
	#reindexation des colonnes
	file_patient_reorder = file_patient_popmax.reindex(columns =reorder)
	# Write of file
	file_patient_reorder.to_csv(out, sep = ';')

# Creation de la colonne FreqDatabase
def freq_database(file_patient,file_database):
	"""
		IN: : matrix of annotation
		OUT: matrix with information of freq_database
	"""
	file_patient["Freq_Database"] = 0
	list_index = list(file_patient.index)
	#Colonne Database
	data_base = open(file_database,"r")
	for ligne in data_base:
		line = ligne.split("\t")
		ID = line[0]	
		# Si la ligne dans dictionnary statistic est reconnu 
		# ecriture dans le fichier d'annotation		
		if ID.startswith("chr") == True and ID in list_index:
			percent = float(line[4])
			# Ecriture de la valeur
			file_patient.loc[ID,"Freq_Database"] = percent
			continue
	return (file_patient)

# Transcript
def transcript_selection(annotation_file,dico_annotation):
	"""
		INPUT: File of annotation
		OUTPUT: Creation of new column AAChange.refGene_transcript with only one transcript
	""" 

	# Creation column
	annotation_file['AAChange.refGene_transcript'] = "."
	# Parcours Fichier patient
	for ID in annotation_file.itertuples():
		# Parcours dictionnaire
		for gene in dico_annotation["Transcript"].keys():
			# Recherche du gene
			if re.search(gene,annotation_file.loc[ID.Index,'Gene.refGene']):
				# Creation d'une liste
				liste_AAchange = annotation_file.loc[ID.Index,'AAChange.refGene'].split(",")
				# Parcours liste
				for AAchange in liste_AAchange:
					# Recherche du bon transcript
					if re.search(dico_annotation["Transcript"][gene],AAchange):
						if AAchange == "ASXL1:NM_015338:exon12:c.1927dupG:p.G646Wfs*10":
							AAchange = "ASXL1:NM_015338:exon12:c.1934dupG:p.G646Wfs*12"

						elif AAchange == "ASXL1:NM_015338:exon12:c.1888_1910del:p.E635Rfs*13":
							AAchange = "ASXL1:NM_015338:exon12:c.1900_1922del:p.E635Rfs*15"

						# Notation
						annotation_file.loc[ID.Index,'Favorite_transcript'] = AAchange

	return(annotation_file)

# Note variant with gnomAD d'une population >1% tn SNP et non des catégorie
def mark_popmax(annotation_file):
	"""
		IN : File annotation
		OUT: File annnotation en notant les SNP de population
	"""
	# Liste des colonnes à chercher
	gnomad_liste = ["AF_afr","AF_sas","AF_amr","AF_eas","AF_nfe","AF_fin","AF_asj","AF_oth"]
	# Parcours du tableau
	filtre_SNP = 0.01
	for ID in annotation_file.itertuples():
		# Parcours de la liste
		for population in gnomad_liste:
			if annotation_file.loc[ID.Index,population] != ".":
				if float(annotation_file.loc[ID.Index,population]) >= filtre_SNP:
					annotation_file.loc[ID.Index,"Classification_mutation"] = "SNP"
					break
	return(annotation_file)


# Initialisation du dictionnaire
def initialize_dico():
	"""
	Initialise database and rename si un changement dans les fichiers de biologiste
	"""
	dico_annotation = {}
	dico_annotation["Artefact"] = {}
	dico_annotation["Transcript"] = {}
	return(dico_annotation)

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
		dico_annotation = json.load(infile)
	return(dico_annotation) 

def Insertion_dico_transcrit(dico_annotation,file_artefact,file_transcript):
	"""
	IN: file artefact and file transcript
	OUT: Dictionnary with information
	"""
	# Artefact
	version = file_artefact.split("_")[2].split(".")[0]
	print("version Artefect : ",version)
	# Artefact

	with open(file_artefact,"r") as liste_artefact:

		tmp = ""
		for art in liste_artefact:
			artefact = art.replace("\n","")
			if artefact.startswith("chr") == True:
			
				if artefact not in dico_annotation["Artefact"]:
					dico_annotation["Artefact"][artefact] = "OUI"

				elif artefact in dico_annotation["Artefact"]:
					print("doublon",artefact)
				else:
					print("error",artefact)
				
			tmp = artefact

	# Transcript	
	with open(file_transcript,"r") as liste_transcript:
		for line in liste_transcript:
			# On ne lit lit pas l'en tête
			if line.startswith("GENE") == False:
				gene_transcript = line.split(",")
				gene = gene_transcript[0]
				transcript = gene_transcript[1].split("\n")[0]
				if gene not in dico_annotation["Transcript"]:
					#dico_annotation["Transcript"] = gene
					dico_annotation["Transcript"][gene] = transcript
	
	return (dico_annotation)

# *****************************************************
# ********************* Main **************************
# *****************************************************
if __name__ == '__main__':

	parser = argparse.ArgumentParser(prog="Treatment_of_Annotation_v2")
	# Create sub_parsers (one for Treatment_Annotation, other for Combine)
	# Directory
	parser.add_argument("-d", "--directory",
		dest="directory",
		required=False,
		type=str,
		help="directory of result annotation for annovar"
	)
	# File
	parser.add_argument("-f", "--fileresult",
		dest="file",
		required=False,
		type=str,
		help="File of result annotation for annovar"
	)
	# Output
	parser.add_argument("-o", "--outout",
		dest="out",
		required=False,
		type=str,
		help="Name of new file in repertory of directory"
	)
	# Call method
	parser.add_argument("-m", "--method",
		dest="method",
		required=False,
		type=str,
		help="Method of calcul of variant calling or fusion : GATK Mutect2 Varscan Pindel or All"
	)
	parser.add_argument("-i", "--info_dico_annotation",
		dest="info_dico_annotation",
		required=False,
		type=str,
		help="Open of dico annotation"
	)
	# *****************************
	# Create Dictionnary argument
	parser.add_argument("-c", "--create_dico_annotation",
		dest="dico_annotation",
		required=False,
		type=str,
		help="Creation of dico annotation"
	)

	# File Artefact
	parser.add_argument("-a", "--artefact",
		dest="artefact",
		required=False,
		type=str,
		help="Last file of base artefact"
	)
	# File transcript
	parser.add_argument("-t", "--trancript",
		dest="transcript",
		required=False,
		type=str,
		help="Last file of trancript information"
	)
	# Add new ifo in same dictionnary
	parser.add_argument("-n", "--newartefact",
		dest="newartefact",
		required=False,
		type=str,
		help="Integration of new artefact information in same database"
	)
	
	# ******************
	#  End Final Statistic Database
	parser.add_argument("-s", "--databasestat",
		dest="database",
		required=False,
		type=str,
		help="Exit of statistic dictionnary for frequence information"
	)


	args = parser.parse_args()
	# Liste option
	Variant_calling = ["GATK","Mutect2","Varscan","Pindel"]
	# Launch code
	# Launch Dictionnary of information in annotation
	if args.dico_annotation:
		
		# Initialize dico
		dico = initialize_dico()
		#Load Dictionnary
		#dico = load_dictionnary(dico_annotation)
		# Implementation of dico
		print("Insertion_dico_transcrit")
		dico_annot = Insertion_dico_transcrit(dico,args.artefact,args.transcript)
		# Write dictionnary
		write_dictionnary(dico_annot,args.dico_annotation)

		sys.exit(0)
	
	# Deplacement dans repertoire
	os.chdir(args.directory) 

	# Preparation annotation file by method
	if args.method in Variant_calling:
		# Preparation du fichier d'annotation et calcul du VAF
		data_annotation = preparation_file(args.file,args.out,args.method)
		# Réalisation des filtres simple 
		filter_annotation_dico(data_annotation,args.out,args.method)
	#  Combination of file
	elif args.method == "All":

		# Combine figure and annotation
		annotation_combine = Combine_figure(args.file,args.out)
		
		# Filter
		data = filter_annotation(annotation_combine,args.out)
		# Remove variant artefact
		dico = load_dictionnary(args.info_dico_annotation)
		remove_artefact(data,dico,args.out)
	# Final statistic^et création du rapport final
	elif args.method == "Statistic":
		dico_annotation = load_dictionnary(args.info_dico_annotation)
		final_statistic_database(args.file,args.database,dico_annotation,args.out)
	# Error
	else:
		print("error or call of function")
		sys.exit(1)
	print("Le programme a traité le fichier {:30s}.".format(args.file))
	sys.exit(0)


