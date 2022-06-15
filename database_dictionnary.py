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
def dictionnary_init():
    """
	Initialise database and rename if
	if a significant change in the pipeline is made
	"""
    dict = {}
    dict["patient_sample"] = []

    return(dict)

# Ecriture du dictionnaire
def dictionnary_write(write_dico, name):
    """
    Writes the dictionnary in JSON formmat
    """ 
    with open(name, 'w') as outfile:
        json.dump(write_dico, outfile)

# Chargement du dictionnaire
def dictionnary_load(name):
    """
    Gets the dictionnary file name and loads it as a JSON
    """

    with open(name,"r") as infile:
        dict = json.load(infile)

    return(dict)

# Insertion des valeurs dans le dictionnaire
def dictionnary_insert(dict,name_file):
    """
        IN : Database
        OUT: Dictionnary with insertion of data
    """

    # Preparation des elements pour bien ordonner le dictionnaire
    # Partie "Patient"
    VAF = ["Mean_VAF","Method_Variant_GATK","Total_DP_GATK","VAF_GATK",
	"Method_Variant_Mutect2","Total_DP_Mutect2","VAF_Mutect2",
	"Method_Variant_Varscan","Total_DP_Varscan","VAF_Varscan",
	"Method_Variant_Pindel","Total_DP_Pindel","VAF_Pindel"]
    sample = name_file.split(".")[0].split("_")[3]

    # Declaration du nom de l'sample
    # Non ecriture de TEM HORIZON dans le dictionnaire == Passage de l'insertion
    liste_exclure = ["TEM-HO","EEQ"]


    if sample.startswith(liste_exclure[0]) or sample.startswith(liste_exclure[1]):
        print('No match sample of control',sample)
        sys.exit(0)
    else:
        print("Nom de l sample en cours d insertion du dictionnaire : " , sample)


    patients = ['0409-GJ-09112010','0422-RC-19072006',"0013-SJ-30092010","0025-CR-13012010","0034-GA-29102009",
                    "0083-WP-07012011","0161-WF-07102011","0165-ND-30062011","0224-BC-10062009","0737-KJ-16112011",
                    "0113-GM-17092013","0117-GP-25012010","0294-CR-14112008","0041-LB-21032011","0628-BJ-07062016",
                    "0859-DG-16042015"]
    if sample in patients:
        print("Patient à repasser")
        sys.exit(0)
    # Case of a patient not in dictionnary
    if sample not in dict["patient_sample"]:
        dict["patient_sample"].append(sample)
    # Reading the input file
    with open(name_file,"r+") as annotation:

        for line in annotation:
			# If the line does not match the header
            if line[0:2] != "ID":
                # Getting variant ID
                info=line.split(";")
                ID=info[0]
				# Case of a new variant that is not a horizon control
                if ID not in dict:
                    # Creating a sub-dictionnary
                    dict[ID] = {}
                    # Creating a value-key for each patient with the new ID key
                    dict[ID][sample] = {}
                    # Getting each sample with a matching ID
                    dict[ID]["patient_IDs"] = []
                    dict[ID]["patient_IDs"].append(sample)
                    # Getting new information 
                    dict[ID]["Annotation"] = {}
                    # Ecriture dans le dictionnaire
                    i = 1
                    while i < len(header):    
                        attribute = header[i]
                        value = info[i]
                        # Element ID à la base de ce sous dictionnaire
                        if attribute == "ID":
                            continue
                        # Insertion des valeurs de VAF
                        elif attribute in VAF:
                            dict[ID][sample][attribute] = value
                        # Insertion des valeurs d'Annotation
                        else:    
                            dict[ID]["Annotation"][attribute] = value
                        i+=1
                
                # Implementation de valeur d'un patient pour un ID reconnu
                elif ID in dict and sample not in dict[ID]["patient_IDs"] :
                    dict[ID][sample]= {}
                    dict[ID]["patient_IDs"].append(sample)
                    i = 1
                    # Ecriture seulement des VAF pour chaque patient
                    while i < len(header):
                        attribute = header[i]
                        value = info[i]
                        # Insertion des valeurs de VAF pour le patient
                        if attribute in VAF:
                            dict[ID][sample][attribute] = value
                        i+=1
                
                # ID déja connu pour ce patient: Signalement
                elif ID in dict and sample in dict[ID]["patient_IDs"] :
                    print("Ce patient {:20s} et cet ID {:20s} a déja été stocké.".format(sample,ID))

                # Other error to understand
                else :
                    print("Other error")
                    sys.exit(1)
            # Recuperation du header du tableau
            else:
                header = line.split(";")

    return(dict)

# Statistics on the dictionnary
def statistic(dict,out):
    # outptu file: CSV of statistic
    # Setting tabs
    # patients to remove
    """
        IN: Dictionnary
        OUT: Output file with main informations
    """
    patients = ['0409-GJ-09112010','0422-RC-19072006',"0013-SJ-30092010","0025-CR-13012010","0034-GA-29102009",
                    "0083-WP-07012011","0161-WF-07102011","0165-ND-30062011","0224-BC-10062009","0737-KJ-16112011",
                    "0113-GM-17092013","0117-GP-25012010","0294-CR-14112008","0041-LB-21032011","0628-BJ-07062016",
                    "0859-DG-16042015"]
    tab = "\t" 
    with open(out,"w") as stat_out:
        
        print(" Patient",dict["patient_sample"])

        count_total = len(dict["patient_sample"])
        print("Number Patient", count_total) 
        # Mutation per gene
        stat_out.write("ID\tGene\tIGV\tcytoBand\tFreq\tPatient\tRedundancy\tExonicFunc\tNotation\n")

        # Temporary variable to spot double entries
        tmp = 0
        # Deleting partient; do not take into account
        dict["patient_sample"] = [i for i in dict["patient_sample"] if i not in patients]

        # Search in dictionnary
        for key in sorted(dict.keys()):
            
            # Double entry delete
            if tmp == key:
                print("ERROR : duplicated ID")
                sys.exit(1)
            # If it is a variant    
            if key.startswith("chr") == True:
                tmp = key
				# Get variant information
                #if key == "chr17_74732935_CGGCGGCTGTGGTGTGAGTCCGGGG_C":
                #    print(dict[key])
                
                # If the key is valid
                if len(dict[key]["patient_IDs"]) > 0:
                    # Deleting patient ID 
                    dict[key]["patient_IDs"] = [i for i in dict[key]["patient_IDs"] if i not in patients]
                    for sample in patients: 
                        if sample in dict[key]:
                            # Deleting samples
                            del dict[key][sample]    

                    # Operation in dictionnary
                    countID = len(dict[key]["patient_IDs"])
                    freq_ID = round(float(len(dict[key]["patient_IDs"]) / count_total),4)
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
                # Getting patients with this ID
                    line_patient = ""
                    for patient in dict[key]["patient_IDs"]:
                        if line_patient != "":
                            line_patient = line_patient + "," + patient
                        # First line
                        else:
                            line_patient = line_patient + patient
                    line = key + tab + str(dict[key]["Annotation"]["Gene.refGene"]) + tab + str(dict[key]["Annotation"]["IGV"]) + tab \
                    + str(dict[key]["Annotation"]["cytoBand"]) + tab  + str(freq_ID) + tab + line_patient + tab \
                    + Note_freq + tab + str(dict[key]["Annotation"]["ExonicFunc.refGene"]) + tab + str(dict[key]["Annotation"]["AAChange.refGene"]) + "\n" 
                    stat_out.write(line)
                    
                    
    return(dict)

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
        help="Name of new dictionnary and its location"
    )
    parser.add_argument("-c", "--newdictionnary",
        dest="dictionnary",
        required=False,
        type=bool,
        help="""Argument to indicate creation of new dictionnary 
                -c TRUE to trigger new dictionary creation 
                False by default
            """
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
	# If the args.dictionnary parameter is TRUE, the dictionnary is initialised
    if args.dictionnary:
        # If some notable change is made on the pipeline, run this command
        dict = dictionnary_init()
    # If the dictionnary needs no update, it is loaded
    else:
        # Load dictionnary
        dict = dictionnary_load(args.out)
        
    # if args.stat == True Statistic on dictionnary
    if args.stat:
    # Calling statistics funtion
        dico = statistic(dict,args.outstat)
		# Writing changes following a database element deletion
        dictionnary_write(dico,args.out)
    # Inserting values in the dictionnary  
    else:
        # Deplacement dans le repertoire
        os.chdir(args.directory)
        dico = dictionnary_insert(dict,args.file)
        # Save new dictionnary
        dictionnary_write(dico,args.out)
        print("Data insertion ({:20s}) in the database went well.".format(args.file))
    
    
    sys.exit(0)
