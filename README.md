# Pipeline SMPHD_hemato
# *********************
#  Pipeline bioinformatique SMPHD                 
#	Données NGS                   

# *********************

# Tutoriel
Code principale pipeline_SMPHD.sh
appel des codes python

# Versioning tools
## Tools
- samtools 1.9
- samtools 1.10 en soutien coverage non rendu
- bwa mem  
## Call variant
- Varscan 2.4.3
- GATK Haplotype caller et Mutect 4.1.5
- Pindel
## Annotation
- ajout avsnp138
- Cosmic 89 et Cosmic 92
- cytoBand
- gnomad211_exome
- clinvar_20200316
- dbnsfp35a (score)
- IARC,icgc21

# Amelioration
# **********
## 3.0
- Amélioration du menu de lancement du pipeline
- Inclusion de cosmic 91 à la place de cosmic 90
- Remplacement de la database clinvar 20200316 à la place de clinvar20190305
# 3.1
- Cosmic 92
- NPM1
- Bed target %

- Ajout de Cosmic 92
- Remplacement de la database clinvar 20200316 à la place de clinvar20190305
- Affichage du % of target dans R
# 3.1
- NPM1
## 3.3
- Analyse of target
- Inclusion de cosmic 92 à la place de cosmic 91
- Ajout de dbsnp138
- Affichage du % of target dans R
- Ajout dans le fichier csv

