# Pipeline SMPHD_hemato
# *********************
#  Pipeline bioinformatique SMPHD                 
#	Données NGS                   
#   Version 3.1, version stable
# *********************

# Tutoriel
Code principale pipeline_SMPHD_V3.0.sh
appel des codes python

# Versioning tools
## Tools
- samtools 1.9
- bwa mem  
## Call variant
- Varscan 2.4.3
- GATK Haplotype caller et Mutect 4.1.5
- Pindel
## Annotation
- Cosmic 89 et Cosmic 91
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
- Ajout de Cosmic 92
- Remplacement de la database clinvar 20200316 à la place de clinvar20190305
- Affichage du % of target dans R
# 3.1
- NPM1
