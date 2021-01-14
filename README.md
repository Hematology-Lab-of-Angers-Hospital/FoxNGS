# Pipeline SMPHD_hemato
# *********************
#  Pipeline bioinformatique SMPHD                 
#	Données NGS                   
#   Version 3.2, version essai
# *********************

# Tutoriel
Code principale pipeline_SMPHD.sh
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
- ajout dbsnp
- Cosmic 89 et Cosmic 92s
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
