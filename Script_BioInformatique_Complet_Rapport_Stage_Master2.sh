#!/bin/bash

set -euo pipefail
trap 'echo "? Erreur ligne $LINENO : commande \"$BASH_COMMAND\""; exit 1' ERR

#SBATCH --output=/work/user/efoucaud/My_Data_Rance/Rapport_De_Script/merge_%j.out
#SBATCH --error=/work/user/efoucaud/My_Data_Rance/Rapport_De_Script/merge_%j.err
#SBATCH -p unlimitq
#SBATCH -t 00-20:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=800G

########################################################################################################################################################################
########################################################################################################################################################################
## Script 1 - Appariement et fusion des reads paired-end (Merging R1 & R2) - Tools : VSEARCH v2.22.1

module load bioinfo/VSEARCH/2.22.1

DATA_DIR="/work/user/efoucaud/My_Data_Rance"
cd "$DATA_DIR" || exit 1
mkdir -p ./output/merged_data

echo "Début du merging..."

merged_count=0
failed_count=0

for R1_file in *_R1.fastq; do
    R2_file=${R1_file/_R1.fastq/_R2.fastq}
    if [[ -f "$R2_file" ]]; then
        full_base_name=$(basename "$R1_file" _R1.fastq)
        vsearch --fastq_mergepairs "$R1_file" \
                --reverse "$R2_file" \
                --fastqout "./output/merged_data/${full_base_name}_merged.fastq" \
                --fastq_minovlen 10 \
                --fastq_maxdiffs 10 \
                --fastq_allowmergestagger \
                --fastq_qmax 50 \
                --threads 8 \
                --quiet
        if [ $? -eq 0 ]; then
            merged_count=$((merged_count + 1))
        else
            failed_count=$((failed_count + 1))
            echo "Problème merge : $full_base_name"
        fi
    else
        echo "R2 manquant pour $R1_file"
        failed_count=$((failed_count + 1))
    fi
done

echo "Merging terminé"
echo "Réussis : $merged_count"
echo "Échecs : $failed_count"

########################################################################################################################################################################
########################################################################################################################################################################
## Script 2 - Definition of molecular marker variables - Exemple : COImg2 ! But identical for MiFish (12S) - Tools : OBITools 4.4.15

DATA_DIR="/work/user/efoucaud/My_Data_Rance/output"
WORK_DIR="$DATA_DIR/COImg2"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

module load devel/Miniconda/Miniconda3
module load bioinfo/OBITools/4.4.15

#Définition des variables
BarcodeName="CO1mg2"  ####Metabarcode Name
 Fwd1="TCHACHAAYCAYAARGAYATYGG" ####Fwd primer sequence
 Rev1="ACYATRAARAARATYATDAYRAADGCRTG" ####Rev primer sequence
 RevC1="CAYGCHTTYRTHATRATYTTYTTYATRGT" ####Rev primer sequence reverse complemented
MinSeqLenght="120" #### Minimum sequence length 
MaxSeqLenght="200" #### Minimum sequence length 

########################################################################################################################################################################
########################################################################################################################################################################
## SCRIPT 3 - Removal of primers and selection of sequences of interest + removal of poor-quality reads - Tools : Cutadapt v4.3

module load bioinfo/Cutadapt/4.3

echo "Début du Cutadapt" 

cutadapt -g $Fwd1 -O 10 -e 0.1 --discard-untrimmed /work/user/efoucaud/My_Data_Rance/output/all_samples_annoted.fastq -o /work/user/efoucaud/My_Data_Rance/output/COImg2/cutadapt_intermediaire.fastq --revcomp &&
cutadapt -a $RevC1 -O 10 -e 0.1 --discard-untrimmed cutadapt_intermediaire.fastq -o /work/user/efoucaud/My_Data_Rance/output/COImg2/RANCE_${BarcodeName}_trim_ali.fastq --revcomp --match-read-wildcards 

echo "Fin du Cutadapt"

########################################################################################################################################################################
########################################################################################################################################################################
## SCRIPT 4 - Removal of empty sequences - Tools SeqKit v2.12.0

module load bioinfo/SeqKit/2.12.0 

echo "Début Seqkit"

seqkit seq -m 1 RANCE_CO1mg2_trim_ali.fastq > RANCE_CO1mg2_trim_ali_noEmpty.fastq

echo "Fin Seqkit"

########################################################################################################################################################################
########################################################################################################################################################################
## SCRIPT 5 - Read dereplication - Tools : Obiuniq

echo "Début Obiuniq"

obiuniq --input-OBI-header -m sample RANCE_${BarcodeName}_trim_ali_noEmpty.fastq > RANCE_${BarcodeName}_trim_ali_uniq.fasta

echo "Fin Obiuniq"

########################################################################################################################################################################
########################################################################################################################################################################
## SCRIPT 6 - Dereplicated annotation sequences - Obiannotate

echo "Début Obiannotate"

obiannotate -k count -k merged_sample /work/user/efoucaud/My_Data_Rance/output/COImg2/RANCE_${BarcodeName}_trim_ali_uniq.fasta > RANCE_${BarcodeName}_trim_ali_uniq_simple.fasta

echo "Fin Obiannotate"

########################################################################################################################################################################
########################################################################################################################################################################
## SCRIPT 7 - Obigrep 1 - Sequences Length Selection - COImg2 : 100 pb – 200 pb | MiFish : 138 pb – 203 pb - Tools : Obigrep

echo "Début Obigrep"

obigrep -l 100 -L 200 RANCE_${BarcodeName}_trim_ali_uniq_simple.fasta > RANCE_${BarcodeName}_trim_ali_uniq_simple_grep.fasta

echo "Fin Obigrep"

########################################################################################################################################################################
########################################################################################################################################################################
## SCRIPT 8 - Obigrep 2 - Filtering by minimum abundance (threshold: 500 occurrences) - Tools : Obigrep 

echo "Début Obigrep"

obigrep -c 500 RANCE_${BarcodeName}_trim_ali_uniq_simple_grep.fasta > RANCE_${BarcodeName}_trim_ali_uniq_simple_grep_500.fasta

echo "Fin Obigrep"

########################################################################################################################################################################
########################################################################################################################################################################
## SCRIPT 9 - Denoising and chimera detection — Tools : obiclean 

echo "Début Obiclean"

obiclean -s sample -r 0.05 --detect-chimera -H RANCE_${BarcodeName}_trim_ali_uniq_simple_grep_500.fasta > RANCE_${BarcodeName}_trim_ali_uniq_simple_grep_500_clean.fasta

echo "Fin Obiclean"

########################################################################################################################################################################
########################################################################################################################################################################
## SCRIPT 10 : Final Sequences Annotation - Tools : obiannotate

echo "Début Obiannotate"

obiannotate --number RANCE_${BarcodeName}_trim_ali_uniq_simple_grep_500_clean.fasta | obiannotate --set-id 'sprintf("seq%04d",annotations.seq_number)' > RANCE_${BarcodeName}_trim_ali_uniq_simple_grep_500_clean_final.fasta

echo "Début Obiannotate"

########################################################################################################################################################################
########################################################################################################################################################################
## SCRIPT 11 : Abundance Table Exportation (CSV format)

echo "Début ObiMatrix et ObiCSV"

obimatrix --map obiclean_weight RANCE_${BarcodeName}_trim_ali_uniq_simple_grep_500_clean_final.fasta > RANCE_${BarcodeName}_trim_ali_uniq_simple_grep_500_clean_final.csv

echo "Fin ObiMatrix et ObiCSV"

########################################################################################################################################################################
########################################################################################################################################################################
## Pipeline complete — Abundance table ready 
## Sequences file ready for taxonomic assignment using Kraken / WoRMS / GBIF and analyses using R (metabaR / phyloseq)





















