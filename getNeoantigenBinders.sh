# ------------------------------------------------------------------------------------------------ #
# Claire Margolis
# 23 February 2017
# getNeoantigenBinders.sh

# Summary: Shell script to run a series of python scripts to go through the pipeline from 
# KMA to netMHC for all patients in sample
# *NOTE*: If you want to run this script, go through and verify that the paths to relevant files
# are in the correct format for your cohort. You will need to change out_dir.txt, among other 
# things, to make the script specific to your cohort. You can also change preferences for running
# netMHCIpan vs. netMHCIIpan. 

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Specify shell / UGER preferences

#!/bin/bash

#$ -cwd
#$ -q long
#$ -m e
#$ -l h_vmem=30g 
#$ -t 1-39

# ----------------------------------------------------------------------------------------------- #
# Use statements

source /broad/software/scripts/useuse
reuse Python-2.7
use MySQL-5.6
use EMBOSS

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Set directory paths

patient_dir=/xchip/cga_home/margolis/retainedIntron/VA_Mel_ipi/out_dir.txt
PAT_DIR=$(cat $patient_dir | head -n $SGE_TASK_ID | tail -n 1)

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run splitKMA.py on each patient
# ( Splits aggregate KMA output into patient-specific KMA output )

echo 'Running splitKMA.py.'
python /xchip/cga_home/margolis/retainedIntron/goldStandard/splitKMA.py /xchip/cga_home/asmart/KMA/VA_Mel_IPI/160822_VA_Mel_IPI_KMA_IR_flat_filtered_ex41_v1.csv $PAT_DIR $PAT_DIR/kma_results.txt

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run kmaToPeptideSeqs.py on each patient
# ( Converts kma output into FASTA peptide files for each retained intron )

echo 'Running kmaToPeptideSeqs.py.'
python /xchip/cga_home/margolis/retainedIntron/goldStandard/kmaToPeptideSeqs.py $PAT_DIR/kma_results.txt 9 $PAT_DIR

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run runNetMHCpan.py on each patient
# ( Runs netMHCpan with retained-intron peptides and HLA-alleles specific to each patient )
# NOTE:: Use "1" to run netMHCpanI program and "2" to run netMHCpanII

echo 'Running runNetMHCpan.py.'
python /xchip/cga_home/margolis/retainedIntron/goldStandard/runNetMHCpan.py $PAT_DIR/peptideSeqsFASTA.txt ../$PAT_DIR/hla_alleles.txt '1' $PAT_DIR

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run postprocessOutput.py on each patient
# ( Processes netMHCpan output file to a more user-friendly, relevant format )

echo 'Running postprocessOutput.py.'
python /xchip/cga_home/margolis/retainedIntron/goldStandard/postprocessOutput.py $PAT_DIR/NETMHCpan_out.xls $PAT_DIR/headermap.txt $PAT_DIR $PAT_DIR

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run aggregateSampleInfo.py once to aggregate all patient results
# ( Takes previous step output for all patients and aggregates into one document for whole cohort )

#echo 'Running aggregateSampleInfo.py.'
#python /xchip/cga_home/margolis/retainedIntron/goldStandard/aggregateSampleInfo.py patientnamesandresults.txt /xchip/cga_home/____/patientDirs CohortName /xchip/cga_home/____/outfilepath

# ----------------------------------------------------------------------------------------------- #
