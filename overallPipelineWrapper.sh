# ------------------------------------------------------------------------------------------------ #
# Claire Margolis
# Last updated: 31 March 2017
# overallPipelineWrapper.sh

# Summary: Shell script to run entire retained intron-neoantigen pipeline from start to finish. 

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
use Java-1.8
use Bowtie2
use Samtools

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Set relevant paths

patient_dir=pat_dirs.txt
PAT_DIR=$(cat $patient_dir | head -n $SGE_TASK_ID | tail -n 1)
fastq1_fls=fastq1_fls.txt
FASTQ1=$(cat $fastq1_fls | head -n $SGE_TASK_ID | tail -n 1)
fastq2_fls=fastq2_fls.txt
FASTQ2=$(cat $fastq2_fls | head -n $SGE_TASK_ID | tail -n 1)
bam_fls=bam_fls.txt
BAM_FILE=$(cat $bam_fls | head -n $SGE_TASK_ID | tail -n 1)
kma_bam_fls=kma_bam_fls.txt
KMA_BAM_FILE=$(cat $kma_bam_fls | head -n $SGE_TASK_ID | tail -n 1)

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run bamToFastq on each patient
# ( Converts patient RNA-Seq bam files back into fastq files )

echo 'Running bamToFastq.'
nice java -jar /xchip/cga_home/dmiao/software/picard-tools-1.140/picard.jar SamToFastq I=$BAM_FILE FASTQ=$FASTQ1 SECOND_END_FASTQ=$FASTQ2

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run Bowtie on each patient
# ( Realigns patient FASTQ files to an augmented transcriptome including intronic regions )

echo 'Running Bowtie.'
bowtie2 -k 200 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 -x /xchip/cga_home/asmart/KMA/160705_output/trans_and_introns -1 /xchip/Hugo_Mel_PD1/Snyder_ipi/bamtofastq/$FASTQ1 -2 /xchip/Hugo_Mel_PD1/Snyder_ipi/bamtofastq/$FASTQ2 | samtools view -Sb - > $KMA_BAM_FILE

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run eXpress on each patient
# ( Quantifies expression by counting reads )

echo 'Running eXpress.'
/xchip/cga_home/asmart/KMA/express-1.5.1-linux_x86_64/express /xchip/cga_home/asmart/KMA/160705_output/gencode.v19.pc_transcripts_with_introns.fa $KMA_BAM_FILE -o ../$PAT_DIR

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run zero coverage filter on each patient
# ( Pre-processing step necessary for KMA )

echo 'Running zero coverage filter.'
python /xchip/cga_home/asmart/KMA/kma/pre-process/zeroCoverage.py ../$PAT_DIR/results.xprs $KMA_BAM_FILE ../$PAT_DIR/zero_coverage.txt

cp ../$PAT_DIR/zero_coverage.txt ../$PAT_DIR/zero_coverage_2.txt

sed 's/\..*|ENSG.*|//g' ../$PAT_DIR/zero_coverage_2.txt > ../$PAT_DIR/zero_coverage_3.txt

cp ../$PAT_DIR/results.xprs ../$PAT_DIR/results_2.xprs

sed 's/\..*|ENSG.*|//g' ../$PAT_DIR/results_2.xprs > ../$PAT_DIR/results_3.xprs
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run KMA on entire cohort
# ( Identifies retained introns in each patient for whole cohort )

echo 'Running KMA.'
Rscript run_kma.R

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
python /xchip/cga_home/margolis/retainedIntron/goldStandard/postprocessOutput.py $PAT_DIR/NETMHCpan_out.xls $PAT_DIR $PAT_DIR

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Run aggregateSampleInfo.py once to aggregate all patient results
# ( Takes previous step output for all patients and aggregates into one document for whole cohort )

#echo 'Running aggregateSampleInfo.py.'
#python /xchip/cga_home/margolis/retainedIntron/goldStandard/aggregateSampleInfo.py patientnamesandresults.txt /xchip/cga_home/____/patientDirs CohortName /xchip/cga_home/____/outfilepath

# ----------------------------------------------------------------------------------------------- #
