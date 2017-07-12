# ----------------------------------------------------------------------------------------------- #
# Claire Margolis
# runNetMHCpan.py

# Summary: Takes in a fasta file containing all peptide sequences upon which netMHCpan (either I, 
# II, or both) is to be run, runs netMHCpan I or II (or both) depending on the user's selection, 
# writes output to a tab-delimited text file.

# Input format: python runNetMHCpan.py FASTAproteinsequences.txt HLAalleles.txt '1,2' outpath
# 	Options for specifying which netMHCpan version: 
# 	'1,2' = both netMHCIpan and netMHCIIpan
#	'1' = netMHCIpan
#	'2' = netMHCIIpan
# *RELEVANT*: HLA allele input file can be in one of two formats: 
#	1. Polysolver winners_hla.txt output file
# 		example line from file: HLA-A   hla_a_02_01_01_01       hla_a_32_01_01
# 	2. Already processed, one allele per line in netMHC compatible format
#		example line from file: HLA-A02:01

# Output: netMHCpan output .xls file(s)

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Import necessary packages

#!/usr/bin/python
import sys
import numpy as np
import subprocess

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: runNetMHCIpan
# Inputs: FASTA file of peptide sequences, patient HLA alleles (in a format specified above), outpath 
# Returns: None (netMHCIpan will automatically write output to a .xls file
# Summary: Pre-processes patient HLA alleles, runs netMHCIpan. 
def runNetMHCIpan(pepfile, hlafile, outpath):
	# Read in HLA alleles file and process
	with open(hlafile) as f:
		hlalines  = f.read().splitlines()
	hlaalleles = []
	# Determine which input format the hla allele file is in
	if len(hlalines[0].split('\t')) <= 1:  # In already pre-processed format
		hlaalleles = hlalines
	else:  # Polysolver output file
		for line in hlalines:
			split = line.split('\t')
			# Reformat each allele (2 for each type of HLA A, B, and C)
			for i in range(1, 3):
				currallele = 'HLA-'
				allele = split[i]
				components = allele.split('_')
				currallele += components[1].upper() + components[2] + ':' + components[3]
				hlaalleles.append(currallele)
	hlaalleles = list(set(hlaalleles))  # Remove duplicate alleles if there are any
	hlastring = ','.join(hlaalleles)
	# Run netMHCI pan
	command =  'export NHOME=/xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0; export NETMHCpan=/xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/Linux_x86_64; /xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/Linux_x86_64/bin/netMHCpan -a '+hlastring+' -f '+pepfile+' -inptype 0 -l 9,10 -s -xls -xlsfile '+outpath+'/NETMHCpan_out.xls -allname /xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/Linux_x86_64/data/allelenames -hlapseudo /xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/Linux_x86_64/data/MHC_pseudo.dat -t 500 -version /xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/data/version -tdir /xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/scratch/XXXXXX -rdir /xchip/cga_home/margolis/Packages/netMHCPan/netMHCpan-3.0/Linux_x86_64/ > '+outpath+'/netMHCpanout.txt'
	subprocess.call(command, shell=True)
	
	return

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Main function
def main():
	 # Check to make sure we have the right number of inputs
	if len(sys.argv) != 5:
		print 'Error: incorrect number of inputs.'
		print 'Please input a FASTA file, a HLAalleles.txt file, a netMHCpan version, and an outpath.'
		sys.exit()
	# Parse inputs
	fasta = sys.argv[1]
	alleles = sys.argv[2]
	versionchoice = sys.argv[3]
	outpath = sys.argv[4]
	# Split versionchoice
	versions = versionchoice.split(',')
	if len(versions) == 1:
		if versions[0] == '1':
			runNetMHCIpan(fasta, alleles, outpath)
		else:
			runNetMHCIIpan(fasta, alleles, outpath)
	else:
		runNetMHCIpan(fasta, alleles, outpath)
		runNetMHCIIpan(fasta, alleles, outpath)
			
	return

if __name__ == '__main__':
    main()

# ----------------------------------------------------------------------------------------------- #

