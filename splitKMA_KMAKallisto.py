# ----------------------------------------------------------------------------------------------- #
# Claire Margolis
# splitKMA_KMAKallisto.py

# Summary: Takes in KMA file, patient ID, and outfile path, and outputs a modified 
# KMA file containing only the retained intron locations that correpond to the patient of 
# interest. This is meant to be run in a shell script that will batch run this for each patient 
# in the cohort. Specifically designed for kma-kallisto version.  

# Input format: python splitKMA.py kmaOutfile.csv patientID patientIDoutfile.csv

# Output format: patientIDoutfile.csv

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Import necessary packages

#!/usr/bin/python
import sys
import numpy as np
import subprocess

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: splitFile
# Inputs: original kma file, patient ID, outfile
# Returns: none (writes to file)
# Summary: Takes in original KMA file and includes only rows with patient ID matching our 
# patient of choice and intron counts > 0 AND TPM > 1.
def splitFile(kmafile, patient, outfile):
	# Read in KMA file
	with open(kmafile) as f:
	    lines = f.read().splitlines()
	# Open output file for writing
	out = open(outfile, 'w')
	# Write header to file
	out.write(lines[0]+'\n')
	# Loop through to get only lines that belong to our patient, write those to outfile
	for i in range(1, len(lines)):
		line = lines[i]
		currpatient = line.split(',')[2]
		currTPM = line.split(',')[3]
		currcounts = line.split(',')[7]
		currTPMfilter = line.split(',')[8]
		currPSIfilter = line.split(',')[9]
		currcountsfilter = line.split(',')[10]
		if currpatient == '"'+patient+'"' and float(currcounts) > 0 and float(currTPM) > 1:
			if str(currTPMfilter) in "TRUE" and str(currPSIfilter) in "TRUE" and str(currcountsfilter) in "TRUE":
				out.write(line+'\n')
	# Close outfile
	out.close()

	return

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Main function
def main():
	# Check to make sure we have the right number of inputs
	if len(sys.argv) != 4:
		print 'Error: incorrect number of inputs.'
		print 'Please input a KMA .csv file, valid patient ID, and outfile path'
		sys.exit()
	# Read in inputs
	kmafile = sys.argv[1]
	patient = sys.argv[2]
	outfile = sys.argv[3]
	# Split kma file
	splitFile(kmafile, patient, outfile)

	return

if __name__ == '__main__':
	main()

# ----------------------------------------------------------------------------------------------- #

