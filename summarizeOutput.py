# ----------------------------------------------------------------------------------------------- #
# Claire Margolis
# 21 October 2016
# summarizeOutput.py
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Import necessary packages

#!/usr/bin/python
import sys
import numpy as np
import subprocess
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: createFile
def createFile(fname, patID, response, out, validationset):
	array = np.loadtxt(fname, dtype=str, delimiter='\t', skiprows=1)
	with open(validationset) as f:
		validset = f.read().splitlines()
	strongbinders = []
	for i in range(0, len(array[:,6])):
		if float(array[i,6]) < 50.0:
			strongbinders.append(array[i,1])
	numstrongbinders = len(list(set(strongbinders)))
	allpeps = array[:,1]
	alleles = list(set((array[:,4])))
	numtotalbinders = len(list(set(allpeps)))
	numweakbinders = numtotalbinders - numstrongbinders
	pepsinvalidationset = list(set(allpeps) & set(validset))
	out = open(patID+'_SummaryStats.txt','w')
	out.write('Patient ID\tResponseStatus\tNumStrongBinders\tNumWeakBinders\tNumTotalBinders\tPepsInValidationSet\tAlleleString\n')
	allelestring = ','.join(alleles)
	pepstring = []
	if not pepsinvalidationset:
		pepstring = 'FALSE'
	else:
		pepstring = ','.join(pepsinvalidationset)
	out.write(patID+'\t'+response+'\t'+str(numstrongbinders)+'\t'+str(numweakbinders)+'\t'+str(numtotalbinders)+'\t'+pepstring+'\t'+allelestring+'\n')
	out.close()
	

	return


# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def main():
	filename = sys.argv[1]
	patientID = sys.argv[2]
	responderstatus = sys.argv[3]
	outpath = sys.argv[4]
	validatedpepfile = sys.argv[5]
	createFile(filename, patientID, responderstatus, outpath, validatedpepfile)

	return

if __name__ == '__main__':
    main()
# ----------------------------------------------------------------------------------------------- #
