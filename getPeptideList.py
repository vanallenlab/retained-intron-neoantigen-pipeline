# ----------------------------------------------------------------------------------------------- #
# Claire Margolis
# getPeptideList.py
#
# Summary: 
#       Creates dataframe of patient, response status, peptide, chromosome location, HLA type for 
#	each cohort (will eventually be combined for multiple cohorts into a master list.
#
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Import necessary packages

#!/usr/bin/python
import sys
import numpy as np
import subprocess
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: makePepDF
# Inputs: File containing patient names and response statuses, path to patient directories, out
# Returns: None (writes to output file)
# Summary: Creates a dataframe with every RI binder peptide for each patient in cohort plus
# additional information.
def makePepDF(patsandstatus, pathtopats, outpath):
	# Read in array and isolate patient names and response statuses
        array = np.loadtxt(patsandstatus, dtype=str, delimiter='\t', ndmin=2)
        pats = array[:,0]
        status = array[:,1]
	# Initialize relevant lists (which will be combined eventually into an array)
	IDs = []
	statuses = []
	peptides = []
	chromlocs = []
	hlas = []
	nMs = []
	ranks = []
	for i in range(0, len(pats)):
                # Read in each patient's results file and add information to master lists
                results = np.loadtxt(pathtopats+'/'+pats[i]+'/'+pats[i]+'processedNETMHCpan_out.txt', dtype=str, delimiter='\t', skiprows=1, usecols=(1,2,4,6,7))
		for row in results:
			IDs.append(pats[i])
			statuses.append(status[i])
			peptides.append(row[0])
			chromlocs.append(row[1])
			hlas.append(row[2])
			nMs.append(row[3])
			ranks.append(row[4])
	# Once all patients' information is contained in lists, combine into np array
	overall_array = np.column_stack(tuple([IDs, statuses, peptides, chromlocs, hlas, nMs, ranks]))
	# Write to output file
	headerstring = 'PatientID\tResponseStatus\tPeptideBinder\tGenomicLocation\tHLA\tnM\trank'
	np.savetxt(outpath+'/aggregateRIBindersArray.txt', overall_array, fmt='%s', delimiter='\t', comments='', header=headerstring)	

	return
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: main
def main():
        # Check to make sure we have the right number of inputs
	if len(sys.argv) != 4:
		print 'Error: incorrect number of inputs.'
		print 'Please input a .txt file containing patient directory names (path not included) + response status, a path to patient directories, and an out path.'
		sys.exit()
        # Read in inputs
	filelist = sys.argv[1]
	pathstring = sys.argv[2]
	outfilepath = sys.argv[3]
        # Gather and aggregate patient summary information 
	makePepDF(filelist, pathstring, outfilepath)

	return

if __name__ == '__main__':
	main()
# ----------------------------------------------------------------------------------------------- #

