# ----------------------------------------------------------------------------------------------- #
# Claire Margolis
# aggregateSampleInfo.py
#
# Summary: 
#       Script to aggregate summary information (number of RI locations, number of binders, etc.)
#	across all samples in a cohort. 
#
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Import necessary packages

#!/usr/bin/python
import sys
import numpy as np
import subprocess
from collections import Counter
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: aggregateData
# Inputs: Patient directory/name and response status list, path to patient directory location, cohort name, out path
# Returns: None (writes to outfile)
# Summary: For every patient in list, retrieves number of total RIs, number of total binders, 
# number of strong binders, number of weak binders, and HLA type. Writes to a table in output file.
def aggregateData(patsandstatus, path, cohort, out):
	# Read in array and isolate patient names and response statuses
	array = np.loadtxt(patsandstatus, dtype=str, delimiter='\t', ndmin=2)
	pats = array[:,0]
	status = array[:,1]
	# Initialize lists to contain each of the relevant pieces of information
	numRIs = []
	numRIsyieldingbinders = []
	numtotalbinders = []
	numstrongbinders = []
	numweakbinders = []
	hlatypes = []
	# Loop through each patient and gather summary statistics
	for i in range(0, len(pats)):
		# Get number of RIs and add to list
		command = 'wc -l < '+path+'/'+pats[i]+'/uniqueIntronList.txt'
		numlines = subprocess.check_output(command, shell=True).strip()
		numRIs.append(numlines)
		# Get number of binders and number RIs yielding binder(s)
		netmhcarray = np.loadtxt(path+'/'+pats[i]+'/'+pats[i]+'processedNETMHCpan_out.txt', dtype=str, delimiter='\t', skiprows=1)
		tempstrong = []
		numRIsyieldingbinders.append(len(set(netmhcarray[:,2])))
		for i in range(0, len(netmhcarray[:,6])):
                	if float(netmhcarray[i,6]) < 50.0: # May want to change this to if float(netmhcarray[i,7]) < 0.5: to use rank versus affinity
                        	tempstrong.append(netmhcarray[i,1])
        	numstrongbinders.append(len(list(set(tempstrong))))
		numtotalbinders.append(len(list(set(netmhcarray[:,1]))))
		numweakbinders.append(len(list(set(netmhcarray[:,1])))-len(list(set(tempstrong))))
		alleles = sorted(list(set(netmhcarray[:,4])))
		hlatypes.append(','.join(alleles))
	# Calculate mean stats for overall, and each response category, and write directly to file
	outputfile = open(out+'/RI_cohort_summary.txt', 'a')
	outputfile.write('# --------------------------- Retained Intron Summary Document for '+cohort+' --------------------------- #\n\n')
	outputfile.write('Overall stats:: Mean RI load: '+str(np.mean(np.array(numRIs).astype(np.float)))+'; Mean num RIs yielding binder(s): '+str(np.mean(np.array(numRIsyieldingbinders).astype(np.float)))+'; Mean num total binders: '+str(np.mean(np.array(numtotalbinders).astype(np.float)))+'; Mean num strong binders: '+str(np.mean(np.array(numstrongbinders).astype(np.float)))+'; Mean num weak binders: '+str(np.mean(np.array(numweakbinders).astype(np.float)))+'\n\n')	
	uniqueresponses = list(set(status))
	responsecounts = []
	for response in uniqueresponses:
		indices = [i for i, x in enumerate(status) if x == response]
		responsecounts.append(str(response)+': '+str(len(indices)))
        	meanRI = meanpepRI = meanTB = meanSB = meanWB = 0
		for i in range(0, len(indices)):
			meanRI += float(numRIs[indices[i]])
			meanpepRI += float(numRIsyieldingbinders[indices[i]])
			meanTB += float(numtotalbinders[indices[i]])
			meanSB += float(numstrongbinders[indices[i]])
			meanWB += float(numweakbinders[indices[i]])
		meanRI = meanRI/len(indices)
		meanpepRI = meanpepRI/len(indices)
		meanTB = meanTB/len(indices)
		meanSB = meanSB/len(indices)
		meanWB = meanWB/len(indices)
		outputfile.write(response+' statistics:: Mean RI load: '+str(meanRI)+'; Mean num RIs yielding binder(s): '+str(meanpepRI)+'; Mean num total binders: '+str(meanTB)+'; Mean num strong binders: '+str(meanSB)+'; Mean num weak binders: '+str(meanWB)+'\n\n')
	# Write counts of response types to file
	outputfile.write('Response distribution: '+'; '.join(responsecounts)+'\n\n')
	# Calculate distribution of HLA types and write to file
	hlalist = []
	for allelestring in hlatypes:
		indivalleles = allelestring.split(',')
		for ind in indivalleles: 
			hlalist.append(ind)
	alleletable = Counter(hlalist)
	outputfile.write('HLA type distribution amongst patients: '+str(alleletable)+'\n\n')
	# Convert lists into single numpy array and append to outfile
	summaryarray = np.column_stack(tuple([pats, status, numRIs, numRIsyieldingbinders, numtotalbinders, numstrongbinders, numweakbinders, hlatypes]))
	headerstring = 'Individual stats:\nPatient\tResponseStatus\tTotalRIload\tNumberRIsYieldingBinder(s)\tNumberTotalPeptideBinders\tNumberTotalStrongBinders\tNumberTotalWeakBinders\tHLAType'
	np.savetxt(outputfile, summaryarray, fmt='%s', delimiter='\t', comments='', header=headerstring)

	return
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Main function
def main():
        # Check to make sure we have the right number of inputs
        if len(sys.argv) != 5:
                print 'Error: incorrect number of inputs.'
                print 'Please input a .txt file containing patient directory names (path not included) + response status, a path to patient directories, the cohort name, and an out path.'
                sys.exit()
        # Read in inputs
        filelist = sys.argv[1]
	pathstring = sys.argv[2]
	cohortname = sys.argv[3]
        outfilepath = sys.argv[4]
        # Gather and aggregate patient summary information 
        aggregateData(filelist, pathstring, cohortname, outfilepath)

        return

if __name__ == '__main__':
        main()
# ----------------------------------------------------------------------------------------------- #

