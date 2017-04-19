# ----------------------------------------------------------------------------------------------- #
# Claire Margolis
# 11 October 2016
# postprocessOutput.py
#
# Summary: Takes in NETMHC_out.xls file and does postprocessing to create a more user-friendly 
# output format.
# Input format: python postprocessOutput.py NETMHCpan_out.xls patientID outpath
# Output format: processedNETMHCpan_out.txt
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Import necessary packages
#!/usr/bin/python
import sys
import numpy as np
import subprocess
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Function: processFile
# Inputs: netMHCpan output file, patient ID, outpath
# Returns: None (writes to file) 
# Summary: Postprocesses the netMHCpan output to eliminate useless rows, change data format from 
# wide to long, add allele name columns. Writes new data to a file in outpath. 
def processFile(filename, patID, outpath):
	# Read in first line of file to get number and names of alleles
	with open(filename, 'r') as f:
	    alleles = f.readline().strip().split('\t')
	alleles = filter(None, alleles)
	data = np.loadtxt(filename, dtype=str, delimiter='\t', skiprows=2)
	nrow = data.shape[0]
	ncol = data.shape[1]
	# Remove all rows with last column (NB) == 0
	for i in range(0, nrow):
		data[i,-1] = data[i,-1].strip()
	data = data[data[:,-1] != '0']
	# Move columns so data is in long form
	nrow = data.shape[0]
	listofarrays = []  # Will store all allele-specific arrays
	initcols = data[:,0:3]  # Initial three columns that are common to all HLA alleles
	for i in range(0, len(alleles)):
		currstartcol = (3*(i+1))+i
		currendcol = currstartcol+4
		currarray = data[:,currstartcol:currendcol]
		listofarrays.append(currarray)
	datalong = np.vstack(tuple(listofarrays))
	# Add initial columns and allele column into data frame
	# Allele column
	allelevec = []
	for i in range(0, len(alleles)):
		currnewcol = [alleles[i]]*nrow
		allelevec.extend(currnewcol)
	datalong = np.insert(datalong, 1, allelevec, axis=1)  # Add allele column to datalong
	# Initial columns
	initcollist = []
	for i in range(0, len(listofarrays)):
		initcollist.append(initcols)
	initcolstoappend = np.vstack(tuple(initcollist))
	updateddata = np.concatenate((initcolstoappend, datalong), axis=1)
	# Eliminate any columns that have a rank above 2
	toremove = []
	updatednrows = updateddata.shape[0]
	for i in range(0, updatednrows):
		if float(updateddata[i,7]) > 2:
			toremove.append(i)
	updateddata = np.delete(updateddata, toremove, 0)
	# Write updated data to new file
	outfilepath = outpath+'/'+patID+'processedNETMHCpan_out.txt'
	np.savetxt(outfilepath, updateddata, fmt='%s', delimiter='\t', header='Pos\tPeptide\tID\tcore\tallele\t1-log50k\tnM\tRank') 

	return
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Main function
def main():
        # Check to make sure we have the right number of inputs
        if len(sys.argv) != 4:
                print 'Error: incorrect number of inputs.'
                print 'Please input a netMHCpan output file, patient ID, and valid outfile path'
                sys.exit()
        # Read in inputs
        netmhcfile = sys.argv[1]
	patientID = sys.argv[2]
        outfilepath = sys.argv[3]
        # Process netMHCpan file 
        processFile(netmhcfile, patientID, outfilepath)

        return

if __name__ == '__main__':
        main()
# ----------------------------------------------------------------------------------------------- #

