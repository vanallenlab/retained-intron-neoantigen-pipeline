# ----------------------------------------------------------------------------------------------- #
# Claire Margolis
# 1 March 2017
# addTPMValsToOutput.py
#
# Summary: Function to add TPM values for each retained intron peptide to the final output 
# 	   for each patient so that we can look at expression of RIs that yield neoantigens.
#
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Import necessary packages

#!/usr/bin/python
import sys
import numpy as np
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: addTPM
# Inputs: netMHC file, kma file, and outfile path
# Returns: None (writes output to a file)
# Summary: Takes in processed netMHC file and individual patient KMA file, finds the TPM for each
# neoantigen by doing a bit of manipulation and cross-referencing the kma file, and returns 
# an augmented matrix containing netMHC cleaned output + TPM values for each patient.
def addTPM(netmhc, kma, path): 
	# Read in netMHC results file
	netmhcarray = np.loadtxt(netmhc, dtype=str, delimiter='\t', skiprows=1)
	# Read in KMA file and create dictionary matching RI location to TPM value
	with open(kma) as f:
		kmalines = f.read().splitlines()
	del kmalines[0]  # Remove header line
	tpmdict = {}
	for line in kmalines:
		key = line.split(',')[1]
		value = float(line.split(',')[4])
		# Manipulate key a bit to get it in the same format as the netMHC output
		chrom = key.split(':')[0]
		coord = key.split(':')[1].split('-')[0]
		coord = int(coord)-30
		newkey = chrom[1:len(chrom)] + '_' + str(coord)
		tpmdict[newkey] = value
	# Loop through netmhcarray and retrieve TPM values for each neoantigen by referencing dict
	rowstodelete = []
	tpmlist = []
	for i in range(0, netmhcarray.shape[0]):
		row = netmhcarray[i,:]
		key = row[2].split('-')[0]
		tpm = tpmdict[key]
		tpmlist.append(tpm)
	augmentedarray = np.delete(netmhcarray, rowstodelete, axis=0)
	# Add expression column to end of array
	tpmcol = np.asarray(tpmlist).reshape(len(tpmlist),1)
	augmentedarray = np.hstack((augmentedarray, tpmcol))
	# Write thresholded array to output file
	headerstring = 'Pos\tPeptide\tID\tcore\tallele\t1-log50k\tnM\trank\texp_TPM'
	np.savetxt(path+'/final_neoantigen_catalog.txt', augmentedarray, fmt='%s', delimiter='\t', header=headerstring, comments='')

	return
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Main function
def main():
        # Check to make sure we have the right number of inputs
        if len(sys.argv) != 4:
                print 'Error: incorrect number of inputs.'
                print 'Please input a processed netMHC output file, an individual patient kma results file, and a valid outfile path'
                sys.exit()
        # Read in inputs
        netmhcfile = sys.argv[1]
        kmafile = sys.argv[2]
        outfilepath = sys.argv[3]
        # Process netMHCpan file 
        addTPM(netmhcfile, kmafile, outfilepath)

        return

if __name__ == '__main__':
        main()
# ----------------------------------------------------------------------------------------------- #

