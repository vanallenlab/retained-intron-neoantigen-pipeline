# ----------------------------------------------------------------------------------------------- #
# Claire Margolis
# kmaToPeptideSeqs.py

# Summary: Reads in KMA output .flat_filtered.csv file, extracts chromosome locations for each 
# unique intron, then focuses on the most biologically relevant scenario: A nucleotide sequence 
# that starts in preceding exon and has at least 1 AA in intron (we already know where these ORFs 
# start and their orientation wrt the intron start site). Translates to protein and then stores the 
# sequence output in a file, as well as the list of unique introns. 

# Input format: python kmaToFasta.py ____.flat_filtered.csv 10 outdirpath
# 	"9" is default for netMHCI (9 AA window = 27 bases before intron start)
# 	"15" is default for netMHCII (15 AA window = 45 bases before intron start)

# Output format: 
# 	uniqueIntronList.txt (list of unique introns in .flat_filtered.csv file) 
#	peptideSeqsFASTA.txt (list of peptide sequences in FASTA format)

# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Import necessary packages

#!/usr/bin/python
import sys
import numpy as np
import subprocess
from Bio.Seq import Seq
import bisect
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: createUniqueIntronList()
# Inputs: kma .flat_filtered.csv file, outfile path
# Returns: List of unique intron chromosomal locations, list of corresponding TPM values
# Summary: Reads in intron chromosome locations from KMA output .csv file, extracts unique 
# sequences and their TPM values, saves them in an array for use in subsequent functions, 
# writes them to .txt file.
def createUniqueIntronList(csvfile, outpath):
	# Read in chromosome locations and TPM values (columns 1 and 4) from KMA output file
	chromlocs = np.loadtxt(csvfile, dtype=str, delimiter=',', skiprows=1, usecols=[1,3])
	# Only extract unique chromosomal locations
	_, indices = np.unique(chromlocs[:,0], return_index=True)
	uniquelocs = chromlocs[indices,:]
	uniquelocs = np.core.defchararray.strip(uniquelocs, '"')  # Strip "s from locations
	# Write contents of uniquelocs array to file (for reference later)
	filepath = outpath+'/'+'uniqueIntronList.txt'
	np.savetxt(filepath, uniquelocs, fmt='%s', delimiter='\t')
	# Return unique list of intron locations and corresponding TPM values
	return list(uniquelocs[:,0]),list(uniquelocs[:,1])
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: manualTranslate
# Inputs: FASTA sequence to translate to protein
# Returns: List of protein sequences corresponding to specific FASTA sequence
# Summary: Uses manual codon-->AA dictionary to walk through and translate FASTA sequence. 
def manualTranslate(fastasequence):
        # Initialize codon table and list of peptides
        codontable = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

        # Translate full length fasta sequence all the way through
        fulllengthprotein = ''
        for i in xrange(0, len(fastasequence), 3):
                codon = fastasequence[i:i+3]
                # Account for bizarre edge cases that should really never happen
                if len(codon) != 3 or codon not in codontable:
                        break
                if codontable[codon] == '*':  # Stop translating when we hit a stop codon
                        break
                AA = codontable[codon]
                fulllengthprotein += AA
	
	return fulllengthprotein
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Function: getSeqs() 
# Inputs: list of unique intron chromosome locations and their expression values, window size (
# number of AAs around start site)
# Returns: None
# Summary: Parses out chromosome and intron start location, then expands window size according to 
# cases (1) and (2) requirements, calls UCSC table browser twoBitToFa() function to get 
# corresponding FASTA sequences, writes to output file. Also writes corresponding file mapping 
# headers 
def getSeqs(intronlocs, tpms, nAAs, outpath):
	outfile = open(outpath + '/peptideSeqsFASTA.txt','a')
	headermapfile = open(outpath + '/headermap.txt','a')
	for l in range(0, len(intronlocs)):
		loc = intronlocs[l]
		tpmval = tpms[l]
		chrom = loc.split(':')[0]
		intronstart = loc.split(':')[1].split('-')[0]
		intronend = loc.split(':')[1].split('-')[1]
		# Query UCSC table browser to find whether sequence is on minus or plus strand and get ORF orientation
		sqlcommand = "SELECT strand,exonStarts,exonEnds,exonFrames FROM wgEncodeGencodeBasicV19 WHERE chrom='"+chrom+"' AND txStart<"+intronstart+" AND txEnd>"+intronend
		fullcommand = 'mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -A --connect_timeout=60 -e "'+sqlcommand+'"'
		tablebrowserout = subprocess.check_output(fullcommand, shell=True)
		# Catch instance where table browser doesn't think there's actually a gene in this region
		if len(tablebrowserout) < 1:
			continue
		# Find the transcript with the most exons/introns if there is more than one and use that 
		strand = ''
		exonstarts = []
		exonends = []
		exonframes = []
		tablebrowserlist = filter(None,tablebrowserout.split('\n'))
		infostring = ''
		lengthholder = 0
		for i in range(1,len(tablebrowserlist)):
			if (len(filter(None,tablebrowserlist[i].split(','))) > lengthholder):
				infostring = tablebrowserlist[i]
				lengthholder = len(filter(None,infostring.split(',')))
		strand = infostring.split('\t')[0]
		exonstarts = filter(None,infostring.split('\t')[1].split(','))
		exonstarts = map(int, exonstarts)
		exonends = filter(None,infostring.split('\t')[2].split(','))
		exonends = map(int, exonends)
		exonframes = filter(None,infostring.split('\t')[3].split(','))
		exonframes = map(int, exonframes)
		# Check to make sure exon starts, ends, and frames list are the same length, and if there is an error, skip this intron
		if not len(exonstarts) == len(exonends) == len(exonframes):
			continue
		# Get ORF orientation at the start of the intron
		# Handle + and - strand cases separately (need to look for different positions)
		ORForientation = 0
		frame = 0
		if strand == ('+'):
			index = bisect.bisect_left(exonstarts, int(intronstart)) - 1
			frame = exonframes[index]
			if frame == 0:
				ORForientation = (int(intronstart)-exonstarts[index]) % 3
			elif frame == 1:
				ORForientation = (int(intronstart)-2-exonstarts[index]) % 3
			elif frame == 2:
				ORForientation = (int(intronstart)-1-exonstarts[index]) % 3
			else: # If frame == -1 (meaning no translation takes place according to table browser)
				continue
		else: 
			index = bisect.bisect_left(exonstarts, int(intronend))
			frame = exonframes[index]
			if frame == 0:
				ORForientation = (exonends[index]-int(intronend)) % 3
			elif frame == 1: 
				ORForientation = (exonends[index]-2-int(intronend)) % 3
			elif frame == 2: 
				ORForientation = (exonends[index]-1-int(intronend)) % 3 
			else: 
				continue
		# Determine nucleotide window around which to get sequence
		wholeseqstart = 0
		wholeseqend = 0
		if strand == ('+'):
			wholeseqstart = int(intronstart) - ORForientation - (nAAs*3-3)
			wholeseqend = int(intronend) + (nAAs*3-3)
		else:
			wholeseqstart = int(intronstart) - (nAAs*3-3)
			wholeseqend = int(intronend) + ORForientation + (nAAs*3-3)
		# Get genomic sequence
		loc = chrom + ':' + str(wholeseqstart) + '-' + str(wholeseqend)
		loc = '-seq='+loc
		twobitoutput = (subprocess.check_output(['/xchip/cga_home/margolis/Packages/tableBrowser/twoBitToFa',
                        loc,'/xchip/cga_home/margolis/General/hg19.2bit','stdout']))
                # Parse output
		seqlist = twobitoutput.split('\n')
		headerline = seqlist[0]+"|"+tpmval
		seqlist = seqlist[1:len(seqlist)-1]
		sequence = ''.join(str(elem) for elem in seqlist)
		sequence = sequence.upper()
		# Reverse complement if it's on the negative strand
		if strand == ('-'):
			sequence = str(Seq(sequence).reverse_complement())
		# Manually translate sequence
		peptide = manualTranslate(sequence)
		# Check to make sure peptide is at least "length" AAs long, and if so, write to output file
		if len(peptide) < nAAs:
			continue
		else:
			newheaderline = '>seq'+str(l)
                	outfile.write(newheaderline+'\n')
                	outfile.write(peptide+'\n')
			headermapfile.write(newheaderline+'\t'+headerline+'\n')
	
	return

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Main function that processes command line input and calls other functions
def main():
	# Check to make sure we have the right number of inputs
	if len(sys.argv) != 4:
		print 'Error: incorrect number of inputs.'
		print 'Please input a KMA .csv file, the AA window you want, and an outfile path.'
		sys.exit()
	# Store inputs
	kmafile = sys.argv[1]
	window = int(sys.argv[2])
	outpath = sys.argv[3]
	# Create unique intron output file
	uniqueIntrons,TPMvals = createUniqueIntronList(kmafile, outpath)
	# Create nucleotide sequences file
	getSeqs(uniqueIntrons, TPMvals, window, outpath)

if __name__ == '__main__':
    main()
# ----------------------------------------------------------------------------------------------- #

