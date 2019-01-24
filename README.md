# retained-intron-neoantigen-pipeline

This pipeline calls RNA-based neoantigens from intron retention events derived from RNA-Seq data and identified through the KMA package (see run instructions below for further detail on this).

To run: 
- Download NetMHCPan-3.0 (http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan) and change paths in runNetMHCpan.py file (line 62).
- Download twoBitToFa utility from UCSC genome browser (https://genome.ucsc.edu/goldenpath/help/twoBit.html) and change paths in kmaToPeptideSeqs.py file (line 173).
- Download MySQL (you will use it to query the UCSC table browser via public servers).
- Download and run KMA package (https://github.com/pachterlab/kma). The output from the KMA package (ir$flat file) will be the direct input to this pipeline.
- Change paths in shell script getNeoantigenBinders.sh (notes in file comments).
- Run getNeoantigenBinders.sh from command line as an SGE Array Job. This script is a wrapper and will call all other relevant Python scripts.

Additional notes:
- Detailed execution instructions and functionality descriptions can be found in each script header, as well as for each individual function.
- Feel free to create an Issue if errors arise.
