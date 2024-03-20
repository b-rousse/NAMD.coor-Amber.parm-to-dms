# NAMD.coor-Amber.parm-to-dms
This repo assists in converting a system simulated in NAMD with Amber parameters into an Anton/Desmond- compatible .dms file. The process involves two steps: 1. Converting the NAMD output into an Amber simulation output, and 2. combines the new Amber output file with the Amber parameters used in the NAMD simulation to form a .cms format file, which can be used by the mae2dms function in the Anton supercomputer, as detailed on the Pittsburgh Supercomputing Center (PSC) website.

Step 1: Convert the NAMD output into an Amber simulation output
Usage:
vmd -dispdev text -python -e convertNAMDtoRST-final-system.py -args -p amber-parmfile.top -c system-from-namd.restart.coor -v system-from-namd.restart.vel -x system-from-namd.xsc -o outfile -s -S "protein"

You must modify convertNAMDtoRST-final-system.py: on line 250, enter the number of atoms in your system

Note: this conversion was only tested on version 1.9 of VMD.

Step 2: Combine this new amber-type .rst file with the amber FF files to generate the .cms file
Usage: 

module load desmond/31023-bin

$SCHRODINGER/run -FROM mmshare ~/amber_prm2cms_v.py -p amber-parmfile.top -c output.rst -o system.cms

Note: details on this latter step are provided on the Anton/Anton2 website via the PSC.
