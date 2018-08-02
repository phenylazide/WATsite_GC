#! /usr/bin/python

#############################################################
# prot_lig.py
# run MD simulation for computing solvation effects 
#
#############################################################

import os, sys, getopt

ligand_code = "SUB"
#############################################################
def main(argv):
	inlig = argv[0]
	print inlig
	bak = "lig_GMX.itp_bak"
	os.system("mv %s %s" % (inlig, bak ) )
	fi = open(bak, "r")
	fo = open(inlig, "w")
	flag = 0
	for line in fi:
		if line.find("[ moleculetype ]") > -1:
			flag = 1
		if flag == 1:
			if line.find("lig") > -1:
				line_1 = "SUB   " + line[6:]
			else:
				line_1 = line
			fo.write(line_1)
		else:
			fo.write(line)
        fi.close()
        fo.close()

main(sys.argv[1:])
