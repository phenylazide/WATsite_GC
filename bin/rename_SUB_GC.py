#! /usr/bin/python

#############################################################
# prot_lig.py
# run MD simulation for computing solvation effects 
#
#############################################################

import os
import os.path
import sys
import getopt
import math
import re
import time
import glob
import shutil
import string
import fnmatch
import operator

ligand_code = "SUB"
#############################################################
def main(argv):
	inlig = argv[0]
	print (inlig)
	fi = open(inlig, "r")
	fo = open("lig_tmp.mol2", "w")
	line = fi.readline()
	while line != '':
		if line.find("@<TRIPOS>MOLECULE") > -1:
			fo.write(line)
			line = fi.readline()
			fo.write(line)
			line = fi.readline()
			fo.write(line)
			numatoms = int(line.split()[0])
		elif line.find("@<TRIPOS>ATOM") > -1:
			fo.write(line)
			hcount = 1
			for i in range(numatoms):
				line = fi.readline()
				ss = line.split()
				if ss[1][0] == 'H':
					atom_name = 'H%d'%hcount
					hcount += 1
				else:
					atom_name = ss[1]
				newline = "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 SUB       % 8.4f\n"%(int(ss[0]), atom_name, float(ss[2]), float(ss[3]), float(ss[4]), ss[5], float(ss[-1]))
				fo.write(newline)
		else:
			fo.write(line)
		line = fi.readline()
	fi.close()
	fo.close()		
	os.system("mv lig_tmp.mol2 ligand.mol2")
	#os.system("mv %s lig_bak.mol2"%inlig)
	#os.system("mv tmp.mol2 %s"%inlig)
main(sys.argv[1:])
