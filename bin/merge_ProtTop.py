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

	fi = open("prot.top", "r")
	fo = open("prot_lig.top", "w")
	flag = 1
	sol_count = 0
        for line in fi:
                if line.find("SOL") > -1:
                        sol_count += 1
	fi.seek(0)
	sol_flag = 1
	flag_chain = 0
	for line in fi:
		if line.find("; Include chain topologies") > -1:
			flag_chain += 1
			fo.write(line)
			fo.write("#include \"lig_GMX.itp\"\n")
		elif line.find("; Compound        #mols") > -1:
			fo.write(line)
			fo.write("SUB         1\n")
		elif line.find("SOL") > -1:
			if sol_flag == 1:
				if sol_count == 1:
					fo.write(line)
				elif sol_count > 1:
					fo.write("SOL          %d\n"%sol_count)
			sol_flag = 0
		else:
			fo.write(line)
	fi.close()
	fo.close()		

	fi = open("prot_lig.top", "r")
	filelines = fi.readlines()
	fi.close()
	print "\nchain flag"
	print flag_chain
        if flag_chain < 1:
                fi = open("prot_lig.top", "w")
                for i in range(len(filelines)):
			if filelines[i].find("; Include forcefield parameters") > -1:
				numLine = i
				break
		for i in range(numLine):
			fi.write(filelines[i])
		for i in range(numLine, numLine+2):
			fi.write(filelines[i])
                fi.write('\n\n; Include chain topologies\n')
                fi.write('#include \"lig_GMX.itp\"\n\n')
		
  		for i in range(numLine+2, len(filelines) ):
                        fi.write(filelines[i])
                fi.close()
        else:
                pass
	
	os.system("mv prot.top prot_bak.top")
	os.system("cp prot_lig.top prot.top")
main(sys.argv[1:])
