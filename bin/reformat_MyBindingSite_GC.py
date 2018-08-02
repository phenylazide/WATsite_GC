import os, sys, re, getopt
#WATsite_home = os.environ['WATSITEHOME']
MyBindingSite_pdb = os.environ['MyBindingSite_pdb']
def reformat_ligand():
	#shift ligand
	if not os.path.exists(MyBindingSite_pdb): raise IOError('can not find the MyBindingSite file')
	Ri = open(MyBindingSite_pdb, "r")
	newfile = open("MyBindingSite_formatted.pdb", "w")
	for line in Ri:
		if ((line.find("ATOM") > -1) or (line.find("HETATM") > -1)):
			newline = line[:17] + "BDS" + line[20:]
			newfile.write(newline)
		elif line.find("TER") > -1: continue
		else:
			newfile.write(line)
	newfile.close()
	Ri.close()
reformat_ligand()
