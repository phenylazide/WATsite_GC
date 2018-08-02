import os, sys, re, getopt
WATsite_home = os.environ['WATSITEHOME']
myligand = os.environ['myligand']
def reformat_ligand():
	#shift ligand
	try:
		Ri = open(myligand, "r")
		try:
			pass
		finally:
			Ri.close()
	except IOError:
			print("\nError: can not find the ligand file\n")
			sys.exit(0)
	
	Ri = open(myligand, "r")
	newfile = open("MyBindingSite_formatted.pdb", "w")
	for line in Ri:
		if ((line.find("ATOM") > -1) or (line.find("HETATM") > -1)):
			newline = line[:17] + "BDS" + line[20:]
			newfile.write(newline)
		elif line.find("TER") > -1:
                        continue
                else:
			newfile.write(line)
	newfile.close()
	Ri.close()
reformat_ligand()
