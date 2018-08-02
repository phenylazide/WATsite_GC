import os, sys, re, getopt
WATsite_home = os.environ['WATSITEHOME']
def getCharge():
	# This is a check to see if the reduce_prep.pdb file is there
	try:
		Ri = open("charge.log", "r")
		try:
			pass
		finally:
			Ri.close()
	except IOError:
			print("\nError: can not find the file charge.log\n")
			sys.exit(0)
	
	# Open the reduce_prep.pdb file, obtain HIS protonation states, sort by res#, write out logfile
	Ri = open("charge.log", "r")
	
	for line in Ri:
		if line.find("System has non-zero total charge:") > -1:
			charge = float(line.split(':')[1])
			int_charge = int(round(charge))
			break
	print int_charge
getCharge()
