#############################################################
# pdb2mol2.py by Bingjie Hu 04/25/2015
#############################################################

import os, os.path, sys, getopt
def main(argv):
    try:
	opts, args = getopt.getopt(argv, "hi:o:") 
    except getopt.GetoptError:          
	sys.exit(2)                     
    for opt, arg in opts:
	if opt == '-h':
	    print usage                     
	    sys.exit()                  
	elif opt == '-i':
	    mypdb = arg
	elif opt == '-o':
            mymol2 = arg

    try:
	Ri = open(mypdb, "r")
	try:
            pass
        finally:
            Ri.close()
    except IOError:
	print("\nError: can not find the file tmp_shift.pdb\n")
	sys.exit(0)
    Ri = open(mypdb, "r")
    Ro = open(mymol2, "w")

    num_atoms = 0
    for line in Ri:
        if line.find("ATOM") > -1:
            num_atoms += 1
    Ri.seek(0)

    Ro.write("@<TRIPOS>MOLECULE\n")
    Ro.write("****\n")
    Ro.write("%5d %5d %5d %5d %5d\n"%(num_atoms, 0, 0, 0, 0))
    Ro.write("SMALL\n")
    Ro.write("USER_CHARGES\n\n\n")
    Ro.write("@<TRIPOS>ATOM\n")
    for line in Ri:
        if line.find("ATOM") > -1:
            x = float(line[31:38])
            y = float(line[39:46])
            z = float(line[47:54])
            occupancy = float(line[55:60])  #deltaS
            B_fact = float(line[61:66])     #deltaG
            newline = "      1 WAT  %8.3f %8.3f %8.3f O         1 <0>        %6.2f  %6.2f\n"%(x, y, z, B_fact, occupancy)
            Ro.write(newline)
    Ri.close()
    Ro.close()
    
main(sys.argv[1:])
