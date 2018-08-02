#! /usr/bin/python

#
# hydrodock.py
#

import os, re
import os.path
import sys
import getopt
import math
import glob

md_output_dir = ".."
#########################################################################################
#this script seperates the wtrindex_????.ndx files under the "wtrindex_cutoff_0.00" folder
#into seperate subfolders and generates the RunScripts for each subfolder. The RunScripts
#call the "prep_genergy.py" to generate the *xvg files for hydroenergy analysis. By doing
#this, we can paralelly run the whole process on the super-clusters.
#########################################################################################
def main(argv):  

	mpi = 0

	try:                                
		opts, args = getopt.getopt(argv, "hi:q:m") 
	except getopt.GetoptError:          
		sys.exit(2)                     
	for opt, arg in opts:
		if opt == '-h':
			print usage                     
			sys.exit()                  
		elif opt == '-i':
			folder_in = arg                        
		elif opt == '-m':
			mpi = 1                       
		elif opt == '-q':
			qsub = arg                       

	print "MPI = %d " % mpi
			           
	#os.chdir(folder_in)
	name_lst = glob.glob('%s/*.ndx'%folder_in)
	name_lst.sort()
	
	for ndx_file in name_lst:
		file_index = ndx_file[-8:-4]
		in_dir = "Enthalpy_output_%s" % file_index
		os.mkdir(in_dir)
		os.system("cp %s %s" % (ndx_file, in_dir))
		out_dir = "Enthalpy_output_%s" % file_index
		runscript_file = "Enthalpy_RunScript_%s" % file_index
		Energy_calc(in_dir, out_dir, runscript_file, mpi )

	os.system("chmod 777 Enthalpy_RunScript_????")
	name_lst = glob.glob('Enthalpy_RunScript_????')
	name_lst.sort()
	super_runscript = open("SuperRunScript", "w")
	i = 0
	for file in name_lst:
		if mpi == 1:
			comm = "qsub -q %s %s | tee myjob.pbsid%d\n" % (qsub, file, i)
		elif mpi == 0:
			comm = "./%s\n" % (file)
		super_runscript.write(comm)
		i += 1
	super_runscript.close()

def Energy_calc(folder_in, folder_out, runscript_file, mpi): 
	
	root = os.getcwd()
	print root
	md_output_dir = root
	E_out = folder_out
	if os.path.isdir(E_out) == 0:
		os.mkdir(E_out)
	index_folder = os.path.join(root,folder_in)
	print index_folder

	#os.chdir(index_folder)
	# count = number of index files
	count = -1
	file_list = []
	num_wtr = []
	name_lst = []
	# os.path.join (this is a system indepent way to get the directory)
	E_out = os.path.join(root,E_out)
	print E_out
	# glob.glob (way to list the file in the directory that match the .ndx ending)

	name_lst = glob.glob('%s/*.ndx'%index_folder)
	name_lst.sort()
	for fname in name_lst:
		file_list.append(fname)
		water_list = []
		path_in = os.path.join(index_folder,fname)
		file_in = open(path_in,'r')

		# reading the .ndx file and searching for Water, store to water_list, bump counter       
		for line in file_in:
			#if line.find("Water") > -1:
			if re.search('Water_[0-9]+', line):
				water_list.append(line.split()[1])
		count += 1
		num_wtr.append(len(water_list)) 
	
		# write out the pr.mdp file
		# make sure the mdp file is consistent with the original run!
		fout = open('%s/pr_%d.mdp' %(index_folder,count),'w')
		fout.write(';\n')
		fout.write('; position restraint\n')
		fout.write(';\n')
		fout.write(';\n')
		fout.write(';\n')
		fout.write('title		    =  position_restraint\n')
		fout.write('cpp                 =  /usr/bin/cpp\n')
		#becareful with this define of position restraint! make sure we used this in the 1st run of md
		fout.write('define              =  -DPOSRES\n')
		fout.write('constraints         =  hbonds\n')
		fout.write('lincs_order         =  8\n')
		fout.write('lincs_iter          =  5\n')
		fout.write('lincs_warnangle     =  40\n')
		fout.write('integrator          =  md\n')
		fout.write('dt		    =  0.002\n')
		fout.write('nsteps              =  500000\n')
		fout.write(';\n')
		fout.write('; Center-of-mass specifications\n')
		fout.write(';\n')
		fout.write('comm_mode   	    = Linear\n')
		#fout.write('comm_grps           = protein non-protein\n')
		fout.write('nstcomm		    = 1\n')
		fout.write(';\n')
		fout.write('; Output specifications\n')
		fout.write(';\n')
                fout.write('cutoff-scheme      = group\n')
		fout.write('nstxout             = 10\n')
		fout.write('nstvout             = 10\n')
		fout.write('nstfout             = 0\n')
		fout.write('nstlog              = 10\n')
		fout.write('nstenergy           = 5000\n')
		fout.write('energygrps          = ')
		# Writing all the Water names to the energygrps line
		for i in range(len(water_list)):
			fout.write('%s '%water_list[i])
		fout.write('Environment\n')
		fout.write('energygrp_excl     = Environment Environment\n')
		fout.write(';\n')
		fout.write('; Force field specifications\n')
		fout.write(';\n')
		fout.write('nstlist             = 10\n')
		fout.write('ns_type             = grid\n')
		fout.write('rlist               = 1.0\n')
		fout.write('coulombtype         = PME\n')
		fout.write('rcoulomb            = 1.0\n')
		fout.write('rvdw                = 1.4\n')
		fout.write('vdwtype             = Cut-off\n')
		fout.write('fourierspacing      = 0.12\n')
		fout.write('fourier_nx          = 0\n')
		fout.write('fourier_ny          = 0\n')
		fout.write('fourier_nz          = 0\n')
		fout.write('pme_order           = 4\n')
		fout.write('ewald_rtol          = 1e-5\n')
		fout.write('optimize_fft        = yes\n')
		fout.write(';\n')
		fout.write('; Berendsen temperature coupling is on in four groups\n')
		# make sure the Tcoupl is the same method as the original run!
		fout.write('Tcoupl 		    = nose-hoover\n')
		fout.write('tau_t 		    = 0.1 	  0.1\n')
		fout.write('tc_grps 	    = protein non-protein\n')
		fout.write('ref_t 		    = 300 	  300\n')
		fout.write(';\n')
		fout.write('; Pressure coupling is on\n')
		fout.write('Pcoupl 		    = Parrinello-Rahman\n')
		fout.write('pcoupltype 	    = isotropic\n')
		fout.write('tau_p 		    = 0.5\n')
		fout.write('compressibility     = 4.5e-5\n')
		fout.write('ref_p 		    = 1.0\n')
		fout.write(';\n')
		fout.write('; Generate velocites is on at 300 K.\n')
		fout.write('gen_vel             = yes\n')
		fout.write('gen_temp            = 300.0\n')
		fout.write('gen_seed            = 123456\n')
		fout.close()
	
        for ndx2_file in name_lst:
		fsubname = ndx2_file.split("/")
                file2_index = fsubname[-4] + fsubname[-3][:2] + "-" + ndx2_file[-6:-4]
 
    	if mpi == 1:
		print "Running GROMACS in parallel\n"
		# Run Gromacs commands for each .mdp file
		rfile = open(runscript_file, "w") 
		rfile.write("#!/bin/tcsh\n")
		rfile.write("#PBS -l nodes=1:ppn=16,walltime=03:59:00\n")
		rfile.write("#PBS -N %s\n" % file2_index)
		rfile.write("#\n")
		rfile.write("\n")
		rfile.write("cd $PBS_O_WORKDIR\n")
		rfile.write("module load gromacs\n")
		rfile.write("\n")
		rfile.write("cd %s\n"%E_out)
		for i in range(0,count+1):		
			pr_file = os.path.join(index_folder,'pr_%d.mdp' %i)
			ndx_file = os.path.join(index_folder,'%s'% file_list[i])
			comm1 = "$GROMACS_HOME/bin/grompp -f %s -c %s/prot_md.gro -p %s/prot.top -o test.tpr -maxwarn 50 -n %s -quiet\n" % (pr_file,md_output_dir,md_output_dir, ndx_file)
			print comm1
			rfile.write(comm1)

			#make sure the *xtc file is the same one as used to generate the mdsnaps.pdb file
			comm1 = "$GROMACS_HOME/bin/mdrun -ntmpi 16 -rerun %s/prot_md.xtc -s test.tpr -o test.trr -c test.gro -g test.log -e test.edr -quiet\n" % (md_output_dir)
			print comm1
			rfile.write(comm1)
			rfile.write("rm test.log test.trr \#*\n")
			
			#print 'number of water',num_wtr[i]
			rfile.write("chmod a+x g_energy*.txt\n")
			for j in range(0,num_wtr[i]):
				wtr_index = int(water_list[j].split('_')[1])
				fin = open("%s/g_energy_%d.txt"%(E_out,j), "w")
				fin.write("g_energy -f test.edr -s test.tpr -sum -o Energy_output_%d -quiet<< EOF\n" % (wtr_index)) 
				for k in range(0,j):
					fin.write("Coul-SR:%s-%s\nLJ-SR:%s-%s\nLJ-LR:%s-%s\nCoul-14:%s-%s\nLJ-14:%s-%s\n"%(water_list[k],water_list[j], water_list[k],water_list[j], water_list[k],water_list[j], water_list[k],water_list[j], water_list[k],water_list[j]))
				#for k in range(j,num_wtr[i]):
				fin.write("Coul-SR:%s\nLJ-SR:%s\nLJ-LR:%s\nCoul-14:%s\nLJ-14:%s\n"%(water_list[j],water_list[j],water_list[j],water_list[j],water_list[j]))

				fin.write("\n")

				fin.close()
				# Need to change the permissons of this file so that it can be ran % run it
				rfile.write("./g_energy_%d.txt\n"%j)

		rfile.write("rm test* \#*\n")
		rfile.close()


	elif mpi == 0:
                print "Running GROMACS rerun one by one\n"
		# Run Gromacs commands for each .mdp file
		rfile = open(runscript_file, "w") 
		#os.chdir(E_out)	
		rfile.write("#!/bin/tcsh\n")
		#rfile.write("#PBS -l procs=8,walltime=04:00:00\n")
		#rfile.write("#\n")
		#rfile.write("#\n")
		#rfile.write("\n")
		#rfile.write("cd $PBS_O_WORKDIR\n")
		#rfile.write("module load gromacs\n")
		#rfile.write("\n")
		rfile.write("cd %s\n"%E_out)
		for i in range(0,count+1):		
			pr_file = os.path.join(index_folder,'pr_%d.mdp' %i)
			ndx_file = os.path.join(index_folder,'%s'% file_list[i])
			comm1 = "$GROMACS_HOME/bin/grompp -f %s -c %s/prot_md.gro -p %s/prot.top -o test.tpr -maxwarn 50 -n %s  -quiet\n" % (pr_file,md_output_dir,md_output_dir, ndx_file)
			print comm1
			rfile.write(comm1)

			#make sure the *xtc file is the same one as used to generate the mdsnaps.pdb file
			comm1 = "$GROMACS_HOME/bin/mdrun -rerun %s/prot_md.xtc -s test.tpr -o test.trr -c test.gro -g test.log -e test.edr -quiet\n" % (md_output_dir)
	#		comm1 = "mpirun -np $ncpus $GROMACS_HOME/bin/mdrun_mpi -rerun %s/prot_md.xtc -s test.tpr -o test.trr -c test.gro -g test.log -e test.edr\n" % (md_output_dir)
			print comm1
			rfile.write(comm1)
			rfile.write("rm #*\n")
			
			#print 'number of water',num_wtr[i]
			rfile.write("chmod a+x g_energy*.txt\n")
			for j in range(0,num_wtr[i]):
				wtr_index = int(water_list[j].split('_')[1])
				fin = open("%s/g_energy_%d.txt"%(E_out,j), "w")
				fin.write("$GROMACS_HOME/bin/g_energy -f test.edr -s test.tpr -sum -o Energy_output_%d -quiet<< EOF \n" % (wtr_index)) 
				for k in range(0,j):
					fin.write("Coul-SR:%s-%s\nLJ-SR:%s-%s\nLJ-LR:%s-%s\nCoul-14:%s-%s\nLJ-14:%s-%s\n"%(water_list[k],water_list[j], water_list[k],water_list[j], water_list[k],water_list[j], water_list[k],water_list[j], water_list[k],water_list[j]))
				#for k in range(j,num_wtr[i]):
				fin.write("Coul-SR:%s\nLJ-SR:%s\nLJ-LR:%s\nCoul-14:%s\nLJ-14:%s\n"%(water_list[j],water_list[j],water_list[j],water_list[j],water_list[j]))

				fin.write("\n")

				fin.close()
				# Need to change the permissons of this file so that it can be ran % run it
				rfile.write("./g_energy_%d.txt\n"%j)

		rfile.write("rm test* #*\n")
		rfile.close()


main(sys.argv[1:])
