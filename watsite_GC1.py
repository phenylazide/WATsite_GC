# coding:utf-8
import os,glob,sys,re

#MyBindingSite = 'MyBindingSite.pdb'
os.environ['MyBindingSite_pdbc'] = 'MyBindingSite.pdb'
os.environ['ligand_mol2'] = 'ligand.mol2'
os.environ['protein_pdb'] = 'protein.pdb'
estimate_enthalpy = 1


reduce_exe_dir = '/home/phzd/test/amber16/bin/to_be_dispatched'
ncpus = 3
watermodel = 'spc'
withligand = 1
ligcharge = 1

os.environ['WATSITEHOME'] = '/home/phzd/g09E/WATsite_GC'
os.environ['GROMACS_HOME'] = '/opt/gmx5.12'
#os.environ['GROMACS_HOME'] = '/opt/gmx2016.4' 
os.environ['AMBERHOME'] = '/home/phzd/test/amber16'
os.environ['PYMOL_EXE'] = '/home/phzd/soft/pymol/./pymol'




base = os.getcwd()
os.chdir(base)

#bashline = "python /home/phzd/g09E/WATsite_GC/bin/run_reduce.py -r {reduce_exe_dir} -p {myprotein}".format(reduce_exe_dir=reduce_exe_dir,myprotein=myprotein)
#bashline1_0 = "$GROMACS_HOME/bin/gmx pdb2gmx -ff amber99sb-ildn -f protein.pdb -o proteinH.pdb -p protein.top -ignh -water {watermodel}".format(watermodel=watermodel)
## using amber99sb-ildn, why not -o proteinH.gro ??


for f in ["protein.pdb","ligand.mol2","MyBindingSite.pdb"]:
	if not os.path.exists(f):	
		os.system("cp backup/{f} ./{f}".format(f=f))
	else: continue
	if not os.path.exists(f):raise Exception('error: no {} file was generate in pdb2gmx step'.format(f))

if not os.path.exists('backup'):
	os.mkdir('backup')
	if not os.path.exists('backup/ligand.mol2'):
		os.system("cp protein.pdb ligand.mol2 MyBindingSite.pdb backup/") # backup
bashline1_0 = "$GROMACS_HOME/bin/gmx pdb2gmx -f protein.pdb -o proteinH.pdb -p protein.top -ignh <<EOF\n6\n1\nEOF\n"
os.system(bashline1_0)

for f in ["proteinH.pdb","protein.top"]:
	if not os.path.exists(f):raise Exception('error: no {} file was generate in pdb2gmx step'.format(f))

bashline1_1 = "sed -i 's/HO4/SOL/g' proteinH.pdb" 
bashline1_2 = "sed -i 's/HOH/SOL/g' proteinH.pdb"
bashline1_3 = "cp proteinH.pdb prot_only.pdb"
os.system(bashline1_1 + " && " + bashline1_2 + " && " + bashline1_3)

for f in ["prot_only.pdb"]:
	if not os.path.exists(f):raise Exception('error: no {} file was generate in pdb2gmx step'.format(f))
print('******* gmx pdb2gmx has done ************')

# 2b. prepare ligand. 
#  result in two major file: proteinH.pdb and protein.top which are included lig
if withligand == 1:	
	print('begining prepare ligand, will merge ligand coordinate and top itp to protein.pdb and protein.top')
	bashline2 = "python /home/phzd/g09E/WATsite_GC/bin/rename_SUB_GC.py ligand.mol2"
	os.system(bashline2)
	if (not os.path.exists('ligand.acpype/ligand_GMX.itp')) or (not os.path.exists('ligand.acpype/ligand_GMX.gro')):
		bashline2_1 = "python /home/phzd/g09E/WATsite_GC/bin/acpype.py -i ligand.mol2 -n {ligcharge}".format(ligcharge=ligcharge)
		os.system("source /home/phzd/test/amber16/amber.sh && "+ bashline2_1)
	if os.path.exists('ligand.acpype/ligand_GMX.itp') and os.path.exists('ligand.acpype/ligand_GMX.gro'):
		bashline2_2 = "cp ligand.acpype/ligand_GMX.itp ligand.acpype/ligand_GMX.gro ./"
		os.system(bashline2_2)
	else:raise IOError("acpype run failed... Please check acpype error...")

	bashline2_3 = "python /home/phzd/g09E/WATsite_GC/bin/mod_SUB_GC.py ligand_GMX.itp"
	## mod_SUB_GC.py  was modified by GC
	os.system(bashline2_3)

  #generate ligand position restraint file
	bashline2b = "echo '0' |$GROMACS_HOME/bin/gmx genrestr -f ligand_GMX.gro -o posre_ligand.itp"
	os.system(bashline2b)
	for f in ["posre_ligand.itp"]:
		if not os.path.exists(f):raise Exception('error: no {} file was generate in ligand gmx genrestr step'.format(f))

  # fuse ligand and protein files
	bashline2c = "mv proteinH.pdb prot_only.pdb"

	os.system("$GROMACS_HOME/bin/gmx editconf -f ligand_GMX.gro -o ligand_GMX.pdb")
	os.system("cat ligand_GMX.pdb prot_only.pdb > proteinH.pdb")
	os.system("sed -i '/TITLE/d' proteinH.pdb")
	os.system("sed -i '/MODEL/d' proteinH.pdb") 
	os.system("sed -i '/ENDMDL/d' proteinH.pdb") 
	os.system("python /home/phzd/g09E/WATsite_GC/bin/merge_ProtTop_GC.py") 
  ## after merge_ProtTop_GC.py, protein.top included ligand
	f = open('ligand_GMX.itp','a')
	f.write('\n\n\n; Include Position restraint file \n#ifdef POSRES                     \n#include "posre_ligand.itp"          \n#endif                            \n\n\n')
	f.close()
## 
	print('job of preparing ligand has done, coordinate and top have been merged to proteinH.pdb and protein.top')


# 2c. solvate and neutralize system

print('beginning 2c: solvate and neutralize system') 
if not os.path.exists('proteinH.pdb'):
	raise IOError("no proteinH.pdb was found for solvate and neutralize step")

if not os.path.exists("protein.gro"):
	bashline2c = "$GROMACS_HOME/bin/gmx editconf -bt cubic -f proteinH.pdb -o protein.gro -c -d 1.0"
	os.system(bashline2c)
for f in ["protein.gro"]:
	if not os.path.exists(f):raise Exception('error: no {} file was generate bashline2c gmx editconf -bt cubic step'.format(f))

#bashline2c_1 = "$GROMACS_HOME/bin/gmx genbox -cp protein.gro -cs spc216.gro -o prot_4ion.gro -p protein.top" 
### genbox:This tool has been split to gmx solvate and gmx insert-molecules.
if not os.path.exists("prot_4ion.gro"):
	bashline2c_1 = "$GROMACS_HOME/bin/gmx solvate -cp protein.gro -o prot_4ion.gro -p protein.top"
	os.system(bashline2c_1)
for f in ["prot_4ion.gro"]:
	if not os.path.exists(f):raise Exception('error: no {} file was generate in bashline2c_1 step'.format(f))

if not os.path.exists("prot_4ion.tpr"):
	bashline2c_2 = "$GROMACS_HOME/bin/gmx grompp -f /home/phzd/g09E/WATsite_GC/mdp_files/posre/em.mdp -c prot_4ion.gro -p protein.top -o prot_4ion.tpr -maxwarn 50"
	os.system(bashline2c_2)
for f in ["prot_4ion.tpr"]:
	if not os.path.exists(f):raise Exception('error: no {} file was generate in bashline2c_2 step'.format(f))

#get the ion charge from previous step
if not os.path.exists("prot_em.gro"):
	bashline2c_3 = "echo 'SOL' | $GROMACS_HOME/bin/gmx genion -s prot_4ion.tpr -o prot_em.gro -p protein.top -pname NA -nname CL -neutral"
	os.system(bashline2c_3)
for f in ["prot_em.gro"]:
	if not os.path.exists(f):raise Exception('error: no {} file was generate in bashline2c_3 step'.format(f))

# 3. run simulation

# 3a. minimization
print('now beginning minimization ')
if not os.path.exists("prot_min.trr"):
	bashline3_0 = "$GROMACS_HOME/bin/gmx grompp -f /home/phzd/g09E/WATsite_GC/mdp_files/posre/em.mdp -c prot_em.gro -p protein.top -o prot_min.tpr -maxwarn 50"
	bashline3_1 = "$GROMACS_HOME/bin/gmx mdrun -v -s prot_min.tpr -o prot_min.trr -c prot_pr.gro -g min.log -e min.edr"
	bashline3_2 = "rm \#*"
	os.system(bashline3_0 + " && " +  bashline3_1 + " && " + bashline3_2)

for f in ["prot_min.trr","min.log"]:
	if not os.path.exists(f):raise Exception('error: no {} file was generate in minimization step'.format(f))


# minimization check point
if os.path.exists('min.log') and (not os.path.exists('prot_pr.gro')):
	txt = open('min.log','r').read();	pattern = re.compile('no domain decomposition')
	match = re.findall(pattern, txt)
	if match:    ## using one cpu only
		os.system("$GROMACS_HOME/bin/gmx mdrun -v -s prot_min.tpr -o prot_min.trr -c prot_pr.gro -g min.log -e min.edr -nt 1")
	os.system("rm \#*")
 

# 3b. md, for getting prot_md.trr

print('now beginning pr of md ')
bashline3b = "$GROMACS_HOME/bin/gmx grompp -f /home/phzd/g09E/WATsite_GC/mdp_files/posre/pr.mdp -c prot_pr.gro -p protein.top -o prot_pr.tpr -maxwarn 50"
bashline3b_1 = "$GROMACS_HOME/bin/gmx mdrun -v -s prot_pr.tpr -o prot_pr.trr -c prot_equ.gro -g pr.log -e pr.edr"
os.system(bashline3b + " && " + bashline3b_1 + " && " + "rm \#*")


print('now beginning equ of md ')
bashline3b_2 = "$GROMACS_HOME/bin/gmx grompp -f /home/phzd/g09E/WATsite_GC/mdp_files/posre/equ.mdp -c prot_equ.gro -p protein.top -o prot_equ.tpr -maxwarn 50"
bashline3b_3 = "$GROMACS_HOME/bin/gmx mdrun -v -s prot_equ.tpr -o prot_equ.trr -c prot_md.gro -g equ.log -e equ.edr"
os.system(bashline3b_2 + " && " + bashline3b_3 + " && " + "rm \#*")

bashline3b_4 = "$GROMACS_HOME/bin/gmx grompp -f /home/phzd/g09E/WATsite_GC/mdp_files/posre/md_new.mdp -c prot_md.gro -p protein.top -o prot_md.tpr -maxwarn 50"
os.system(bashline3b_4)

if not os.path.exists('prot_md.tpr'):
	bashline3b_4 = "$GROMACS_HOME/bin/gmx grompp -f /home/phzd/g09E/WATsite_GC/mdp_files/posre/md.mdp -c prot_md.gro -p protein.top -o prot_md.tpr -maxwarn 50"
	os.system(bashline3b_4)

bashline3b_5 = "$GROMACS_HOME/bin/gmx mdrun -v -s prot_md.tpr -o prot_md.trr -c prot_fin.gro -g md.log -e md.edr"
os.system(bashline3b_5 + " && " + "rm \#*")


