base=`pwd`
cd $base

export WATSITEHOME='/home/phzd/g09E/WATsite_GC'
export MyBindingSite_pdb='MyBindingSite.pdb'
export estimate_enthalpy=1
export GROMACS_HOME='/opt/gmx2016.4'
myligandmol2='ligand.mol2'
myprotein='protein.pdb'
PYMOL_EXE='/home/phzd/2_software/pymol2/./pymol'


# 4. post processing 
if [ -e prot_md.trr ]&&[ ! -e mdsnaps.pdb ];then
   #echo 'Protein' | trjconv -f prot_md.gro -s prot_md.tpr -o ref.pdb
   echo 'System' | /opt/gmx2016.4/bin/gmx trjconv -f prot_md.trr -s prot_md.tpr -o mdsnaps.pdb  
   /opt/gmx2016.4/bin/gmx trjconv -f prot_md.trr -o prot_md.xtc
   #rm prot_min.trr prot_pr.trr prot_equ.trr prot_md.trr
   rm min.log pr.log equ.log md.log
   rm *.edr
   rm \#*
   rm ANTECHAMBER* 
   rm step*
elif [ -e mdsnaps.pdb ]; then
   echo "mdsnaps.pdb generated"
else
   echo "md simulation failed"
   exit 4
fi
(echo '0'; echo 'q';) |/opt/gmx2016.4/bin/gmx make_ndx -f prot_em.gro -o index.ndx
echo 'successful done # 4. post processing '


# 5. move 

export MyBindingSite_pdb='MyBindingSite.pdb'

echo 'Protein' | /opt/gmx2016.4/bin/gmx trjconv -f prot_md.gro -s prot_md.tpr -o ref.pdb 
python /home/phzd/g09E/WATsite_GC/bin/reformat_MyBindingSite_GC.py
/home/phzd/2_software/pymol2/./pymol -c /home/phzd/g09E/WATsite_GC/bin/shift_ligand.pml 
rm prot_only_MyBindingSite.pdb MyBindingSite_formatted.pdb

# 6. calculate hydration sites
/home/phzd/g09E/WATsite_GC/bin/hydroentropy -c /home/phzd/g09E/WATsite_GC/dat/HydroEntropy_mod.bcf -o OUTPUT 
if [ ! -e OUTPUT/HydrationSites.mol2 ];then
   echo "fail to generate hydration sites"
   exit 4
fi
echo 'succesful done # 5. move '

# 6a. estimate enthalpy
source activate python27
export GROMACS_HOME='/opt/gmx2016.4'

if [ $estimate_enthalpy == 1 ]; then
   python /home/phzd/g09E/WATsite_GC/bin/estimate_enthalpy_GC.py -i OUTPUT/wtrindex_folder/ 
   ##chmod 777 SuperRunScript
   ###./SuperRunScript    ### this is step is big keng!!!!!!!!!!!!!!!, one run will generate 44 GB, let alone bunch of file lists. 
fi
