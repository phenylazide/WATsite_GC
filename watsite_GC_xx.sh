
export WATSITEHOME='/home/phzd/g09E/WATsite_GC'
export MyBindingSite_pdb='MyBindingSite.pdb'
export estimate_enthalpy=1
#export GROMACS_HOME='/opt/gmx2016.4'
export GROMACS_HOME='/opt/gmx5.12'
export ligand_mol2='ligand.mol2'
export protein_pdb='protein.pdb'
export PYMOL_EXE='/home/phzd/2_software/pymol2/./pymol'


# 6a. estimate enthalpy
if ( $estimate_enthalpy == 1 ) then
   #python /home/phzd/g09E/WATsite_GC/bin/estimate_enthalpy.py -i OUTPUT/wtrindex_folder/ 
   chmod 777 SuperRunScript
   ./SuperRunScript
   rm ENTHALPY_OUTPUT -r
   mkdir ENTHALPY_OUTPUT
   echo "mv files"
   cp Enthalpy_output_*/*.xvg ENTHALPY_OUTPUT/
   /home/phzd/g09E/WATsite_GC/bin/hydroenthalpy -c /home/phzd/g09E/WATsite_GC/dat/HydroEnthalpy.bcf
   cd ENTHALPY_OUTPUT
   tar -cf Energy_output.tar *.xvg
   rm *.xvg
   cd ../

   #shift hydration site
   mv HydrationSites.pdb md_HydrationSites.pdb
   mv HydrationSites.mol2 md_HydrationSites.mol2
   #shift based on full protein
   $PYMOL_EXE -c /home/phzd/g09E/WATsite_GC/bin/shift_hydrationsites.pml
   python /home/phzd/g09E/WATsite_GC/bin/pdb2mol2.py -i HydrationSites.pdb -o HydrationSites.mol2
   rm ref_md_HydrationSites.pdb

   #shift based on binding site residues (w. 5A to MyBindingSite.pdb)
   $PYMOL_EXE -c /home/phzd/g09E/WATsite_GC/bin/shift_hydrationsites_byBS.pml
   python /home/phzd/g09E/WATsite_GC/bin/pdb2mol2.py -i HydrationSites_byBS.pdb -o HydrationSites_byBS.mol2
   rm bs_res_ref_HydrationSites.pdb bs_res.pdb 

   rm myjob.pbsid*
   rm Enthalpy_output_* -r
endif

# 7. clean up the files
if ( -e HydrationSites.mol2 && -e cluster.egy ) then
   #rm *.itp *.top *.gro *.tpr *.xtc *.log *.cpt *.trr RunScript_* SuperRunScript
   #rm mdsnaps.pdb protH.pdb prot.pdb lig.pdb prot_amber.pdb tmp_shift.pdb
   #rm mdsnaps.pdb *trr
   rm Enthalpy_* 
   rm \#* 
endif 

	echo "Protein             " $myprotein > WATsite.out
	echo "Ligand              " $myligandmol2 >> WATsite.out
	echo "HydrationSiteMol2   HydrationSites.mol2" >> WATsite.out
	echo "HydrationSitePDB    HydrationSites.pdb" >> WATsite.out
	echo "WaterTraj           OUTPUT/WATinside.mol2" >> WATsite.out
	echo "EnergyFile          cluster.egy" >> WATsite.out

