#!/bin/tcsh

setenv myligand 'MyBindingSite.pdb'
setenv myligandmol2 'MyLigand.mol2'
setenv myprotein 'MyProtein.pdb'
setenv estimate_enthalpy 1

# 4. post processing 
if ( -e prot_md.trr && ! -e mdsnaps.pdb ) then
   #echo 'Protein' | trjconv -f prot_md.gro -s prot_md.tpr -o ref.pdb
   echo 'System' | /opt/gmx2016.4/bin/gmx trjconv -f prot_md.trr -s prot_md.tpr -o mdsnaps.pdb  
   /opt/gmx2016.4/bin/gmx trjconv -f prot_md.trr -o prot_md.xtc
   #rm prot_min.trr prot_pr.trr prot_equ.trr prot_md.trr
   rm min.log pr.log equ.log md.log
   rm *.edr
   rm \#*
   rm ANTECHAMBER* 
   rm step*
else if ( -e mdsnaps.pdb ) then
   echo "mdsnaps.pdb generated"
else
   echo "md simulation failed"
   exit 4
endif
(echo 0; echo 'q';) | /opt/gmx2016.4/bin/gmx make_ndx -f prot_em.gro -o index.ndx 
   
# 5. move 
echo 'Protein' | /opt/gmx2016.4/bin/gmx trjconv -f prot_md.gro -s prot_md.tpr -o ref.pdb 
python /home/phzd/g09E/WATsite_GC/bin/reformat_MyBindingSite.py
$PYMOL_EXE -c /home/phzd/g09E/WATsite_GC/bin/shift_ligand.pml 
rm prot_only_MyBindingSite.pdb MyBindingSite_formatted.pdb

# 6. calculate hydration sites
/home/phzd/g09E/WATsite_GC/bin/hydroentropy -c /home/phzd/g09E/WATsite_GC/dat/HydroEntropy_mod.bcf -o OUTPUT 
if ( ! -e OUTPUT/HydrationSites.mol2 ) then
   echo "fail to generate hydration sites"
   exit 4
endif

# 6a. estimate enthalpy
if ( $estimate_enthalpy == 1 ) then
   python /home/phzd/g09E/WATsite_GC/bin/estimate_enthalpy.py -i OUTPUT/wtrindex_folder/ 
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

