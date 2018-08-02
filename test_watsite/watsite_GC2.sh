# 4. post processing , for getting mdsnaps.pdb


echo 'System' | /opt/gmx2016.4/bin/gmx trjconv -f prot_md.trr -s prot_md.tpr -o mdsnaps.pdb  

/opt/gmx2016.4/bin/gmx trjconv -f prot_md.trr -o prot_md.xtc
#rm prot_min.trr prot_pr.trr prot_equ.trr prot_md.trr
rm min.log pr.log equ.log md.log
rm *.edr
rm \#*
rm ANTECHAMBER* 
rm step*
(echo 0; echo 'q';) | /opt/gmx2016.4/bin/gmx make_ndx -f prot_em.gro -o index.ndx 


 
# 5. move 
echo 'Protein' | /opt/gmx2016.4/bin/gmx trjconv -f prot_md.gro -s prot_md.tpr -o ref.pdb 
python /home/phzd/soft/WATsite2.0/bin/reformat_MyBindingSite.py
$PYMOL_EXE -c /home/phzd/soft/WATsite2.0/bin/shift_ligand.pml 
rm prot_only_MyBindingSite.pdb MyBindingSite_formatted.pdb

# 6. calculate hydration sites
/home/phzd/soft/WATsite2.0/bin/hydroentropy -c /home/phzd/soft/WATsite2.0/dat/HydroEntropy_mod.bcf -o OUTPUT 
if not os.path.exists("OUTPUT/HydrationSites.mol2"): 
	raise IOError("fail to generate hydration sites")


# 6a. estimate enthalpy

python /home/phzd/soft/WATsite2.0/bin/estimate_enthalpy.py -i OUTPUT/wtrindex_folder/ 
chmod 777 SuperRunScript
./SuperRunScript
rm ENTHALPY_OUTPUT -r
mkdir ENTHALPY_OUTPUT
echo "mv files"
cp Enthalpy_output_*/*.xvg ENTHALPY_OUTPUT/
/home/phzd/soft/WATsite2.0/bin/hydroenthalpy -c /home/phzd/soft/WATsite2.0/dat/HydroEnthalpy.bcf
cd ENTHALPY_OUTPUT
tar -cf Energy_output.tar *.xvg
rm *.xvg
cd ../

#shift hydration site
mv HydrationSites.pdb md_HydrationSites.pdb
mv HydrationSites.mol2 md_HydrationSites.mol2
#shift based on full protein
$PYMOL_EXE -c /home/phzd/soft/WATsite2.0/bin/shift_hydrationsites.pml
python /home/phzd/soft/WATsite2.0/bin/pdb2mol2.py -i HydrationSites.pdb -o HydrationSites.mol2
rm ref_md_HydrationSites.pdb

#shift based on binding site residues (w. 5A to MyBindingSite.pdb)
$PYMOL_EXE -c /home/phzd/soft/WATsite2.0/bin/shift_hydrationsites_byBS.pml
python /home/phzd/soft/WATsite2.0/bin/pdb2mol2.py -i HydrationSites_byBS.pdb -o HydrationSites_byBS.mol2
rm bs_res_ref_HydrationSites.pdb bs_res.pdb 

rm myjob.pbsid*
rm Enthalpy_output_* -r


# 7. clean up the files
if ( -e HydrationSites.mol2 && -e cluster.egy ) then
   #rm *.itp *.top *.gro *.tpr *.xtc *.log *.cpt *.trr RunScript_* SuperRunScript
   #rm mdsnaps.pdb proteinH.pdb protein.pdb ligand.pdb prot_amber.pdb tmp_shift.pdb
   #rm mdsnaps.pdb *trr
   rm Enthalpy_* 
   rm \#* 
endif 

print("Protein    " + myprotein)
print("Ligand  " + myligandmol2)
print("HydrationSiteMol2   HydrationSites.mol2")
print("HydrationSitePDB    HydrationSites.pdb")
print("WaterTraj           OUTPUT/WATinside.mol2")
print("EnergyFile          cluster.egy")

