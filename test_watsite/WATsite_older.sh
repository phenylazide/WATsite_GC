#!/bin/tcsh

setenv myligand 'MyBindingSite.pdb'
setenv myligandmol2 'MyLigand.mol2'
setenv myprotein 'MyProtein.pdb'
setenv estimate_enthalpy 1
# WATSITEHOME='/home/phzd/soft/WATsite2.0'
python $WATSITEHOME/bin/run_reduce.py -r $reduce_exe_dir -p $myprotein

if ( ! -e prot_em.gro ) then 
  if ( $watermodel != 'tip4pew' ) then
    $GROMACS_HOME/bin/gmx pdb2gmx -ff amber99sb-ildn -f prot.pdb -o protH.pdb -p prot.top -water $watermodel >&/dev/null
  else
    echo '3' | $GROMACS_HOME/bin/gmx pdb2gmx -ff amber99sb-ildn -f prot.pdb -o protH.pdb -p prot.top >&/dev/null
  endif
  sed -i 's/HO4/SOL/g' protH.pdb 
  sed -i 's/HOH/SOL/g' protH.pdb 
  if ( ! -e protH.pdb ) then
    echo "failed to add hydrogen by gromacs"
    echo "try amber prep now..."
    #sed -i '/HOH/d' prot.pdb
    #sed -i '/HO4/d' prot.pdb
    $AMBERHOME/bin/tleap -f $AMBERHOME/dat/leap/cmd/oldff/leaprc.ff99SBildn -f $WATSITEHOME/dat/addmissingAtom.leapP.in >&/dev/null
    if ( $watermodel != 'tip4pew' ) then
      $GROMACS_HOME/bin/gmx pdb2gmx -ff amber99sb-ildn -f prot_amber.pdb -o protH.pdb -p prot.top -water $watermodel -ignh -chainsep id >&/dev/null
    else
      echo '3' | $GROMACS_HOME/bin/gmx pdb2gmx -ff amber99sb-ildn -f prot_amber.pdb -o protH.pdb -p prot.top -ignh -chainsep id >&/dev/null
    endif
    sed -i 's/HO4/SOL/g' protH.pdb 
    sed -i 's/HOH/SOL/g' protH.pdb
  endif
endif

if ( ! -e protH.pdb || ! -e prot.top ) then
  echo "protein preparation has failed. exit..."
  exit 1
endif
cp protH.pdb prot_only.pdb

# 2b. prepare ligand
if ( $withligand == 1) then
  python $WATSITEHOME/bin/rename_SUB.py MyLigand.mol2
  if ( ! -e lig.acpype/lig_GMX.itp || ! -e lig.acpype/lig_GMX.gro ) then
    python $WATSITEHOME/bin/acpype.py -i lig.mol2 -o gmx -n $ligcharge
  if ( -e lig.acpype/lig_GMX.itp && -e lig.acpype/lig_GMX.gro ) then
    cp lig.acpype/lig_GMX.itp lig.acpype/lig_GMX.gro .
  else
    echo "\n\nAcpype run failed... Please check acpype error...\n\n"
    exit
  endif
  python $WATSITEHOME/bin/mod_SUB.py lig_GMX.itp

  #generate ligand position restraint file
  echo '0' |$GROMACS_HOME/bin/gmx genrestr -f lig_GMX.gro -o posre_lig.itp >&/dev/null

  # fuse ligand and protein files
  mv protH.pdb prot_only.pdb
  $GROMACS_HOME/bin/gmx editconf -f lig_GMX.gro -o lig_GMX.pdb
  cat lig_GMX.pdb prot_only.pdb > protH.pdb
  sed -i '/TITLE/d' protH.pdb 
  sed -i '/MODEL/d' protH.pdb 
  sed -i '/ENDMDL/d' protH.pdb 
  python $WATSITEHOME/bin/merge_ProtTop.py

  echo '' >> lig_GMX.itp
  echo '' >> lig_GMX.itp
  echo "; Include Position restraint file " >> lig_GMX.itp
  echo "#ifdef POSRES                     " >> lig_GMX.itp
  echo '#include "posre_lig.itp"          ' >> lig_GMX.itp
  echo '#endif                            ' >> lig_GMX.itp
  echo '' >> lig_GMX.itp
  echo '' >> lig_GMX.itp

endif #endif prepare ligand


# 2c. solvate and neutralize system
if ( -e protH.pdb ) then
   $GROMACS_HOME/bin/gmx editconf -bt cubic -f protH.pdb -o prot.gro -c -d 1.0 >&/dev/null
   if ( $watermodel == 'spc' || $watermodel == 'spce' ) then
      #$GROMACS_HOME/bin/gmx genbox -cp prot.gro -cs spc216.gro -o prot_4ion.gro -p prot.top >&/dev/null
  ### genbox:This tool has been split to gmx solvate and gmx insert-molecules.
      $GROMACS_HOME/bin/gmx solvate -cp prot.gro -cs spc216.gro -o prot_4ion.gro -p prot.top >&/dev/null
   else if ( $watermodel == 'tip4p' || $watermodel == 'tip4pew' ) then
      $GROMACS_HOME/bin/gmx genbox -cp prot.gro -cs tip4p.gro -o prot_4ion.gro -p prot.top >&/dev/null
   endif
   if ( ! -e prot_4ion.gro ) then
      if ( $watermodel == 'spc' || $watermodel == 'spce' ) then
         $GROMACS_HOME/bin/gmx solvate -cp prot.gro -cs spc216.gro -o prot_4ion.gro -p prot.top >&/dev/null
      else if ( $watermodel == 'tip4p' || $watermodel == 'tip4pew' ) then
         $GROMACS_HOME/bin/gmx solvate -cp prot.gro -cs tip4p.gro -o prot_4ion.gro -p prot.top >&/dev/null
      endif
   endif
   $GROMACS_HOME/bin/gmx grompp -f $WATSITEHOME/mdp_files/posre/em.mdp -c prot_4ion.gro -p prot.top -o prot_4ion.tpr -maxwarn 50 >& charge.log
   #get the ion charge from previous step
   set charge=`python $WATSITEHOME/bin/getCharge.py`
   if ( $charge > 0 ) then
      echo 'SOL' | $GROMACS_HOME/bin/gmx genion -s prot_4ion.tpr -o prot_em.gro -nname CL -nn $charge -p prot.top >&/dev/null
   else if ( $charge < 0 ) then
      set abscharge=`expr $charge \* -1`
      echo 'SOL' | $GROMACS_HOME/bin/gmx genion -s prot_4ion.tpr -o prot_em.gro -nname NA -nn $abscharge -p prot.top >&/dev/null
   else
      cp prot_4ion.gro prot_em.gro
   endif 
else 
   echo "failed to add hydrogens"
   exit 2
endif
endif

# 3. run simulation

if ($ncpus == 1) then
  # 3a. minimization
  if ( ! -e prot_pr.gro ) then 
    $GROMACS_HOME/bin/gmx grompp -f $WATSITEHOME/mdp_files/posre/em.mdp -c prot_em.gro -p prot.top -o prot_min.tpr -maxwarn 50 >&/dev/null
    $GROMACS_HOME/bin/gmx mdrun -s prot_min.tpr -o prot_min.trr -c prot_pr.gro -g min.log -e min.edr
    rm \#*
  endif
  # minimization check point
  if ( -e min.log && ! -e prot_pr.gro ) then
    set exist = `grep -c "no domain decomposition" min.log`
    if ( $exist > 0 ) then
      echo "test3"
      $GROMACSHOME/bin/mdrun -s prot_min.tpr -o prot_min.trr -c prot_pr.gro -g min.log -e min.edr -nt 1
      rm \#*
    endif
  endif   
  # 3b. md
  if ( ! -e prot_equ.gro ) then
    $GROMACS_HOME/bin/gmx grompp -f $WATSITEHOME/mdp_files/posre/pr.mdp -c prot_pr.gro -p prot.top -o prot_pr.tpr -maxwarn 50 >&/dev/null
    $GROMACS_HOME/bin/gmx mdrun -s prot_pr.tpr -o prot_pr.trr -c prot_equ.gro -g pr.log -e pr.edr
    rm \#*
  endif
  if ( ! -e prot_md.gro ) then
    $GROMACS_HOME/bin/gmx grompp -f $WATSITEHOME/mdp_files/posre/equ.mdp -c prot_equ.gro -p prot.top -o prot_equ.tpr -maxwarn 50 >&/dev/null
    $GROMACS_HOME/bin/gmx mdrun -s prot_equ.tpr -o prot_equ.trr -c prot_md.gro -g equ.log -e equ.edr
    rm \#*
  endif
  if ( ! -e prot_fin.gro ) then
    if ( -e $WATSITEHOME/mdp_files/posre/md_new.mdp ) then
      $GROMACS_HOME/bin/gmx grompp -f $WATSITEHOME/mdp_files/posre/md_new.mdp -c prot_md.gro -p prot.top -o prot_md.tpr -maxwarn 50 >&/dev/null
    else
      $GROMACS_HOME/bin/gmx grompp -f $WATSITEHOME/mdp_files/posre/md.mdp -c prot_md.gro -p prot.top -o prot_md.tpr -maxwarn 50 >&/dev/null
    endif
    $GROMACS_HOME/bin/gmx mdrun -s prot_md.tpr -o prot_md.trr -c prot_fin.gro -g md.log -e md.edr
    rm \#*
  endif
  rm \#*
else if ( $ncpus > 1 ) then
  # 3a. minimization
  if ( ! -e prot_pr.gro ) then 
    $GROMACS_HOME/bin/gmx grompp -f $WATSITEHOME/mdp_files/posre/em.mdp -c prot_em.gro -p prot.top -o prot_min.tpr -maxwarn 50 >&/dev/null
    mpirun -np $ncpus $GROMACS_HOME/bin/gmx mdrun_mpi -s prot_min.tpr -o prot_min.trr -c prot_pr.gro -g min.log -e min.edr
    rm \#*
  endif
  # minimization check point
  if ( -e min.log && ! -e prot_pr.gro ) then
    set exist = `grep -c "no domain decomposition" min.log`
    if ( $exist > 0 ) then
      echo "test3"
      mpirun -np $ncpus $GROMACS_HOME/bin/gmx mdrun_mpi -s prot_min.tpr -o prot_min.trr -c prot_pr.gro -g min.log -e min.edr -nt 1
      rm \#*
    endif
  endif   
 # 3b. md
  if ( ! -e prot_equ.gro ) then
    $GROMACS_HOME/bin/gmx grompp -f $WATSITEHOME/mdp_files/posre/pr.mdp -c prot_pr.gro -p prot.top -o prot_pr.tpr -maxwarn 50 >&/dev/null
    mpirun -np $ncpus $GROMACS_HOME/bin/gmx mdrun_mpi -s prot_pr.tpr -o prot_pr.trr -c prot_equ.gro -g pr.log -e pr.edr
    rm \#*
  endif
  if ( ! -e prot_md.gro ) then
    $GROMACS_HOME/bin/gmx grompp -f $WATSITEHOME/mdp_files/posre/equ.mdp -c prot_equ.gro -p prot.top -o prot_equ.tpr -maxwarn 50 >&/dev/null
    mpirun -np $ncpus $GROMACS_HOME/bin/gmx mdrun_mpi -s prot_equ.tpr -o prot_equ.trr -c prot_md.gro -g equ.log -e equ.edr
    rm \#*
  endif
  if ( ! -e prot_fin.gro ) then
    if ( -e $WATSITEHOME/mdp_files/posre/md_new.mdp ) then
      $GROMACS_HOME/bin/gmx grompp -f $WATSITEHOME/mdp_files/posre/md_new.mdp -c prot_md.gro -p prot.top -o prot_md.tpr -maxwarn 50 >&/dev/null
    else
      $GROMACS_HOME/bin/gmx grompp -f $WATSITEHOME/mdp_files/posre/md.mdp -c prot_md.gro -p prot.top -o prot_md.tpr -maxwarn 50 >&/dev/null
    endif
    mpirun -np $ncpus $GROMACS_HOME/bin/gmx mdrun_mpi -s prot_md.tpr -o prot_md.trr -c prot_fin.gro -g md.log -e md.edr
    rm \#*
  endif
  rm \#*
endif

# 4. post processing 
if ( -e prot_md.trr && ! -e mdsnaps.pdb ) then
   #echo 'Protein' | trjconv -f prot_md.gro -s prot_md.tpr -o ref.pdb
   echo 'System' | $GROMACS_HOME/bin/gmx trjconv -f prot_md.trr -s prot_md.tpr -o mdsnaps.pdb  >&/dev/null
   $GROMACS_HOME/bin/gmx trjconv -f prot_md.trr -o prot_md.xtc
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
(echo 0; echo 'q';) | $GROMACS_HOME/bin/gmx make_ndx -f prot_em.gro -o index.ndx >&/dev/null
   
# 5. move 
echo 'Protein' | $GROMACS_HOME/bin/gmx trjconv -f prot_md.gro -s prot_md.tpr -o ref.pdb >&/dev/null
python $WATSITEHOME/bin/reformat_MyBindingSite.py
$PYMOL_EXE -c $WATSITEHOME/bin/shift_ligand.pml 
rm prot_only_MyBindingSite.pdb MyBindingSite_formatted.pdb

# 6. calculate hydration sites
$WATSITEHOME/bin/hydroentropy -c $WATSITEHOME/dat/HydroEntropy_mod.bcf -o OUTPUT 
if ( ! -e OUTPUT/HydrationSites.mol2 ) then
   echo "fail to generate hydration sites"
   exit 4
endif

# 6a. estimate enthalpy
if ( $estimate_enthalpy == 1 ) then
   python $WATSITEHOME/bin/estimate_enthalpy.py -i OUTPUT/wtrindex_folder/ 
   chmod 777 SuperRunScript
   ./SuperRunScript
   rm ENTHALPY_OUTPUT -r
   mkdir ENTHALPY_OUTPUT
   echo "mv files"
   cp Enthalpy_output_*/*.xvg ENTHALPY_OUTPUT/
   $WATSITEHOME/bin/hydroenthalpy -c $WATSITEHOME/dat/HydroEnthalpy.bcf
   cd ENTHALPY_OUTPUT
   tar -cf Energy_output.tar *.xvg
   rm *.xvg
   cd ../

   #shift hydration site
   mv HydrationSites.pdb md_HydrationSites.pdb
   mv HydrationSites.mol2 md_HydrationSites.mol2
   #shift based on full protein
   $PYMOL_EXE -c $WATSITEHOME/bin/shift_hydrationsites.pml
   python $WATSITEHOME/bin/pdb2mol2.py -i HydrationSites.pdb -o HydrationSites.mol2
   rm ref_md_HydrationSites.pdb

   #shift based on binding site residues (w. 5A to MyBindingSite.pdb)
   $PYMOL_EXE -c $WATSITEHOME/bin/shift_hydrationsites_byBS.pml
   python $WATSITEHOME/bin/pdb2mol2.py -i HydrationSites_byBS.pdb -o HydrationSites_byBS.mol2
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

