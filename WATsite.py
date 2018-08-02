#
# WATsite.py
#
from Tkinter import *
from tkFileDialog import *
from pymol import cmd, stored, util
from pymol.wizard import Wizard
from ScrolledText import *

import tkSimpleDialog
import tkMessageBox
import os
import sys
import math
import shutil
import subprocess
import colorsys,re
import Tkinter,Pmw
import Pmw

monitor_list = []

delta = u"\u0394" #unicode (somehow does not work on cheeta)
#delta = 'd'



#######################################################################################
def __init__(self):
    cmd.set("retain_order", 1)
    cmd.set("pdb_use_ter_records", 1)

    # Simply add the menu entry and callback
    self.menuBar.addmenu('WATsite', 'Hydration Site Analysis',tearoff=TRUE)
    
    self.menuBar.addmenuitem('WATsite', 'command', 'Modify Paths and Parameters',
                label = 'Modify Program Location and Simulation Parameters',command = lambda s=self : modify_settings(s))
    #self.menuBar.addmenuitem('WATsite', 'command', 'Modify Parameters of MD Simulation',label = 'Modify parameters of simulation',command = lambda s=self : modify_parameters(s))
    self.menuBar.addmenuitem('WATsite', 'separator')

    self.menuBar.addmenuitem('WATsite', 'command', 'Prepare target',
                            label = 'Prepare system and submit hydration site analysis',command = lambda s=self : prepare_protein(s))
    self.menuBar.addmenuitem('WATsite', 'separator')

    self.menuBar.addmenuitem('WATsite', 'command', 'Import/Monitor Results',
                            label = 'Import Results',command = lambda s=self : import_results(s))
    self.menuBar.addmenuitem('WATsite', 'command', 'Estimate desolvation energy for ligand',
                            label = 'Estimate desolvation energy for ligand',command = lambda s=self : estimate_desolvation(s))



#######################################################################################
def read_settings():
    global username
    global WATsite_home
    global amber_home
    global gromacs_home
    global pymol_exe_dir
    global reduce_exe_dir
    global ncpus
    
    path1 = os.environ.get('HOME')
    filename = "%s/WATsite_Settings.txt" % (path1)
    
    try:
        fi = open(filename, "r")
    except:
        tkMessageBox.showwarning("Settings file not found", "%s not found. Please, copy download default file and copy to $HOME directory."% (filename) )
        sys.exit()

    for i in fi:
        if i.find("username") >= 0:
            tmpc, dat = i.split(None, 1)
            dat2 = str(dat)
            username = dat2.rstrip()

        if i.find("WATsite_home") >= 0:
            tmpc, dat = i.split(None, 1)
            dat2 = str(dat)
            WATsite_home = dat2.rstrip()
        if i.find("amber_home") >= 0:
            tmpc, dat = i.split(None, 1)
            dat2 = str(dat)
            amber_home = dat2.rstrip()
        if i.find("gromacs_home") >= 0:
            tmpc, dat = i.split(None, 1)
            dat2 = str(dat)
            gromacs_home = dat2.rstrip()
        if i.find("pymol_exe_dir") >= 0:
            tmpc, dat = i.split(None, 1)
            dat2 = str(dat)
            pymol_exe_dir = dat2.rstrip()
        if i.find("reduce_exe_dir") >= 0:
            tmpc, dat = i.split(None, 1)
            dat2 = str(dat)
            reduce_exe_dir = dat2.rstrip()
        if i.find("ncpus") >= 0:
            tmpc, dat = i.split(None, 1)
            ncpus = int(dat)
    fi.close()
    """
    print username
    print WATsite_home
    print amber_home
    print gromacs_home
    print pymol_exe_dir
    print reduce_exe_dir
    print ncpus
    """


#######################################################################################
class modifySettingsFile:
    def __init__(self, top):
        self.dialog = Pmw.Dialog(top,
                                buttons = ('Exit','Save','Reset to defaults'),
                                defaultbutton = 'Exit',
                                title = 'Modify Settings and Parameters',
                                command = self.apply)

        self.curdir = os.getcwd()
        master = self.dialog.interior()
        Tkinter.Label(master, text='Modify Paths to Installed Programs & Parameters for WATsite Analysis').pack(expand = 1, fill = 'both', padx = 5, pady = 5)

        nb1 = Pmw.NoteBook(master)

        # username page
        user_page = nb1.add('Username')

        Tkinter.Label(user_page, text='username:').grid(row=0, column=0, sticky=W, padx = 5, pady=1)
        self.e_username = Tkinter.Entry(user_page, width=70)
        self.e_username.insert(0, username)
        self.e_username.grid(row=0, column=1, sticky=W, padx = 5, pady=1)


        # MM page
        MM_page = nb1.add('Modify Settings')

        Tkinter.Label(MM_page, text='WATsite_home (WATsite directory):').grid(row=0, column=0, sticky=W, padx = 5, pady=1)
        self.e_WATsite_home = Tkinter.Entry(MM_page, width=70)
        self.e_WATsite_home.insert(0, WATsite_home)
        self.e_WATsite_home.grid(row=0, column=1, sticky=W, padx = 5, pady=1)
        b_WATsite_home = Tkinter.Button(MM_page, text='Browse', command = lambda: self.searchFOLDER(WATsite_home,self.e_WATsite_home)).grid(row=0, column=2, sticky=W, padx = 5, pady=1)

        Tkinter.Label(MM_page, text='amber_home (main amber directory = $AMBERHOME):').grid(row=1, column=0, sticky=W, padx = 5, pady=1)
        self.e_amber_home = Tkinter.Entry(MM_page, width=70)
        self.e_amber_home.insert(0, amber_home)
        self.e_amber_home.grid(row=1, column=1, sticky=W, padx = 5, pady=1)
        b_amber_home = Tkinter.Button(MM_page, text='Browse', command = lambda: self.searchFOLDER(amber_home,self.e_amber_home)).grid(row=1, column=2, sticky=W, padx = 5, pady=1)

        Tkinter.Label(MM_page, text='gromacs_home (main gromacs directory = $GROMACS_HOME):').grid(row=2, column=0, sticky=W, padx = 5, pady=1)
        self.e_gromacs_home = Tkinter.Entry(MM_page, width=70)
        self.e_gromacs_home.insert(0, gromacs_home)
        self.e_gromacs_home.grid(row=2, column=1, sticky=W, padx = 5, pady=1)
        b_gromacs_home = Tkinter.Button(MM_page, text='Browse', command = lambda: self.searchFOLDER(gromacs_home,self.e_gromacs_home)).grid(row=2, column=2, sticky=W, padx = 5, pady=1)
        
        Tkinter.Label(MM_page, text='pymol_exe_dir (directory containing pymol executable):').grid(row=4, column=0, sticky=W, padx = 5, pady=1)
        self.e_pymol_exe_dir = Tkinter.Entry(MM_page, width=70)
        self.e_pymol_exe_dir.insert(0, pymol_exe_dir)
        self.e_pymol_exe_dir.grid(row=4, column=1, sticky=W, padx = 5, pady=1)
        b_pymol_exe = Tkinter.Button(MM_page, text='Browse', command = lambda: self.searchFOLDER(pymol_exe_dir,self.e_pymol_exe_dir)).grid(row=4, column=2, sticky=W, padx = 5, pady=1)
        
        Tkinter.Label(MM_page, text='reduce_exe_dir (directory containing reduce executable):').grid(row=5, column=0, sticky=W, padx = 5, pady=1)
        self.e_reduce_exe_dir = Tkinter.Entry(MM_page, width=70)
        self.e_reduce_exe_dir.insert(0, reduce_exe_dir)
        self.e_reduce_exe_dir.grid(row=5, column=1, sticky=W, padx = 5, pady=1)
        b_reduce_exe_dir = Tkinter.Button(MM_page, text='Browse', command = lambda: self.searchFOLDER(reduce_exe_dir,self.e_reduce_exe_dir)).grid(row=5, column=2, sticky=W, padx = 5, pady=1)

        Tkinter.Label(MM_page, text='ncpus (number of processors):').grid(row=6, column=0, sticky=W, padx = 5, pady=1)
        self.e_ncpus = Tkinter.Entry(MM_page, width=70)
        self.e_ncpus.insert(0, ncpus)
        self.e_ncpus.grid(row=6, column=1, sticky=W, padx = 5, pady=1)


        # PM: ParaMeter page
        PM_page = nb1.add('Modify Parameters')

        Tkinter.Label(PM_page, text="Length of production MD run ", relief="sunken", width=100, pady=4).grid(row=0, columnspan=5, sticky=N)
        Tkinter.Label(PM_page, text="Simulation length (ns):").grid(row=1, column=0, sticky=W, padx = 5, pady=1)

        self.eMD = Entry(PM_page, width=70)
        self.eMD.insert(0, "4.0")
        self.eMD.grid(row=1, column=1)


        Tkinter.Label(PM_page, text="  ").grid(row=2, columnspan=5, sticky=N)
        Tkinter.Label(PM_page, text="Four choices of water models (spc, spce, tip4p, tip4pew)", relief="sunken", width=100, pady=4).grid(row=3, columnspan=5, sticky=N)
        Tkinter.Label(PM_page, text="Water model:").grid(row=4, column=0, sticky=W, padx = 5, pady=1)
        self.eWater = Entry(PM_page, width=70)
        self.eWater.insert(0, "spc")
        self.eWater.grid(row=4, column=1)


        Tkinter.Label(PM_page, text="  ").grid(row=5, columnspan=5, sticky=N)
        Tkinter.Label(PM_page, text="Two choices of clustering methods (QT or DBSCAN)", relief="sunken", width=100, pady=4).grid(row=6, columnspan=5, sticky=N)
        Tkinter.Label(PM_page, text="Clustering method:").grid(row=7, column=0, sticky=W, padx = 5, pady=1)
        self.eCluster = Entry(PM_page, width=70)
        self.eCluster.insert(0, "DBSCAN")
        self.eCluster.grid(row=7, column=1)


        nb1.pack(fill='both',expand=1,padx=5,pady=5)
        nb1.setnaturalsize()
        self.dialog.activate(geometry = 'centerscreenalways')

    def reset(self, textcur, filecur):
        textcur.delete(0, END)
        textcur.insert(0, filecur)

    def writeNewFile(self):
        path_home = os.environ.get('HOME')
        filename = "%s/WATsite_Settings.txt" % path_home
            
        # write WATsite_Settings.txt file
        fi = open(filename, "w")
        fi.write("USER:\n")
        fi.write("username                           %s\n" % username)
        fi.write("\nCLIENT:\n")
        fi.write("WATsite_home                       %s\n" % WATsite_home)
        fi.write("amber_home                         %s\n" % amber_home)
        fi.write("gromacs_home                       %s\n" % gromacs_home)
        fi.write("pymol_exe_dir                      %s\n" % pymol_exe_dir)
        fi.write("reduce_exe_dir                     %s\n" % reduce_exe_dir)
        fi.write("ncpus                              %s\n" % ncpus)
        fi.write("\n")
        fi.close()
        
    def searchFOLDER(self, filecur, textcur):
        self.folder_select = ""
        self.folder_select = askdirectory(title="Select directory", initialdir=filecur, mustexist=1)
        if self.folder_select:
            textcur.delete(0, END)
            textcur.insert(0, self.folder_select)
            filecur = self.folder_select

    def apply(self, result):
        global username
        global WATsite_home
        global amber_home
        global gromacs_home
        global pymol_exe_dir
        global reduce_exe_dir
        global ncpus
       
        global watermodel
        global clustering
        global MDlength
 
        if result == 'Exit':
            self.dialog.deactivate()
            self.dialog.withdraw()
        elif result == 'Save':
            username                 = self.e_username.get()
            WATsite_home             = self.e_WATsite_home.get()
            amber_home               = self.e_amber_home.get()
            gromacs_home             = self.e_gromacs_home.get()
            pymol_exe_dir            = self.e_pymol_exe_dir.get()
            reduce_exe_dir           = self.e_reduce_exe_dir.get()
            ncpus                    = self.e_ncpus.get()
            self.writeNewFile()

            MDlength = "%s" % (self.eMD.get())
            watermodel = "%s" % (self.eWater.get())
            clustering = "%s" % (self.eCluster.get())
            print watermodel, clustering, MDlength

            mdp_dir = WATsite_home + "/mdp_files/posre/"
            Ri = open("%s/md.mdp" % mdp_dir, "r")
            while Ri:
                line = Ri.readline()
                if line.find("dt") > -1:
                    dt = float( line.split()[-1] )
                elif line.find('nsteps') > -1:
                    nsteps = float( line.split()[-1] )
                elif line == '':
                    break
            Ri.close()
            nsteps = ( float(MDlength) *1000) / dt

            mdp_dir = WATsite_home + "/mdp_files/posre/"
            Fo = open("%s/md_new.mdp" % mdp_dir, "w")
            Ri = open("%s/md.mdp" % mdp_dir, "r")
            for line in Ri:
                if line.find('nsteps') > -1:
                    Fo.write("nsteps            =  %f\n" % nsteps)
                else:
                    Fo.write(line)
            Fo.close()
            Ri.close()

            dat_dir = WATsite_home + "/dat/"
            Ri = open("%s/HydroEntropy.bcf" % dat_dir, "r")
            Fo = open("%s/HydroEntropy_mod.bcf" % dat_dir, "w")
            for line in Ri:
                if line.find('number of frames') > -1:
                    Fo.write("%.1f                                              |number of frames in the snapshot files\n" % (float(MDlength)*1000) )
                elif line.find('water model') > -1:
                    if watermodel[0:3] == 'spc':
                        Fo.write("0                                                |water model (1: TIP4P; 0: SPC)\n")
                    elif watermodel[0:3] == 'tip':
                        Fo.write("1                                                |water model (1: TIP4P; 0: SPC)\n")
                    else:
                        print "SOMETHING WRONG!\n Choose a water models from the four options: spc, spce, tip4p, tip4pew\n\n"
                        sys.exit(0)
                elif line.find('mothod of clustering') > -1:
                    if clustering == 'QT':
                        Fo.write("2                                                |mothod of clustering in hydration site identification (=1: DBSCAN; =2: QT clustering)\n")
                    elif clustering == 'DBSCAN':
                        Fo.write("1                                                |mothod of clustering in hydration site identification (=1: DBSCAN; =2: QT clustering)\n")
                    else:
                        print "SOMETHING WRONG!\n Choose a clustering method from the two options: QT, DBSCAN\n\n"
                        sys.exit(0)
                else:
                    Fo.write(line)
            Fo.close()
            Ri.close()
            #self.dialog.withdraw()
            #self.destroy()

        else:
            username                 = os.getlogin()
            WATsite_home             = "/home/WATsite2.0/"
            amber_home               = "/usr/local/amber14"
            gromacs_home             = "/usr/local/gromacs"
            pymol_exe_dir            = "/usr/local/pymol/"
            reduce_exe_dir           = "/usr/local/reduce"
            ncpus                    = '8'
            self.reset(self.e_username, username)
            self.reset(self.e_WATsite_home, WATsite_home)
            self.reset(self.e_amber_home, amber_home)
            self.reset(self.e_gromacs_home, gromacs_home)
            self.reset(self.e_pymol_exe_dir, pymol_exe_dir)
            self.reset(self.e_reduce_exe_dir, reduce_exe_dir)
            self.reset(self.e_ncpus, ncpus)
            #self.writeNewFile()


#######################################################################################
def modify_settings(app):
    mfs = modifySettingsFile(app.root)

#######################################################################################
class prepareSystem4HydrationSite(tkSimpleDialog.Dialog):
    def body(self, master):
        self.master.title('PyMOL WATsite Plugin')
        self.applic = master
        self.ok_flag = 0

        Label(master, text="Choose target protein and define binding site for WATsite", font=("Helvetica", "14")).grid(row=0, columnspan=5, sticky=N)

        #Label(master, text="  ").grid(row=1, columnspan=5, sticky=N)
        Label(master, text="Project name", relief="sunken", width=80, pady=4).grid(row=1, columnspan=5, sticky=N)
        Label(master, text="Base directory ").grid(row=2, column=0, sticky=W)
        Label(master, text="Project subdirectory ").grid(row=3, column=0, sticky=W)

        curdir = os.getcwd()
        self.eproj = Entry(master, width=80)
        self.eproj.insert(0, curdir)
        self.eproj.grid(row=2, column=1)
        self.bproj = Button(master, text='Browse', command=self.read_base_dir)
        self.bproj.grid(row=2, column=2)

        self.eproj2 = Entry(master, width=80)
        self.eproj2.insert(0, "project_name")
        self.eproj2.grid(row=3, column=1)

        
        Label(master, text="  ").grid(row=9, columnspan=5, sticky=N)
        Label(master, text="Protein selection", relief="sunken", width=120, pady=4).grid(row=10, columnspan=5, sticky=N)
        
        self.v0 = StringVar()
        self.rb_top1 = Radiobutton(master, text="Protein from file [PDB,MOL2,MOL]", variable=self.v0, value="file_single")
        self.rb_top1.grid(row=11, column=0, sticky=W)
        curdir = os.getcwd()
        self.e1_pdb = Entry(master, width=80)
        self.e1_pdb.insert(0, "Full path to protein file")
        self.e1_pdb.grid(row=11, column=1)
        self.e3_pdb = Button(master, text='Search and Import', command=self.searchPDB)
        self.e3_pdb.grid(row=11, column=2)
        
        self.rb_top3 = Radiobutton(master, text="Protein in current PyMol session", variable=self.v0, value="visual")
        self.rb_top3.grid(row=12, column=0, sticky=W)
        
        self.rb_sub0 = [[] for i in range(10000)]
        # top level radio box to select between ligand only, ligand + protein(fixed), ligand + protein(zone), all
        self.v_sub = StringVar()
        ci = 0
        for na in cmd.get_names("objects"):
            cmd.select('pro', na)
            cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO')
            L = cmd.count_atoms('pro2')
            if L > 1:
                self.rb_sub0[ci] = Radiobutton(master, text=na, variable=self.v_sub, value=na)
                self.rb_sub0[ci].grid(row=ci+12, column=1, sticky=W)
                ci += 1
            cmd.delete("pro2")
            cmd.delete("pro")
            
        self.rb_top3.select()
        if ci > 0:
            self.rb_sub0[0].select()
        
        Label(master, text="  ").grid(row=19, columnspan=5, sticky=N)
        Label(master, text="Define binding site", relief="sunken", width=120, pady=4).grid(row=20, columnspan=5, sticky=N)
        
        Label(master, text="Define binding site by known ligand [PDB]").grid(row=21, column=0, sticky=W)
        
        curdir = os.getcwd()
        self.bs_mol = Entry(master, width=80)
        self.bs_mol.insert(0, "Full path to ligand pdb file")
        self.bs_mol.grid(row=21, column=1, sticky=W)
        self.lig_mol = Button(master, text='Search and Import', command=self.searchBindingSiteFile)
        self.lig_mol.grid(row=21, column=2)
        angstrom = u"\u212B"
        Label(master, text="margin").grid(row=21, column=3)
        self.mb = Entry(master, width=5)
        self.mb.insert(0, "3.0")
        self.mb.grid(row=21, column=4)
                
        '''
        self.e1 = Entry(master, width=80)
        self.e1.insert(0, "Full path to binding site file")
        self.e1.grid(row=21, column=0, columnspan=2)
        self.e3 = Button(master, text='Search and Import', command=self.searchLibrary)
        self.e3.grid(row=21, column=2)
        '''

        Label(master, text="  ").grid(row=30, columnspan=5, sticky=N)
        Label(master, text="Hydration sites prediction with or without ligand", relief="sunken", width=120, pady=4).grid(row=31, columnspan=5, sticky=N)
        self.v1 = IntVar()
        self.wl_top1 = Radiobutton(master, text="With ligand from file [MOL2]", variable=self.v1, value=1)
        self.wl_top1.grid(row=32, column=0, sticky=W)
        curdir = os.getcwd()
        self.wl_pdb = Entry(master, width=80)
        self.wl_pdb.insert(0, "Full path to ligand file")
        self.wl_pdb.grid(row=32, column=1)
        self.wl3_pdb = Button(master, text='Search and Import', command=self.searchLIGAND)
        self.wl3_pdb.grid(row=32, column=2)
	Label(master, text="ligand charge:").grid(row=32, column=3)
	self.lig_charge = Entry(master, width=5)
	self.lig_charge.insert(0, "0")
	self.lig_charge.grid(row=32, column=4)


        self.wl_top3 = Radiobutton(master, text="Hydration sites prediction without ligand", variable=self.v1, value=0)
        self.wl_top3.grid(row=33, column=0, sticky=W)               

        Label(master, text="  ").grid(row=34, columnspan=3, sticky=N)
        self.v3 = StringVar()
        #self.rs_top1 = Radiobutton(master, text="Start job in the background", variable=self.v3, value="start")
        #self.rs_top1.grid(row=35, column=0, sticky=W)
        self.rs_top2 = Radiobutton(master, text="Just prepare files for background job", variable=self.v3, value="prep")
        self.rs_top2.grid(row=36, column=0, sticky=W)
        #self.rs_top1.select()
        self.v3.set("prep")
        '''
        Label(master, text="  ").grid(row=37, columnspan=3, sticky=N)
        Label(master, text="Resource", relief="sunken", width=120, pady=4).grid(row=38, columnspan=7, sticky=N)
        self.rc_top = [[] for i in range(100000)]
        self.ec_top = [[] for i in range(100000)]
        self.vc = StringVar()
        self.rc_top[0] = Radiobutton(master, text="local", variable=self.vc, value="local", indicatoron="false", width=50)
        self.rc_top[0].grid(row=39, column=0, sticky=E)
        #Label(master, text="Number of processors/cores [max. %d]: " % (int(os.sysconf("SC_NPROCESSORS_ONLN")))).grid(row=39, column=1, sticky=W)
        Label(master, text="Number of processors/cores [max. %d]: " % (1)).grid(row=39, column=1, sticky=W)
        self.ec_top[0] = Entry(master, width=5)
        self.ec_top[0].insert(0, "1")
        self.ec_top[0].grid(row=39, column=2)
        ci = 1
        for i in SERVER_nodename:
                na = str(i)
                self.rc_top[ci] = Radiobutton(master, text=na, variable=self.vc, value=na, indicatoron="false", width=50)
                self.rc_top[ci].grid(row=ci+39, column=0, sticky=E)
                Label(master, text="Number of processors/cores [max. %d]: " % (SERVER_max_proc[ci-1])).grid(row=ci+39, column=1, sticky=W)
                self.ec_top[ci] = Entry(master, width=5)
                self.ec_top[ci].insert(0, "1")
                self.ec_top[ci].grid(row=ci+39, column=2)
                    
                ci += 1
        self.vc.set("local")
        '''

    def read_base_dir(self):
        base_dir = ""
        base_dir = askdirectory(title="Select base directory", mustexist=1)
        self.eproj.delete(0, END)
        self.eproj.insert(0, base_dir)   
    def searchLibrary(self):
        ftypes=(('pdb file', '*.pdb'), ('mol2 file', '*.mol2'), ("All files", "*"))
        indir = os.getcwd()
        self.bindingsite_name = askopenfilename(initialdir=indir, filetypes=ftypes)
        if self.bindingsite_name:
            self.e1.delete(0, END)
            self.e1.insert(0, self.bindingsite_name)
            cmd.load(self.bindingsite_name)        
    def searchPDB(self):
        ftypes=(('pdb file', '*.pdb'), ('mol2 file', '*.mol2'), ('mol file', '*.mol'), ("All files", "*"))
        indir = os.getcwd()
        self.g_pdb_name = askopenfilename(initialdir=indir, filetypes=ftypes)
        if self.g_pdb_name:
            self.e1_pdb.delete(0, END)
            self.e1_pdb.insert(0, self.g_pdb_name)
            cmd.load(self.g_pdb_name)
    def searchBindingSiteFile(self):
        ftypes=(('pdb file', '*.pdb'), ("All files", "*"))
        indir = os.getcwd()
        self.lig_mol_name = askopenfilename(initialdir=indir, filetypes=ftypes)
        if self.lig_mol_name:
            self.bs_mol.delete(0, END)
            self.bs_mol.insert(0, self.lig_mol_name)
            cmd.load(self.lig_mol_name)

    def searchLIGAND(self):
        ftypes=(('mol2 file', '*.mol2'), ("All files", "*"))
        indir = os.getcwd()
        self.ligand_name = askopenfilename(initialdir=indir, filetypes=ftypes)
        if self.ligand_name:
            self.wl_pdb.delete(0, END)
            self.wl_pdb.insert(0, self.ligand_name)
            cmd.load(self.ligand_name)
            print self.ligand_name
        
    def searchFOLDER(self):
        self.folder_PDB = ""
        self.folder_PDB = askdirectory(title="Select directory", mustexist=1)
        if self.folder_PDB:
            self.e1_folder.delete(0, END)
            self.e1_folder.insert(0, self.folder_PDB)
    
    def validate(self):
        try:
            return 1
        except ValueError:
            tkMessageBox.showwarning(
                "Bad input",
                "Illegal values, please try again"
            )
            return 0

    def apply(self):
        self.project_dir = ""
        self.project_dir_short = ""
        self.project_dir = "%s/%s" % (self.eproj.get(), self.eproj2.get())
        self.project_dir_short = "%s" % (self.eproj2.get())
        self.object_name = self.v0.get()
        self.withligand = self.v1.get()
        self.run_settings = self.v3.get()
        #self.binding_site_method = self.b0.get()
        self.margin = float(self.mb.get())
	global lig_charge
	lig_charge = int(self.lig_charge.get())
        ci = 0
        self.local_or_server = 0
        '''
        j = str(self.vc.get())
        self.computer_name = j
        if j == "local":
                self.local_or_server = 0
                self.numnodes = int(self.ec_top[0].get())
                self.ssh_port = 2200
                self.home_dir = "."
        else:
                self.local_or_server = 1
                for i in SERVER_nodename:
                        if i.find(j) >= 0:
                                self.numnodes = int(self.ec_top[ci+1].get())
                                self.ssh_port = int(SERVER_ssh_port[ci])
                                self.home_dir = SERVER_home_dir[ci]
                ci += 1
        '''
        self.ok_flag = 1

    def cancel(self, event=None):
        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()


#######################################################################################
def write_mol2(sel1, filename):
# find name and object number of selected molecule
        cmd.select('sele_t', sel1)
        cmd.select('sele2', 'sele_t and index 1')
        stored.resnam = []
        cmd.iterate ('sele2', "stored.resnam.append(resn)")
        ci = 0
        for na in cmd.get_names("objects"):
                for sa in stored.resnam:
                        if na == sa:
                            name = ci
                ci += 1
        cmd.delete("sele_t")
        cmd.delete("sele2")

# save pdb and mol file
        cmd.save("%s.pdb" % (filename), sel1, 0, "pdb")
        cmd.save("%s.mol" % (filename), sel1, 0, "mol")
# convert file using babel
        os.system("babel -imol %s.mol -omol2 %s.mol2\n" % (filename, filename))
# change atom names in mol2 file to those in pdb file
        f_pdb = open("%s.pdb" % (filename), "r")
        atnam = []
        resnam = []
        i = 0
        for line in f_pdb:
                if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
                    # Br atoms are written one column to the left by pymol
                    if line[12:13] == "B":
                        atnam.append(line[12:15])
                    elif line[12:13] == "C":
                        atnam.append(line[12:15])
                    else:
                        atnam.append(line[13:16])
                    resnam.append(line[17:21])
                    i = i + 1
        f_pdb.close()
        

        f_old = open("%s.mol2" % (filename), "r")
        f_new = open("new.mol2", "w")

        
        flag = 0
        for line in f_old:
                if line[0:13] == '@<TRIPOS>BOND':
                        flag = 0
                if flag == 0:
                        f_new.write(line)
                else:
                        new_line = []
                        new_line.append(line[0:8])
                        new_line.append(str(atnam.pop(0)))
                        new_line.append(line[11:58])
                        new_line.append(str(resnam.pop(0)))
                        new_line.append(line[62:])
                        f_new.write(''.join(new_line))
                if line[0:13] == '@<TRIPOS>ATOM':
                        flag = 1

        f_old.close()
        f_new.close()

        os.remove("%s.mol2" % (filename))
        os.rename("new.mol2", "%s.mol2" % (filename))
        os.remove("%s.mol" % (filename))
        os.remove("%s.pdb" % (filename))


#######################################################################################
def prepare_protein(app):
        global g_protein_name_slide
        global object_name
                
        curdir = os.getcwd()
        print curdir
        
        ci = 0
        for na in cmd.get_names("objects"):
            cmd.select('pro', na)
            cmd.select('pro2', 'pro and resn ALA+ARG+ASH+ASN+ASP+CYS+CYX+CY1+GLH+GLN+GLU+GLY+HIS+HIE+HID+HIP+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL+HEM+HEO')
            L = cmd.count_atoms('pro2')
            if L > 1:
                object_name = na
                ci += 1
            cmd.delete("pro2")
            cmd.delete("pro")
        
        nprep = prepareSystem4HydrationSite(app.root)
        if nprep.ok_flag == 0:
            return 0
    
        md_outfile = nprep.project_dir
        md_save_project = 1
        if os.path.isdir(md_outfile):
            os.system("rm %s -r"%md_outfile)
        os.mkdir(md_outfile)
        curdir = os.getcwd()
        # change into project directory
        os.chdir(md_outfile)   
    
        if nprep.object_name == "visual":
            object_name = nprep.v_sub.get()
            for na in cmd.get_names("objects"):
                if na == object_name:
                    # remove hydrogens from protein
                    cmd.select('pro', na)
                    cmd.save("MyProtein.pdb", 'pro', 0, "pdb")
                    cmd.delete("pro")
        elif nprep.object_name == "file_single":
            os.system("cp %s ./MyProtein.pdb" % nprep.g_pdb_name)


        '''
        if nprep.binding_site_method == "box":
                #cmd.load(nprep.g_pdb_name)
                nsize = sizeSearchSpaceDialog(app.root)
                if nsize.ok_flag == 0:
                        return 0
                pharmdock_min_x = nsize.cb_x - nsize.bl_x/2
                pharmdock_min_y = nsize.cb_y - nsize.bl_y/2
                pharmdock_min_z = nsize.cb_z - nsize.bl_z/2
                pharmdock_max_x = nsize.cb_x + nsize.bl_x/2
                pharmdock_max_y = nsize.cb_y + nsize.bl_y/2
                pharmdock_max_z = nsize.cb_z + nsize.bl_z/2
                fconf = open("%s/MyBindingSite.pdb"%md_outfile, "w")
                fconf.write("COMPND    binding_site\n")
                fconf.write("AUTHOR    WATsite\n")
                fconf.write("ATOM      1  C01 BDS A   1     %8.3f%8.3f%8.3f   1.00  0.00           C\n"%(pharmdock_min_x, pharmdock_min_y, pharmdock_min_z))  
                fconf.write("ATOM      2  C02 BDS A   1     %8.3f%8.3f%8.3f   1.00  0.00           C\n"%(pharmdock_max_x, pharmdock_max_y, pharmdock_max_z)) 
                fconf.close()
                bs_margin = 0
        '''
        bs_margin = nprep.margin
        os.system("cp %s %s/MyBindingSite.pdb"%(nprep.lig_mol_name, md_outfile))


        runscript = open("RunScript", "w")
        runscript.write("#!/bin/tcsh\n")
        runscript.write("setenv WATSITEHOME %s\n"%WATsite_home) 
        runscript.write("setenv AMBERHOME %s\n"%amber_home)
        runscript.write("setenv GROMACS_HOME %s\n"%gromacs_home)
        runscript.write("setenv PYMOL_EXE %s/pymol\n"%pymol_exe_dir)
        runscript.write("setenv reduce_exe_dir %s\n"%reduce_exe_dir)
        runscript.write("setenv ncpus %s\n"%ncpus)
        runscript.write("setenv watermodel %s\n"%watermodel)
        runscript.write("setenv withligand %d\n"%nprep.withligand)
	runscript.write("setenv ligcharge %d\n"% lig_charge )
        runscript.write("$WATSITEHOME/bin/WATsite.sh >& watsite.log\n")
        runscript.close()
        #os.system("cp %s ./binding_site.pdb" % (nprep.bindingsite_name))##!!!this need to be changed in HydroEntropy program
        if nprep.withligand:
                os.system("cp %s ./MyLigand.mol2" % (nprep.ligand_name))

            
        if nprep.run_settings == "prep":
            if nprep.local_or_server == 0:
                childp = subprocess.Popen("cd ./%s; chmod u+rwx RunScript" % (nprep.project_dir), shell=True)
                tkMessageBox.showinfo("Preparation of hydration site analysis done","Please enter folder %s and start job by typing 'nohup ./RunScript &'.\nData can be read in by using 'WATsite, Import Results'." % (nprep.project_dir))
            '''
            else:
                childp = subprocess.Popen("ssh -p %d %s@%s 'cd %s/%s; chmod u+rwx RunScript'" % (nprep.ssh_port, username, nprep.computer_name, nprep.home_dir, nprep.project_dir), shell=True)
                tkMessageBox.showinfo("Preparation of WATsite simulation done","Please copy folder %s_client/%s to server and start job by typing 'nohup ./RunScript &'.\nData can be read in by using 'MM, Monitor molecular mechanics results'." % (nprep.project_dir, nprep.project_dir))
            '''
        #else:
        #    if nprep.local_or_server == 0:
        #        childp = subprocess.Popen("cd ./%s; chmod u+rwx RunScript; run RunScript &" % (nprep.project_dir_short), shell=True)
        #        tkMessageBox.showinfo("Preparation of hydration site analysis done","Job has been started on localhost.\nData can be read in by using 'WATsite, Import Results'.")
            '''
            else:
                os.system("ssh -p %d %s@%s 'cd %s; mkdir %s'" % (nprep.ssh_port, username,nprep.computer_name, nprep.home_dir, nprep.project_dir))
                os.system("scp -P %d -r %s/* RunScript %s@%s:%s/%s/" % (nprep.ssh_port, md_outfile, username, nprep.computer_name, nprep.home_dir, nprep.project_dir))
                childp = subprocess.Popen("ssh -p %d %s@%s 'cd %s/%s; chmod u+rwx RunScript; /usr/pbs/bin/qsub -q mlill RunScript'" % (nprep.ssh_port, username, nprep.computer_name,nprep.home_dir, nprep.project_dir), shell=True)
                tkMessageBox.showinfo("Preparation of WATsite simulation done","Job has been started on server.\nData can be read in by using 'MM, Monitor molecular mechanics results'.")
            '''        
        # write monitor file
        fmon = open("%s/WATsite.out" % (md_outfile), "w")
        fmon.write("Protein             MyProtein.pdb\n")
        fmon.write("Ligand              MyLigand.mol2\n")
        fmon.write("HydrationSiteMol2   HydrationSites.mol2\n")
        fmon.write("HydrationSitePDB    HydrationSites.pdb\n")
        fmon.write("WaterTraj           OUTPUT/WATinside.mol2\n")
        fmon.write("EnergyFile          cluster.egy\n")
        fmon.close()
        
        os.chdir(curdir)




#######################################################################################
def import_results(app):

    global num_solutions
    global num_frames
    global import_dir
    global HydrationSiteEnergyFile
    global HydrationSiteMol2
    global num_WATsite
    global TdS
    global dH
    global dG
    global Occupancy
    
    mdp_dir = WATsite_home + "/mdp_files/posre/"
    try:
        Ri = open("%s/md_new.mdp" % mdp_dir, "r")
    except:
        Ri = open("%s/md.mdp" % mdp_dir, "r")
    while Ri:
        line = Ri.readline()
        if line.find("dt") > -1:
            dt = float( line.split()[-1] )
        elif line.find('nsteps') > -1:
            nsteps = float( line.split()[-1] )
        elif line == '':
            break
    #print dt, nsteps
    nframe = dt * nsteps 
    numFrame = float(nframe)
    print "in HydrationSiteEnergy: total number of frame %f\n" % float(numFrame)
    ftypes=(('out file', '*.out'), ("All files", "*"))
    curdir = os.getcwd()
    openfile = askopenfilename(initialdir=curdir, filetypes=ftypes)
    if openfile:
        import_dir = os.path.dirname(openfile)
        os.chdir(import_dir)
        
        fin = open(openfile, "r")
        HSwithLig = 0
        for line in fin:
                if line.find("Protein") > -1:
                        ProtFile = line.split()[1]
                if line.find("Ligand") > -1:
                        LigFile = line.split()[1]
                        HSwithLig = 1
                if line.find("HydrationSiteMol2") > -1:
                        HydrationSiteMol2 = line.split()[1]
                if line.find("HydrationSitePDB") > -1:
                        HydrationSitePDB = line.split()[1]
                if line.find("WaterTraj") > -1:
                        WaterTraj = line.split()[1]
                if line.find("EnergyFile") > -1:
                        HydrationSiteEnergyFile = line.split()[1]
        fin.close()
        
        flag = 0
        for ent in os.listdir(import_dir):
            print ent
            if ent.find(HydrationSiteEnergyFile) >= 0:
                flag = 1
        if flag == 0:
            tkMessageBox.showwarning(
                "Could not find Hydratioin Site Energy files",
                "WATsite might be still in progress or didn't finish successfully"
            )
            return 0

        TdS, dH, dG, Occupancy = [], [], [], []
        sie_file = open(HydrationSiteEnergyFile, 'r')
        for i in sie_file:
                if i.find("@") == -1 and i != "":
                        j = i.split()
                        TdS.append(float(j[1]))
                        dH.append(float(j[2])+18.18)
                        dG.append(float(j[3]))
                        Occupancy.append(float(j[8])/numFrame )
        sie_file.close()
        num_WATsite = len(dG)
        z = loadDialog(app.root)

        if z.ok_flag == 0:
                return
        if tkMessageBox.askokcancel("Read project", "Reading results will delete all exisiting objects in current session.") == 0:
                return 0
        # remove existing objects
        for na in cmd.get_names("objects"):
            try:
                cmd.remove(na)
            except:
                pass
            cmd.delete(na)

        if z.protfile.get() == 1:
                try:
                        fi = open(ProtFile, 'r')
                except:
                        tkMessageBox.showwarning(
                        "Missing protein file",
                        "Please, check the protein file specified in WATsite.out is existing in the directory"
                        )
                        return 0
                cmd.load(ProtFile)
        if z.ligfile.get() == 1 and HSwithLig:
                try:
                        fi = open(LigFile, 'r')
                except:
                        tkMessageBox.showwarning(
                        "Missing protein file",
                        "Please, check the ligand file specified in WATsite.out is existing in the directory"
                        )
                        return 0
                cmd.load(LigFile)
        if z.HSmol2file.get() == 1:
                try:
                        fi = open(HydrationSiteMol2, 'r')
                except:
                        tkMessageBox.showwarning(
                        "Missing protein file",
                        "Please, check the hydration site mol2 file specified in WATsite.out is existing in the directory"
                        )
                        return 0

                cmd.load(HydrationSiteMol2)
                cmd.select("tmp", HydrationSiteMol2.split(".")[0])
                cmd.show_as("nb_spheres","tmp")
                #cmd.hide("nonbonded", "tmp")
                cmd.spectrum("pc", "blue_white_red", "tmp")
                #cmd.label("tmp", )
                cmd.delete("tmp")
        if z.WATtraj.get() == 1:
                cmd.load(WaterTraj)
                
        HydrationSiteEnergy(app.root)

            
#######################################################################################
class loadDialog(tkSimpleDialog.Dialog):
        def body(self, master):
                self.master.title('PyMOL WATsite Plugin')
                self.var_check = []

                self.ok_flag = 0

                Label(master, text="Please select files to load").grid(row=0, column=0, sticky=W)
                self.protfile = IntVar()
                self.a = Checkbutton(master, text="Protein", variable=self.protfile)
                self.a.grid(row=1, column=0, sticky=W, padx=5)
                self.protfile.set(1)
                
                self.ligfile = IntVar()
                self.b = Checkbutton(master, text="Ligand", variable=self.ligfile)
                self.b.grid(row=2, column=0, sticky=W, padx=5)
                self.ligfile.set(0)
                
                self.HSmol2file = IntVar()
                self.c = Checkbutton(master, text="Hydration Site", variable=self.HSmol2file)
                self.c.grid(row=3, column=0, sticky=W, padx=5)
                self.HSmol2file.set(1)
                
                self.WATtraj = IntVar()
                self.d = Checkbutton(master, text="Binding site water trajectory", variable=self.WATtraj)
                #self.d.grid(row=4, column=0, sticky=W, padx=5)
                self.WATtraj.set(0)
        def validate(self):
                try:
                        return 1

                except ValueError:
                        tkMessageBox.showwarning(
                                "Bad input",
                                "Illegal values, please try again"
                        )
                        return 0

        def apply(self):
                self.ok_flag = 1

        def cancel(self, event=None):
                # put focus back to the parent window
                self.parent.focus_set()
                self.destroy()

#######################################################################################
class HydrationSiteEnergy:
        def __init__(self, top):
                        self.dialog = Pmw.Dialog(top, buttons = ('OK','Cancel'), defaultbutton = 'OK', title = 'WATsite results', command = self.apply)
                        parent = self.dialog.interior()

                        group1 = Pmw.Group(parent, tag_text='Energy values in kcal/mol (double click to select hydration site in Pymol)')
                        master = Frame(group1.interior())
                        scrollbar = Scrollbar(master)
                        listbox = Listbox(master, yscrollcommand=scrollbar.set)
                        listbox.config(font=('Courier', 15))
                        #delta = u"\u0394"
                        delta = 'd'
                        TXT='HS#      -T%sS      %sH        %sG    occupancy'%(delta, delta, delta)
                        listbox.insert(END, TXT)
                        sie_file = open(HydrationSiteEnergyFile, 'r')
                        ci = 3
                        hsi = 1
                        pos = 0
                        mdp_dir = WATsite_home + "/mdp_files/posre/"
                        try:
                            Ri = open("%s/md_new.mdp" % mdp_dir, "r")
                        except:
                            Ri = open("%s/md.mdp" % mdp_dir, "r")
                        while Ri:
                            line = Ri.readline()
                            if line.find("dt") > -1:
                                dt = float( line.split()[-1] )
                            elif line.find('nsteps') > -1:
                                nsteps = float( line.split()[-1] )
                            elif line == '':
                                break
                        nframe = dt * nsteps 
                        numFrame = float(nframe)
                        print "in HydrationSiteEnergy: total number of frame %f\n" % float(numFrame)

                        self.S, self.H, self.G, self.O = [], [], [], []
                        for i in sie_file:
                                        if i.find("@") == -1 and i != "":
                                                        j = i.split()
                                                        self.S.append(float(j[1]))
                                                        self.H.append(float(j[2])+18.18)
                                                        self.G.append(float(j[3]))
                                                        self.O.append(float(j[8])/numFrame)
                                                #Label(master, text='%3d    %6.2f    %6.2f    %6.2f    %5.2f' % (hsi, float(j[1]), float(j[2])+18.18, float(j[3]), float(j[8])/numFrame), font=('Courier', 10, 'bold')).grid(row=ci, column=0, sticky=W)
                                                        TXT = '%3d    %6.2f    %6.2f    %6.2f    %6.2f' % (hsi, float(j[1]), float(j[2])+18.18, float(j[3]), float(j[8])/numFrame)
                                                        listbox.insert(END, TXT)
                                                #pos += 1
                                                        ci += 1
                                                        hsi += 1
                        sie_file.close()
                        
                        scrollbar.config(command=listbox.yview)
                        listbox.bind('<Double-1>', self.handleList)
                        self.listbox = listbox
                        scrollbar.pack(side=RIGHT, fill=Y)
                        listbox.pack(side=LEFT, expand=YES, fill=BOTH)
                        master.pack(fill = 'both', expand = 1, padx = 5, pady=1)
                        group1.pack(fill = 'both', expand = 1, padx = 5, pady = 5)
                        
                        
                        group2 = Pmw.Group(parent, tag_text='Select how to color the hydration sites')
                        f2 = Frame(group2.interior())
                        
                        b1 = Button(f2, text='-T%sS'%delta, command=self.colorByS)
                        b1.grid(row=2, column=2)
                        b2 = Button(f2, text='%sH'%delta, command=self.colorByH)
                        b2.grid(row=2, column=4)
                        b3 = Button(f2, text='%sG'%delta, command=self.colorByG)
                        b3.grid(row=2, column=5)
                        b4 = Button(f2, text='Occupancy', command=self.colorByO)
                        b4.grid(row=2, column=6)
                        f2.pack(fill = 'both', expand = 0, padx = 5, pady=1)
                        group2.pack(fill = 'both', expand = 0, padx = 5, pady = 5)
                        
                        self.dialog.activate(geometry = 'centerscreenalways')

        def colorByS(self):
                        cmd.select("colorByS", HydrationSiteMol2.split(".")[0])
                        for i in range(len(self.S)):
                                        cmd.alter("colorByS and id %s"%(i+1), "b=%s"%self.S[i])
                        cmd.spectrum("b", "blue_white_red", "colorByS",minimum=min(self.S), maximum=max(self.S))
                        cmd.delete("colorByS")
        def colorByH(self):
                        cmd.select("colorByH", HydrationSiteMol2.split(".")[0])
                        for i in range(len(self.H)):
                                        cmd.alter("colorByH and id %s"%(i+1), "b=%s"%self.H[i])
                        cmd.spectrum("b", "blue_white_red", "colorByH",minimum=min(self.H), maximum=max(self.H))
                        cmd.delete("colorByH")
        def colorByG(self):
                        cmd.select("colorByG", HydrationSiteMol2.split(".")[0])
                        for i in range(len(self.G)):
                                        cmd.alter("colorByG and id %s"%(i+1), "b=%s"%self.G[i])
                        cmd.spectrum("b", "blue_white_red", "colorByG",minimum=min(self.G), maximum=max(self.G))
                        cmd.delete("colorByG")
        def colorByO(self):
                        cmd.select("colorByO", HydrationSiteMol2.split(".")[0])
                        for i in range(len(self.O)):
                                        cmd.alter("colorByO and ID %s"%(i+1), "b=%s"%self.O[i])
                        cmd.spectrum("b", "blue_white_red", "colorByO",minimum=min(self.O), maximum=max(self.O))
                        cmd.delete("colorByO")
                        
        def handleList(self, event):
                        index = self.listbox.curselection()
                        label = self.listbox.get(index)
                        self.runCommand(label)

        def runCommand(self, selection):
                        hs_index = selection.split()[0]
                        cmd.select("selected_hs", HydrationSiteMol2.split(".")[0] + ' and id '+hs_index)
                        #print HydrationSiteMol2.split(".")[0], hs_index
        def apply(self, result):
                        if result == 'OK':
                                        self.ok_flag = 1
                                        self.dialog.deactivate()
                                        self.dialog.withdraw()
                        else:
                                        self.ok_flag = 0
                                        self.dialog.deactivate()
                                        self.dialog.withdraw()
#######################################################################################
def estimate_desolvation(app):
        
        global liglist
        global num_ligands
        global desol_G, desol_S, desol_H #desolvation energies
        desol_G, desol_S, desol_H = [], [], []
        
        imp_lig = import_ligands(app.root)
        
        sel_dir = imp_lig.project_dir
        liglist = []
        for lf in os.listdir(sel_dir):
                if lf.find("mol2") > -1 or lf.find("pdb") > -1:
                        cmd.load("%s/%s"%(sel_dir,lf))
                        ligname = lf.split(".")[0]
                        cmd.remove("%s and hydrogen"%ligname)
                        liglist.append(ligname)
        num_ligands = len(liglist)
        
        try:
                HydrationSiteMol2
        except NameError:
                tkMessageBox.showwarning(
                "Bad input",
                "Please import hydration site results first"
                )
                return 0

        hs = HydrationSiteMol2.split(".")[0]
        #set b factor values by TdS, dH, dG
        for i in range(num_WATsite):
                cmd.alter('%s and id %s'%(hs,i+1), "b=%s"%dG[i])
                cmd.alter('%s and id %s'%(hs,i+1), "partial_charge=%s"%TdS[i])
                cmd.alter('%s and id %s'%(hs,i+1), 'vdw=%s'%dH[i])
        ligG, ligS, ligH = [], [], []
        for lig in liglist:
                sel_hs = '%s_replacedHS'%lig
                cmd.select(sel_hs, '%s w. %.2f of %s'%(hs, imp_lig.HSSradius, lig))
                #print '%s w. %.2f of %s'%(hs, imp_lig.HSSradius, lig)
                myspace = {'ligG': [], 'ligS': [], 'ligH': []}
                cmd.iterate(sel_hs, 'ligG.append(b)', space=myspace)
                cmd.iterate(sel_hs, 'ligS.append(partial_charge)', space=myspace)
                cmd.iterate(sel_hs, 'ligH.append(vdw)', space=myspace)
                #print myspace['ligG']
                sumG, sumS, sumH = 0.0, 0.0, 0.0
                for i in myspace['ligG']:
                        sumG += float(i)
                for i in myspace['ligS']:
                        sumS += float(i)
                for i in myspace['ligH']:
                        sumH += float(i)
                #sumG = sum(ligG)
                #sumS = sum(ligS)
                #sumH = sum(ligH)
                #print sumG, sumS, sumH
                desol_G.append(sumG)
                desol_S.append(sumS)
                desol_H.append(sumH)
        displayResultsDialog(app.root)
#######################################################################################        
class import_ligands(tkSimpleDialog.Dialog):
        def body(self, master):
                self.master.title('PyMOL WATsite Plugin')
                self.applic = master
                self.ok_flag = 0

                Label(master, text="Estimate desolvation energy for ligands", font=("Helvetica", "12")).grid(row=0, columnspan=5, sticky=N)

                curdir = os.getcwd()
                Label(master, text="Load ligands", font=("Helvetica", "12")).grid(row=2, column=1, sticky=W)

                self.eproj = Entry(master, width=80)
                self.eproj.insert(0, "Full path to ligand library [mol2 or pdb]")
                self.eproj.grid(row=2, column=2)
                self.bproj = Button(master, text='Browse', command=self.read_base_dir)
                self.bproj.grid(row=2, column=3)
                
                Label(master, text="Hydration Site Selection Radius (Angstrom)", font=("Helvetica", "12")).grid(row=3, column=1, sticky=W)

                self.hssr = Entry(master, width=10)
                self.hssr.insert(0, "1.0")
                self.hssr.grid(row=3, column=2, sticky=W)

        def read_base_dir(self):
                base_dir = ""
                base_dir = askdirectory(title="Select base directory", mustexist=1)
                self.eproj.delete(0, END)
                self.eproj.insert(0, base_dir)
        def apply(self):
                self.project_dir = "%s" % (self.eproj.get())
                self.HSSradius = float(self.hssr.get())
                self.ok_flag = 1

        def cancel(self, event=None):
                # put focus back to the parent window
                self.parent.focus_set()
                self.destroy()
                
#######################################################################################
class displayResultsDialog:
        def __init__(self, parent):
                top = self.top = Toplevel(parent)
        
                self.ok_flag = 0
                
                Label(top, text="Desolvation Energy Estimate", font=("Helvetica", "14", "bold")).grid(row=0, columnspan=4, sticky=N)
                Label(top, text="(Unit: kcal/mol)", font=("Courier", "12")).grid(row=1, columnspan=4, sticky=N)
                Label(top, text="  ").grid(row=2, columnspan=4, sticky=N)
                
                Label(top, text="Ligand name        ", font=("Courier", "12")).grid(row=10, column=0, sticky=W)
                Label(top, text="%sG "%delta, font=("Courier", "12")).grid(row=10, column=1, sticky=E)
                Label(top, text="-T%sS "%delta, font=("Courier", "12")).grid(row=10, column=2, sticky=E)
                Label(top, text="%sH "%delta, font=("Courier", "12")).grid(row=10, column=3, sticky=E)

                ci = 0
                for na in range(num_ligands):
                    Label(top, text=liglist[na], font=("Courier", "12")).grid(row=ci+11, column=0, sticky=W)
                    Label(top, text="%.2f " % desol_G[na], font=("Courier", "12")).grid(row=ci+11, column=1, sticky=E)
                    Label(top, text="%.2f " % desol_S[na], font=("Courier", "12")).grid(row=ci+11, column=2, sticky=E)
                    Label(top, text="%.2f " % desol_H[na], font=("Courier", "12")).grid(row=ci+11, column=3, sticky=E)
                    ci += 1
                b = Button(top, text="Close", padx=5, command=self.ok)        
                b.grid(row=num_ligands+30, column=1, sticky=S)
        def ok(self):
                self.top.destroy()

       

read_settings()

