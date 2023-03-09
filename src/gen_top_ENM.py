#!/usr/bin/env python
# original perl script was made
# BY RUSSELL DEVANE 2/2008
###############################################################################
# READS CG PDB FILE AND DUMPS CONNECTIVITY BASED ON MODEL
# WILL GENERATE ALL BACKBONE-BACKBONE BONDS 
# SIDECHAIN BACKBONE AND SIDECHAIN SIDECAIN BONDS BASED ON MAPPINGS BELOW
# WILL DUMP ANGLES IF WANTED BASED ON ALL BONDS
# WILL DUMP SIDECHAIN IMPROPERS IF WANTED
# WILL DUMP BACKBONE TORSIONS IF WANTED
################################################################################
#
# USAGE: gen_elastic_network.py <cg.pdb filename> < cg.top filename (output)>
#
#################################################################################
# CAN HANDLE MULTIPLE CHAINS
# DOES REVERSE ANGLES NOW
#
#################################################################################
# DOES REQUIRE THE CHEMISTRY MODULE
#################################################################################

#################################################################################
# By Shuhei Kawamoto 3/2012
#   Deleted angle for ring of side chain.
#   Added terminal backbone particle GBT and ABT.
#   Fixed bugs.
#################################################################################
# and converted to PYTHON from PERL
# BY YUSUKE MIYAZAKI 2020/3
#################################################################################
# and converted to PYTHON3 from PYTHON2
# BY ISSEI KAWABATA 2022/4
#################################################################################

import sys, math
import numpy as np
from argparse import ArgumentParser
import subprocess
from setup_lmp import get_angle  

def get_option():
    # maximum bond length 
    MAXdr = 9.0
    # force constant for elastic network model 
    kENM  = 1.195
    argparser = ArgumentParser()
    argparser.add_argument('cgpdb', type=str,
                            help='Specify input CG PDB file name.')
    argparser.add_argument('aapdb', type=str,
                            help='Specify input AA PDB file name.')
    argparser.add_argument('output', type=str,
                            help='Specify output topology file name.')
    argparser.add_argument('-maxr', type=float,
                            default=MAXdr,
                            help='Cutoff length of ENM (default: 9.0 A).')
    argparser.add_argument('-kENM', type=float,
                            default=kENM,
                            help='Force constant for ENM (default: 1.195 kcal/A2).')
    argparser.add_argument('-pspica', action='store_true',
                            help='Assign partial charge (0.5990) for pSPICA FF (default: 0.1118, for SPICA FF).')
    argparser.add_argument('-dssp', type=str,
                            default='dssp',
                            help='Specify path to dssp binary')
    return argparser.parse_args()

def get_option_script(argv):
    # maximum bond length 
    MAXdr = 9.0
    # force constant for elastic network model 
    kENM  = 1.195
    argparser = ArgumentParser(usage='ENM [-h] [-maxr MAXR] [-kENM KENM] [-pspica] [-dssp dssp] cgpdb aapdb output',
                               prog ="ENM")
    argparser.add_argument('cgpdb', type=str,
                            help='Specify input CG PDB file name.')
    argparser.add_argument('aapdb', type=str,
                            help='Specify input AA PDB file name.')
    argparser.add_argument('output', type=str,
                            help='Specify output topology file name.')
    argparser.add_argument('-maxr', type=float,
                            default=MAXdr,
                            help='Cutoff length of ENM (default: 9.0 A).')
    argparser.add_argument('-kENM', type=float,
                            default=kENM,
                            help='Force constant for ENM (default: 1.195 kcal/A2).')
    argparser.add_argument('-pspica', action='store_true',
                            help='Assign partial charge (0.5990) for pSPICA FF (default: 0.1118, for SPICA FF).')
    argparser.add_argument('-dssp', type=str,
                            default='dssp',
                            help='Specify path to dssp binary')
    return argparser.parse_args(argv)

ncomp = 7
ncomp_pdb = 15
(    PDB_RECNAME,
    PDB_INDEX,
    PDB_ATMNAME,
    PDB_INDICAT,
    PDB_RESNAME,
    PDB_CHAINID,
    PDB_RESID,
    PDB_CODE, 
    PDB_POSX, 
    PDB_POSY,
    PDB_POSZ,
    PDB_OCCUP,
    PDB_TFACT,
    PDB_SEGID,
    PDB_ELESYM 
) = range(ncomp_pdb)
(    ATMNAME,
    RESNAME,
    RESID,
    CHAINID,
    POSX, 
    POSY,
    POSZ,
) = range(ncomp)

res_map = {}
this_type = {}

# CG map for DNA/RNA #
res_map["ADE"] = ["AD1", "AD2", "AD3", "AD4"]
res_map["GUA"] = ["GU1", "GU2", "GU3", "GU4"]
res_map["CYT"] = ["CY1", "CY2", "CY3"]
res_map["THY"] = ["TH1", "TH2", "TH3"]
res_map["URA"] = ["UR1", "UR2", "UR3"]

this_type["PB"]  = "PH"
this_type["RB1"] = "RB1"
this_type["RB2"] = "RB2"
this_type["DB2"] = "DB2"
this_type["AD1"] = "BPY2"
this_type["AD2"] = "BPU1"
this_type["AD3"] = "BPU2"
this_type["AD4"] = "BPU2"
this_type["GU1"] = "BPY1"
this_type["GU2"] = "BPU1"
this_type["GU3"] = "BPU2"
this_type["GU4"] = "BPYG"
this_type["CY1"] = "BPY2"
this_type["CY2"] = "BPY1"
this_type["CY3"] = "CMD2"
this_type["TH1"] = "BPY0"
this_type["TH2"] = "BPY1"
this_type["TH3"] = "CMD3"
this_type["UR1"] = "BPY0"
this_type["UR2"] = "BPY1"
this_type["UR3"] = "CMD2"

# CG map for protein #
res_map["ALA"] = ["ALA"]
res_map["VAL"] = ["VAL"]
res_map["LEU"] = ["LEU"]
res_map["ILE"] = ["ILE"]
res_map["PRO"] = ["PRO"]
res_map["MET"] = ["MET"]
res_map["SER"] = ["SER"]
res_map["THR"] = ["THR"]
res_map["CYS"] = ["CYS"]
res_map["ASN"] = ["ASN"]
res_map["GLN"] = ["GLN"]
res_map["ASP"] = ["ASP"]
res_map["GLU"] = ["GLU"]
res_map["ARG"] = ["AR1", "AR2"]
res_map["LYS"] = ["LY1", "LY2"]
res_map["PHE"] = ["PH1", "PH2", "PH3", "PH4"]
res_map["TYR"] = ["TY1", "TY2", "TY3", "TY4"]
res_map["TRP"] = ["TR1", "TR2", "TR3"]
res_map["HIS"] = ["HI1", "HI2", "HI3"]
res_map["HSP"] = ["HI1", "HI2", "HI3"]
res_map["HSD"] = ["HI1", "HI2", "HI3"]
res_map["HSE"] = ["HI1", "HI2", "HI3"]

this_type["GBB"] = "GBB"
this_type["GBM"] = "GBM"
this_type["ABB"] = "ABB"
this_type["VAL"] = "VAL"
this_type["LEU"] = "LEU"
this_type["ILE"] = "ILE"
this_type["PRO"] = "PRO"
this_type["MET"] = "MET"
this_type["PH1"] = "XYR"
this_type["PH2"] = "BER"
this_type["PH3"] = "BER"
this_type["PH4"] = "BER"
this_type["TR1"] = "TR1"
this_type["TR2"] = "TR2"
this_type["TR3"] = "TR3"
this_type["SER"] = "SER"
this_type["THR"] = "THR"
this_type["CYS"] = "CYS"
this_type["ASN"] = "ASN"
this_type["GLN"] = "GLN"
this_type["TY1"] = "XYR"
this_type["TY2"] = "BER"
this_type["TY3"] = "PHR"
this_type["TY4"] = "BER"
this_type["ASP"] = "ASP"
this_type["GLU"] = "GLU"
this_type["LY1"] = "LY1"
this_type["LY2"] = "LY2"
this_type["AR1"] = "AR1"
this_type["AR2"] = "AR2"
this_type["HI1"] = "HI1"
this_type["HI2"] = "HI2"
this_type["HI3"] = "HI3"
this_type["GBT"] = "GBT"
this_type["ABT"] = "ABT"

###############################################################################
#  IDENTIFY UNIQUE SIDECHAIN BONDS...I.E. RINGS
#  OR MULTIPLE SITES..
#  IF NOTHING HERE IT ASSUMES STRAIGHTFORWARD BONDING
######################################################
spec_bond={}

rbbnd1 = ["PB",  "RB1"]
rbbnd2 = ["RB1", "RB2"]
dbbnd2 = ["RB1", "DB2"]

adbnd1 = ["DB2", "AD3"]
adbnd2 = ["AD1", "AD2"]
adbnd3 = ["AD2", "AD3"]
adbnd4 = ["AD3", "AD4"]
adbnd5 = ["AD4", "AD1"]
adbnd6 = ["AD2", "AD4"]
spec_bond["ADE"] = [adbnd1, adbnd2, adbnd3, adbnd4, adbnd5]

gubnd1 = ["DB2", "GU3"]
gubnd2 = ["GU1", "GU2"]
gubnd3 = ["GU2", "GU3"]
gubnd4 = ["GU3", "GU4"]
gubnd5 = ["GU4", "GU1"]
gubnd6 = ["GU2", "GU4"]
spec_bond["GUA"] = [gubnd1, gubnd2, gubnd3, gubnd4, gubnd5]

cybnd1 = ["DB2", "CY2"]
cybnd2 = ["CY1", "CY2"]
cybnd3 = ["CY2", "CY3"]
cybnd4 = ["CY3", "CY1"]
spec_bond["CYT"] = [cybnd1, cybnd2, cybnd3, cybnd4]

thbnd1 = ["DB2", "TH2"]
thbnd2 = ["TH1", "TH2"]
thbnd3 = ["TH2", "TH3"]
thbnd4 = ["TH3", "TH1"]
spec_bond["THY"] = [thbnd1, thbnd2, thbnd3, thbnd4]

urbnd1 = ["RB2", "UR2"]
urbnd2 = ["UR1", "UR2"]
urbnd3 = ["UR2", "UR3"]
urbnd4 = ["UR3", "UR1"]
spec_bond["URA"] = [urbnd1, urbnd2, urbnd3, urbnd4]

lysbnd1 = ["GBM", "LY1"]
lysbnd2 = ["LY1", "LY2"]
spec_bond["LYS"] = [lysbnd1, lysbnd2]

argbnd1 = ["GBM", "AR1"]
argbnd2 = ["AR1", "AR2"]
spec_bond["ARG"] = [argbnd1, argbnd2]

hisbnd1 = ["GBM", "HI1"]
hisbnd2 = ["HI1", "HI2"]
hisbnd3 = ["HI1", "HI3"]
hisbnd4 = ["HI2", "HI3"]
spec_bond["HIS"] = [hisbnd1, hisbnd2, hisbnd3, hisbnd4];
spec_bond["HSP"] = [hisbnd1, hisbnd2, hisbnd3, hisbnd4];
spec_bond["HSD"] = [hisbnd1, hisbnd2, hisbnd3, hisbnd4];
spec_bond["HSE"] = [hisbnd1, hisbnd2, hisbnd3, hisbnd4];

trpbnd1 = ["GBM", "TR1"]
trpbnd2 = ["TR1", "TR2"]
trpbnd3 = ["TR1", "TR3"]
trpbnd4 = ["TR2", "TR3"]
spec_bond["TRP"] = [trpbnd1, trpbnd2, trpbnd3, trpbnd4];

phebnd1 = ["GBM", "PH1"]
phebnd2 = ["PH1", "PH2"]
phebnd3 = ["PH2", "PH3"]
phebnd4 = ["PH3", "PH4"]
phebnd5 = ["PH4", "PH1"]
phebnd6 = ["PH2", "PH4"]
spec_bond["PHE"] = [phebnd1, phebnd2, phebnd3, phebnd4, phebnd5];

tyrbnd1 = ["GBM", "TY1"]
tyrbnd2 = ["TY1", "TY2"]
tyrbnd3 = ["TY2", "TY3"]
tyrbnd4 = ["TY3", "TY4"]
tyrbnd5 = ["TY4", "TY1"]
tyrbnd6 = ["TY2", "TY4"]
spec_bond["TYR"] = [tyrbnd1, tyrbnd2, tyrbnd3, tyrbnd4, tyrbnd5];


###############################################################################
#  IDENTIFY RINGS DIHEDRALS
#  THEY CAN EASILY BE SPECIFIED HERE
######################################################
spec_dihed={}

# THIS IS NOT FULLY IMPLEMENTED.  THE CODE JUST ASSUMES
# THE FOUR CG SITE ARE PART OF THE DIHEDRAL...FINE FOR NOW
# YOU STILL NEED THESE LINE!!!!
spec_dihed["ADE"] = ["AD1", "AD2", "AD3", "AD4"]
spec_dihed["GUA"] = ["GU1", "GU2", "GU3", "GU4"]

spec_dihed["PHE"] = ["PH1", "PH2", "PH3", "PH4"]
spec_dihed["TYR"] = ["TY1", "TY2", "TY3", "TY4"]


#######################################################

# CHARGE HASH
charge={}
charge["PB"]  = -0.1118

charge["LY2"] = 0.1118
charge["AR2"] = 0.1118
charge["HI1"] = 0.00
charge["HI2"] = 0.00
charge["ASP"] = -0.1118
charge["GLU"] = -0.1118
charge["GBT"] = 0.1118
charge["ABT"] = 0.1118
# MASS HASH
bd_mass={}
bd_mass["RB1"] = 43.0454 
bd_mass["RB2"] = 56.0644 
bd_mass["DB2"] = 40.0650 
bd_mass["PB"]  = 87.1647 
bd_mass["AD1"] = 42.041 
bd_mass["AD2"] = 39.037 
bd_mass["AD3"] = 26.018 
bd_mass["AD4"] = 27.026 
bd_mass["GU1"] = 43.0254 
bd_mass["GU2"] = 39.037  
bd_mass["GU3"] = 26.018  
bd_mass["GU4"] = 42.041  
bd_mass["CY1"] = 42.041  
bd_mass["CY2"] = 43.0254 
bd_mass["CY3"] = 26.038  
bd_mass["TH1"] = 43.0254 
bd_mass["TH2"] = 43.0254 
bd_mass["TH3"] = 41.073  
bd_mass["UR1"] = 43.0254 
bd_mass["UR2"] = 43.0254 
bd_mass["UR3"] = 26.038  

bd_mass["GBB"] = 56.0385
bd_mass["GBM"] = 56.0385
bd_mass["GBT"] = 56.0385
bd_mass["ABB"] = 71.0732 
bd_mass["ABT"] = 71.0732 
bd_mass["VAL"] = 43.0883
bd_mass["LEU"] = 57.1151
bd_mass["ILE"] = 57.1151
bd_mass["PRO"] = 42.0804
bd_mass["MET"] = 75.1543

bd_mass["TR1"] = 39.0567
bd_mass["TR2"] = 38.0488
bd_mass["TR3"] = 53.0634
bd_mass["SER"] = 31.0287
bd_mass["THR"] = 45.0555
bd_mass["CYS"] = 46.0928
bd_mass["ASN"] = 58.0543
bd_mass["GLN"] = 72.0811
bd_mass["PH1"] = 33.5575
bd_mass["PH2"] = 19.5275 
bd_mass["PH3"] = 19.5275
bd_mass["PH4"] = 19.5275 
bd_mass["TY1"] = 33.5575
bd_mass["TY2"] = 19.5275
bd_mass["TY3"] = 35.5275
bd_mass["TY4"] = 19.5275

bd_mass["ASP"] = 58.0258
bd_mass["GLU"] = 72.0526
bd_mass["LY1"] = 42.0804
bd_mass["LY2"] = 31.0572
bd_mass["AR1"] = 42.0804
bd_mass["AR2"] = 59.0706

bd_mass["HI1"] = 26.0378
bd_mass["HI2"] = 27.0256
bd_mass["HI3"] = 28.0335

def read_pdb(cgpdb, pdb_list, cryst, ters):
    natom = 0
    try:
        f = open(cgpdb,"r")
    except:
        print ("ERROR: FILE",cgpdb,"IS NOT FOUND")
        sys.exit(0)
    line = f.readline()
    while line:
        recname = line[0:6].strip()
        if recname == "CRYST1":
            cryst.append(line.strip())
        elif recname == "ATOM":
            natom += 1
            index    = line[6:11]
            atmname = line[12:16]
            indicat  = line[16:17]
            resname =  line[17:21]
            chainid =  line[21:22]
            resid   =  line[22:26]
            code    =  line[26:27]
            posX = float(line[30:38])
            posY = float(line[38:46])
            posZ = float(line[46:54])
            occup = line[54:60]
            Tfact = line[60:66]
            segid = line[72:76]
            elesym = line[76:78].split("\n")[0]
            #charge = line[78:80]
            pdb_list.append([recname,index,atmname,indicat,resname,chainid,resid,code,posX,posY,posZ,occup,Tfact,segid,elesym])
        elif recname == "TER":
            ters.append(float(resid))
        line = f.readline()
    f.close()
    return natom

def open_file(outfile):                                                                                   
    try:
        fout = open(outfile,"w")
    except:
        print ("ERROR: FILE",outfile,"CANNOT OPEN")
        sys.exit(0)
    return fout

class gen_top_ENM:
    def __init__(self, cgpdb, aapdb, outfile, kENM, MAXdr, pspica, dssp):
        self.cgpdb  = cgpdb
        self.aapdb  = aapdb
        self.outfile = outfile
        self.kENM    = kENM
        self.MAXdr   = MAXdr
        self.pspica  = pspica
        self.dssp    = dssp
        self.nat     = 0
        self.nbb     = 0
        self.bbndx   = []
        self.iat2ibb = []
        self.bbtype  = []
        self.name    = []
        self.resname = []
        self.resid   = []
        self.bBackbone = []
        self.bPH1TY1   = []
        self.coord     = []
        self.pdb_data  = []
        self.ters      = []
        self.structure = []
        self.bb = ['GBM','GBB','GBT','ABB','ABT']
        self.helix  = ['H','G','I']
        self.sheet  = ['B','E','T']
        self.loop = ['S','C']
        cryst    = []
        self.natom = read_pdb(cgpdb, self.pdb_data, cryst, self.ters)
        self.ftop = open_file(outfile)
        self._charge_mod()
        self._set_array()
        self.read_dssp()

    def _charge_mod(self):
        if self.pspica:
            charge_unit   = 0.5590 
            charge["LY2"] = charge_unit
            charge["AR2"] = charge_unit
            charge["ASP"] = -charge_unit
            charge["GLU"] = -charge_unit
            charge["GBT"] = charge_unit
            charge["ABT"] = charge_unit

    def _set_array(self):
        for i in range(self.natom):
            tmp_name    = self.pdb_data[i][PDB_ATMNAME].strip()
            #tmp_name    = self.pdb_data[i][PDB_ATMNAME].strip()[:2]
            tmp_resname = self.pdb_data[i][PDB_RESNAME].strip()
            if tmp_name in ["GBB", "GB", "ABB", "GBM", "PB", "RB1", "RB2", "DB2"]:
                if tmp_resname in ["ALA", "ADE", "GUA", "CYT", "THY", "URA"]:
                    self.name.append(tmp_name)
                elif tmp_resname == "GLY":
                    self.name.append("GBB")
                else:
                    self.name.append("GBM")
                self.bbtype.append(tmp_resname)
                self.resname.append(tmp_resname)
                self.bbndx.append(self.nat)
                self.resid.append(float(self.pdb_data[i][PDB_RESID]))
                self.coord.append([self.pdb_data[i][PDB_POSX], self.pdb_data[i][PDB_POSY], self.pdb_data[i][PDB_POSZ]])
                self.iat2ibb.append(self.nbb)
                self.bBackbone.append(1)
                self.bPH1TY1.append(0)
                self.nat += 1; self.nbb += 1
                if tmp_resname not in ["GLY", "ALA", "ADE", "GUA", "CYT", "THY", "URA"]:
                    for j in res_map[tmp_resname]:
                        self.name.append(j)
                        self.coord.append([self.pdb_data[self.nat][PDB_POSX], self.pdb_data[self.nat][PDB_POSY], self.pdb_data[self.nat][PDB_POSZ]])
                        self.resname.append(tmp_resname)
                        self.resid.append(float(self.pdb_data[i][PDB_RESID]))
                        self.iat2ibb.append(None)
                        self.bBackbone.append(0)
                        if j == "PH2" or j == "TY2" or j == "PH3" or j == "TY3" or j == "PH4" or j == "TY4":
                            self.bPH1TY1.append(1)
                        else :
                            self.bPH1TY1.append(0)
                        self.nat += 1
                if tmp_name in ["RB2", "DB2"]:
                    for j in res_map[tmp_resname]:
                        self.name.append(j)
                        self.coord.append([self.pdb_data[self.nat][PDB_POSX], self.pdb_data[self.nat][PDB_POSY], self.pdb_data[self.nat][PDB_POSZ]])
                        self.resname.append(tmp_resname)
                        self.resid.append(float(self.pdb_data[i][PDB_RESID]))
                        self.iat2ibb.append(None)
                        self.bBackbone.append(0)
                        if j == "AD2" or j == "GU2" or j == "AD3" or j == "GU3" or j == "AD4" or j == "GU4":
                            self.bPH1TY1.append(1)
                        else :
                            self.bPH1TY1.append(0)
                        self.nat += 1
        resid_0 = self.resid[0]
        for i in range(len(self.resid)):
            self.resid[i] = self.resid[i] - resid_0
            print(self.resid[i])

    def read_dssp(self):
        aapdb = self.aapdb
        dssp = self.dssp
        cmd = f"{dssp} {aapdb} dssp.out"
        runcmd = subprocess.call(cmd.split())
        if runcmd != 0:
            print ("ERROR: CANNOT MAKE DSSP FILE")
            sys.exit(0)
        f = open('dssp.out','r')
        lines = f.readlines()
        f.close()
        for line in lines[28:]:
            if line[16] == ' ':
                self.structure.append('C')
            else:
                self.structure.append(line[16])
        self.structure.append('C')

    ######################################################
    #########        PRINT OUT TOP FILE      #############
    ######################################################
    def write_atom(self):
        ftop = self.ftop
        for i in range(self.nat):
            name_i    = self.name[i]
            resname_i = self.resname[i]
            I = i + 1
            if name_i in charge:
                thischrg = charge[name_i]
            else :
                thischrg = 0.0
        
            wtype = this_type[name_i]
            mass  = bd_mass[name_i]

            if resname_i not in ["ADE", "GUA", "CYT", "THY", "URA"]:
                if self.bbndx[0] == i:
                    # backbone N-terminal
                    thischrg = charge["GBT"]
                    if resname_i == "ALA":
                        if self.pspica and thischrg > 0:
                            print ("atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
                                %(I,resname_i,"ABTP","ABTP",mass,thischrg), file=ftop)
                        else:
                            print ("atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
                                %(I,resname_i,"ABT","ABT",mass,thischrg), file=ftop)
                    else :
                        if self.pspica and thischrg > 0:
                            print ("atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
                                %(I,resname_i,"GBTP","GBTP",mass,thischrg), file=ftop)
                        else:
                            print ("atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
                                %(I,resname_i,"GBT","GBT",mass,thischrg), file=ftop)
                elif self.bbndx[self.nbb-1] == i:
                    # backbone C-terminal
                    thischrg = -1*charge["GBT"]
                    if resname_i == "ALA":
                        if self.pspica and thischrg < 0:
                            print ("atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
                                %(I,resname_i,"ABTN","ABTN",mass,thischrg), file=ftop)
                        else:
                            print ("atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
                                %(I,resname_i,"ABT","ABT",mass,thischrg), file=ftop)
                    else :
                        if self.pspica and thischrg < 0:
                            print ("atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
                                %(I,resname_i,"GBTN","GBTN",mass,thischrg), file=ftop)
                        else:
                            print ("atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
                                %(I,resname_i,"GBT","GBT",mass,thischrg), file=ftop)
                else :
                    print ("atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
                        %(I,resname_i,name_i,wtype,mass,thischrg), file=ftop)
            else :
                print ("atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
                    %(I,resname_i,name_i,wtype,mass,thischrg), file=ftop)
            
        print ("", file=ftop)

    ##########################
    # SET UP BONDS NOW
    ##########################
    
    #######################################
    # Elastic network model for Backbone
    ########################################
    def write_ENM(self):
        ftop  = self.ftop
        MAXdr = self.MAXdr
        MINdr = 0.0
        kENM  = self.kENM
        print ()
        print ("******** Elastic network info ********")
        print ("Cutoff distance ... ", MAXdr, "A")
        print ("Force constant ... ", kENM, "kcal/A^2")
        print ("**************************************")
        print ()
        nat = self.nat
        for i1 in range(nat): 
            I1 = i1 + 1
            for i2 in range(i1 + 1,nat):
                I2 = i2 + 1
                if self.bBackbone[i1] == 1 and self.bBackbone[i2] == 1:
                    dx = self.coord[i1][0] - self.coord[i2][0]
                    dy = self.coord[i1][1] - self.coord[i2][1]
                    dz = self.coord[i1][2] - self.coord[i2][2]
                    dr = math.sqrt(dx*dx + dy*dy + dz*dz)
                    if dr < MAXdr and dr > MINdr and self.resid[i2] - self.resid[i1] >= 3:
                        print ("bondparam %5d %5d   %f %f # %s-%s" \
                            %(I1, I2, kENM, dr, self.name[i1], self.name[i2]), file=ftop)
        print ("", file=ftop)

    #######################################
    # PRINT THE BACKBONE BONDS FIRST 
    #######################################
    def write_bond_BB(self):
        ftop  = self.ftop
        self.total_bonds = 0
        bond_index1 = []
        bond_index2 = []
        for i in range(self.nbb - 1):
            if self.resid[self.bbndx[i]] > self.resid[self.bbndx[i+1]]:
                continue
            if self.resid[self.bbndx[i]] in self.ters:
                continue
            self.total_bonds += 1
            print ("bond %5d %5d # %s-%s" \
                %(self.bbndx[i]+1, self.bbndx[i+1]+1, self.name[self.bbndx[i]], self.name[self.bbndx[i+1]]), file=ftop)
            bond_index1.append(self.bbndx[i])
            bond_index2.append(self.bbndx[i+1])
        print ("", file=ftop)
        self.bond_index1 = bond_index1
        self.bond_index2 = bond_index2

    #######################################
    # DO SIDECHAINS. SPEC_BOND ABOVE CAN BE USED 
    # TO DESCRIBE SPECIAL BONDS THAT WOULD BE
    # DIFFICULT TO GUESS 
    #######################################
    def write_bond_SC(self):
        ftop  = self.ftop
        for i in range(self.nbb):
            bbtype_i = self.bbtype[i]
            bbndx_i  = self.bbndx[i]
            if self.name[bbndx_i] in ["PB", "RB1"]:
                continue
            if bbtype_i not in ["GLY", "ALA"]:
                # check to see if it is a special bond
                if bbtype_i in spec_bond:
                    for j in spec_bond[bbtype_i]:
                        bndx1 = 0
                        bndx2 = 0
                        for k in j:
                            if bndx1 == 0:
                                tmp_ndx = bbndx_i
                                while self.name[tmp_ndx] != k:
                                    tmp_ndx += 1
                                bndx1 = tmp_ndx
                            else :
                                tmp_ndx = bbndx_i
                                while self.name[tmp_ndx] != k:
                                    tmp_ndx += 1
                                bndx2 = tmp_ndx
                        print ("bond %5d %5d # %s-%s" \
                                %(bndx1+1, bndx2+1, self.name[bndx1], self.name[bndx2]), file=ftop)
                        self.total_bonds += 1
                        self.bond_index1.append(bndx1)
                        self.bond_index2.append(bndx2)
                else :
                    print ("bond %5d %5d # %s-%s" \
                        %(bbndx_i+1, bbndx_i+1+1, self.name[bbndx_i], self.name[bbndx_i+1]), file=ftop)
                    self.total_bonds += 1
                    self.bond_index1.append(bbndx_i)
                    self.bond_index2.append(bbndx_i+1)
        print ("", file=ftop)

    #########################################
    # NOW TAKE CARE OF ANGLES IF WANTED
    # GENERATE EVERY COMBINATION FROM THE
    # BOND LIST UNLESS SOME ARE NOT SPECIFIED 
    # ON THE COMMAND LINE.
    #########################################
    def write_angle(self):
        ftop  = self.ftop
        name        = self.name
        bond_index1 = self.bond_index1
        bond_index2 = self.bond_index2
        bBackbone   = self.bBackbone
        bPH1TY1     = self.bPH1TY1
        for i1 in range(self.total_bonds):
            for i2 in range(i1 + 1, self.total_bonds):
                if bond_index1[i1] == bond_index1[i2]:
                    andx1 = bond_index2[i1]
                    andx2 = bond_index1[i1]
                    andx3 = bond_index2[i2]
                    # GB-GB-GB is non zero.
                    # GB-GB-SC is non zero.
                    # GB-SC-SC is non zero.
                    # SC-TY1-SC is non zero.
                    # SC-PH1-SC is non zero.
                    # SC-SC-SC is zero.
                    if bBackbone[andx1] + bBackbone[andx2] + bBackbone[andx3] != 0 or bPH1TY1[andx2] == 1:
                        if name[andx2] in ["PH3","TY3","AD3","GU3"]:
                            print("angleparam %5d %5d %5d  0.0  90.0000 # %s %s %s"\
                                %(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3]), file=ftop)
                        elif name[andx2] in ["PH2","PH4","TY2","TY4"]:
                            print("angle      %5d %5d %5d               # %s %s %s" \
                                %(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3]), file=ftop)
                        else:
                            r1 = np.array(self.coord[andx1])
                            r2 = np.array(self.coord[andx2])
                            r3 = np.array(self.coord[andx3])
                            angle_in_pdb = 180.0/np.pi*get_angle(r1,r2,r3)
                            print("angleparam %5d %5d %5d  -1  %8.4f # %s %s %s"\
                                %(andx1+1,andx2+1,andx3+1,angle_in_pdb,name[andx1],name[andx2],name[andx3]),file=ftop)
                elif bond_index1[i1] == bond_index2[i2]:
                    andx1 = bond_index2[i1]
                    andx2 = bond_index1[i1]
                    andx3 = bond_index1[i2]
                    if bBackbone[andx1] + bBackbone[andx2] + bBackbone[andx3] != 0 or bPH1TY1[andx2] == 1:
                        if name[andx2] in ["PH3","TY3","AD3","GU3"]:
                            print("angleparam %5d %5d %5d  0.0  90.0000 # %s %s %s"\
                                %(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3]), file=ftop)
                        elif name[andx2] in ["PH2","PH4","TY2","TY4"]:
                            print("angle      %5d %5d %5d               # %s %s %s" \
                                %(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3]), file=ftop)
                        else:
                            r1 = np.array(self.coord[andx1])
                            r2 = np.array(self.coord[andx2])
                            r3 = np.array(self.coord[andx3])
                            angle_in_pdb = 180.0/np.pi*get_angle(r1,r2,r3)
                            print("angleparam %5d %5d %5d  -1  %8.4f # %s %s %s"\
                                %(andx1+1,andx2+1,andx3+1,angle_in_pdb,name[andx1],name[andx2],name[andx3]),file=ftop)
                elif bond_index2[i1] == bond_index1[i2]:
                    andx1 = bond_index1[i1]
                    andx2 = bond_index2[i1]
                    andx3 = bond_index2[i2]
                    if bBackbone[andx1] + bBackbone[andx2] + bBackbone[andx3] != 0 or bPH1TY1[andx2] == 1:
                        if name[andx2] in ["PH3","TY3","AD3","GU3"]:
                            print("angleparam %5d %5d %5d  0.0  90.0000 # %s %s %s"\
                                %(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3]), file=ftop)
                        elif name[andx2] in ["PH2","PH4","TY2","TY4"]:
                            print("angle      %5d %5d %5d               # %s %s %s" \
                                %(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3]), file=ftop)
                        else:
                            r1 = np.array(self.coord[andx1])
                            r2 = np.array(self.coord[andx2])
                            r3 = np.array(self.coord[andx3])
                            angle_in_pdb = 180.0/np.pi*get_angle(r1,r2,r3)
                            print("angleparam %5d %5d %5d  -1  %8.4f # %s %s %s"\
                                %(andx1+1,andx2+1,andx3+1,angle_in_pdb,name[andx1],name[andx2],name[andx3]),file=ftop)
                elif bond_index2[i1] == bond_index2[i2]:
                    andx1 = bond_index1[i1]
                    andx2 = bond_index2[i1]
                    andx3 = bond_index1[i2]
                    if bBackbone[andx1] + bBackbone[andx2] + bBackbone[andx3] != 0 or bPH1TY1[andx2] == 1:
                        if name[andx2] in ["PH3","TY3","AD3","GU3"]:
                            print("angleparam %5d %5d %5d  0.0  90.0000 # %s %s %s"\
                                %(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3]), file=ftop)
                        elif name[andx2] in ["PH2","PH4","TY2","TY4"]:
                            print("angle      %5d %5d %5d               # %s %s %s" \
                                %(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3]), file=ftop)
                        else:
                            r1 = np.array(self.coord[andx1])
                            r2 = np.array(self.coord[andx2])
                            r3 = np.array(self.coord[andx3])
                            angle_in_pdb = 180.0/np.pi*get_angle(r1,r2,r3)
                            print("angleparam %5d %5d %5d  -1  %8.4f # %s %s %s"\
                                %(andx1+1,andx2+1,andx3+1,angle_in_pdb,name[andx1],name[andx2],name[andx3]),file=ftop)
        print ("", file=ftop)

    ######################################
    #  WRITE OUT ANY DIHEDRALS FOR RINGS
    ######################################
    def write_dihedral(self):
        ftop  = self.ftop
        for i in range(self.nbb):
            bbtype_i = self.bbtype[i]
            bbndx_i  = self.bbndx[i]
            name_i   = self.name[bbndx_i]
            #MAKE SURE WE HAVE A SIDECHAIN
            if bbtype_i != "GLY" and bbtype_i != "ALA":
            #CHECK TO SEE IF IT IS A SPECIAL DIHEDRAL
                if bbtype_i in spec_dihed and name_i not in ["PB", "RB1"]:
                    ndx1 = bbndx_i + 1
                    ndx2 = ndx1 + 1
                    ndx3 = ndx2 + 1
                    ndx4 = ndx3 + 1
                    bbndx_i = self.bbndx[i]
                    print ("dihedralparam  %5d %5d %5d %5d   50.0 1 180 0.0 # %s-%s-%s-%s" \
                        %(ndx1+1, ndx2+1, ndx3+1, ndx4+1, 
                          self.name[bbndx_i+1],self.name[bbndx_i+2],self.name[bbndx_i+3],self.name[bbndx_i+4]), file=ftop)

    def finalize(self):
        self.ftop.close()
        if self.pspica:
            print ('## Protein topology file has been generated for pSPICA FF ##')
            print ('NOTE: Remove "-pspica" option when you use generated protein topology files for SPICA FF.')
        else:
            print ('## Protein topology file has been generated for SPICA FF ##')
            print ('NOTE: Use "-pspica" option when you use generated protein topology files for pSPICA FF.')
        print ()

    def run(self):
        self.write_atom()
        self.write_ENM()
        self.write_bond_BB()
        self.write_bond_SC()
        self.write_angle()
        self.write_dihedral()
        self.finalize()
     
########################################################
# MAIN CODE STARTS HERE..... THE PDB FILE WAS READ ABOVE
# START BY LOOPING OVER EACH RESIDUE (DOMAIN)
########################################################
if __name__ == "__main__":
    args = get_option()
    cgpdb   = args.cgpdb
    aapdb   = args.aapdb
    outfile = args.output
    kENM    = args.kENM
    MAXdr   = args.maxr
    pspica  = args.pspica

    gen = gen_top_ENM(cgpdb, aapdb, outfile, kENM, MAXdr, pspica)
    gen.run()
