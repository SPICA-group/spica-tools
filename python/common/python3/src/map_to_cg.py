#!/usr/bin/env python
# original perl script was made
# BY RUSSELL DEVANE 2/2008
# and modified 
# BY SHUHEI KAWAMOTO
# and converted 
# BY YUSUKE MIYAZAKI 2019/10
# and converted to python3
# BY ISSEI KAWABATA 2022/4
###############################################################################
# MAPS AND ALL ATOM PROTEIN/PEPTIDE/LIPID IN PDB FORMAT
# TO A CG MODEL AND OUTPUTS A PDB FILE
################################################################################
# MOLECULE INFO. OTHER THAN PROTEIN/PEPTIDE IS LOADED FROM JSON FORMAT
################################################################################
#
# USAGE: map_to_cg.py <aa.pdb filename> < cg.pdb filename (output)> 
#
#################################################################################

import sys, json, random
from pathlib import Path
import numpy as np
from argparse import ArgumentParser


def open_file(outfile):
    try:
        fout = open(outfile,"w")
    except:
        print ("ERROR: FILE",outfile,"CANNOT OPEN")
        sys.exit(0)
    return fout

def get_option():
    json = Path(__file__).parents[0] / "spica_top.json"
    argparser = ArgumentParser()
    argparser.add_argument('input', type=str,
                            help='Specify input AA PDB file name.')
    argparser.add_argument('output', type=str,
                            help='Specify output CG PDB file name.')
    argparser.add_argument('-json', type=str,
                            default=json,
                            help='input json file name (default: spica_top.json).')
    argparser.add_argument('-nodelwat', action='store_true',
                            help='not delete excess water due to CG ion mapping (default: off).')
    argparser.add_argument('-verbose', type=int,
                            default=1,
                            help='activate verbose logging, 0 : off, 1 : on (default: 1).')
    return argparser.parse_args()

def get_option_script(argv):
    json = Path(__file__).parents[0] / "spica_top.json"
    argparser = ArgumentParser(usage='map2cg [-h] [-json JSON] [-nodelwat] [-verbose VERBOSE] input output',
                               prog ="map2cg")
    argparser.add_argument('input', type=str,
                            help='Specify input AA PDB file name.')
    argparser.add_argument('output', type=str,
                            help='Specify output CG PDB file name.')
    argparser.add_argument('-json', type=str,
                            default=json,
                            help='input json file name (default: spica_top.json).')
    argparser.add_argument('-nodelwat', action='store_true',
                            help='not delete excess water due to CG ion mapping (default: off).')
    argparser.add_argument('-verbose', type=int,
                            default=1,
                            help='activate verbose logging, 0 : off, 1 : on (default: 1).')
    return argparser.parse_args(argv)

ncomp     = 7
ncomp_pdb = 15
(   PDB_RECNAME,
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
(   ATMNAME,
    RESNAME,
    RESID,
    CHAINID,
    POSX, 
    POSY,
    POSZ,
) = range(ncomp)

atm_mass = {    
        "H":1.008000,
        "C":12.010700,
        "N":14.006740,
        "O":15.999400,
        "P":30.973762,
        "S":32.066000
       }
res_name = {}
res_map  = {}    
res_map["ALA"] = [ [ "CB"] ]
res_map["VAL"] = [ [ "CB", "CG1", "CG2"] ]
res_map["LEU"] = [ [ "CB", "CD1", "CD2", "CG" ] ]
res_map["ILE"] = [ [ "CB", "CD", "CG1", "CG2" ] ]
res_map["PRO"] = [ [ "CB", "CD", "CG" ] ]
res_map["MET"] = [ [ "CB", "CG", "SD", "CE" ] ] 
res_map["SER"] = [ [ "CB", "OG" ] ]
res_map["THR"] = [ [ "CB", "CG2", "OG1" ] ]
res_map["CYS"] = [ [ "CB", "SG" ] ]
res_map["ASN"] = [ [ "CB", "CG", "ND2", "OD1" ] ]
res_map["GLN"] = [ [ "CB", "CD", "CG", "NE2", "OE1" ] ]
res_map["ASP"] = [ [ "CB", "CG", "OD1", "OD2" ] ]
res_map["GLU"] = [ [ "CB", "CD", "CG", "OE1", "OE2" ] ]
PH1 = ["CB", "CG"]
PH2 = ["CE1", "CD1"]
PH3 = ["CZ"]
PH4 = ["CD2", "CE2"]
res_map["PHE"] = [ PH1, PH2, PH3, PH4 ]
TR1 = ["CB", "CG", "CD1"]
TR2 = ["CD2", "CE3", "CZ3"]
TR3 = ["CH2", "CZ2", "CE2", "NE1"]
res_map["TRP"] = [ TR1, TR2, TR3 ]
TY1 = ["CB",  "CG"]
TY2 = ["CD1", "CE1"]
TY3 = ["CZ", "OH"]
TY4 = ["CD2", "CE2"]
res_map["TYR"] = [ TY1, TY2, TY3, TY4 ]
LY1 = ["CB", "CD", "CG"]
LY2 = ["CE", "NZ"]
res_map["LYS"] = [ LY1, LY2 ]
AR1 = ["CB", "CD", "CG"]
AR2 = ["CZ", "NE", "NH1", "NH2"]
res_map["ARG"] = [ AR1, AR2 ]
HI1 = ["CB", "CG"]
HI2 = ["CD2", "NE2"]
HI3 = ["ND1", "CE1"]
res_map["HIS"] = [ HI1, HI2, HI3 ]
res_map["HSP"] = [ HI1, HI2, HI3 ]
res_map["HSD"] = [ HI1, HI2, HI3 ]
res_map["HSE"] = [ HI1, HI2, HI3 ]
# SPECIAL RESIDUES
res_map["SERP"] = [ [ "CB", "OG" ], [ "P", "O1P", "O2P", "OT" ] ]
res_map["DCYS"] = [ [ "CB", "SG" ] ]
res_map["DLYS"] = [ LY1, LY2 ]
res_map["DARG"] = [ AR1, AR2 ]

prot_list = ["GLY","ALA","VAL","LEU","ILE","PRO","MET","PHE","TRP","SER","THR","CYS","ASN","GLN","TYR","ASP","GLU","LYS","ARG","HIS",
         "HSP","HSD","HSE",
         "SERP","DCYS","DLYS","DARG"
         ]
### END OF SIDECHAIN HASH 
### BACKBONE MAPPING FOR GENERAL AND ALANINE
GBB = ["CA"]
ABB = ["CA"]

######################################################
def read_json(jsonfile,mol_list,res_name,res_map):
    try:
        f = open(jsonfile,"r")
    except:
        print ('ERROR: "%s" cannot be opened.' % jsonfile)
        sys.exit(0)
    else:
        jsn = json.load(f)
        f.close()
    for ir in jsn["topo"].keys():
        tmp = []
        mol_list.append(ir)
        res_name[ir] = jsn["topo"][ir]["name"]
        nname = len(res_name[ir])
        for ia in range(nname):
            imap = jsn["topo"][ir]["map"][ia]
            if type(imap) is not list:
                imap = [imap]
            tmp.append(imap)
        res_map[ir] = tmp

def make_domain(pdb_list,domain_list):
    nl = len(pdb_list)
    resid0 = None
    nres = 0
    tmp_domain_list = [[],[],[],[],[],[],[]]
    for i in range(nl):
        if resid0 is None:
            resid0 = pdb_list[i][PDB_RESID]
            nres += 1
        if resid0 == pdb_list[i][PDB_RESID]:
            tmp_domain_list[ATMNAME].append(pdb_list[i][PDB_ATMNAME].strip())
            tmp_domain_list[RESNAME].append(pdb_list[i][PDB_RESNAME].strip())
            tmp_domain_list[RESID].append(pdb_list[i][PDB_RESID])
            tmp_domain_list[CHAINID].append(pdb_list[i][PDB_CHAINID])
            tmp_domain_list[POSX].append(pdb_list[i][PDB_POSX])
            tmp_domain_list[POSY].append(pdb_list[i][PDB_POSY])
            tmp_domain_list[POSZ].append(pdb_list[i][PDB_POSZ])
        else:
            domain_list.append(tmp_domain_list)
            tmp_domain_list = [[],[],[],[],[],[],[]]
            resid0 = pdb_list[i][PDB_RESID]
            nres += 1
            tmp_domain_list[ATMNAME].append(pdb_list[i][PDB_ATMNAME].strip())
            tmp_domain_list[RESNAME].append(pdb_list[i][PDB_RESNAME].strip())
            tmp_domain_list[RESID].append(pdb_list[i][PDB_RESID])
            tmp_domain_list[CHAINID].append(pdb_list[i][PDB_CHAINID])
            tmp_domain_list[POSX].append(pdb_list[i][PDB_POSX])
            tmp_domain_list[POSY].append(pdb_list[i][PDB_POSY])
            tmp_domain_list[POSZ].append(pdb_list[i][PDB_POSZ])
    
    domain_list.append(tmp_domain_list)
    return nres

def read_pdb(infile):
    pdb_list = []
    cryst    = []
    natom    = 0
    try:
        f = open(infile,"r")
    except:
        print ("ERROR: FILE",infile,"IS NOT FOUND")
        sys.exit(0)
    line = f.readline()
    while line:
      recname = line[0:6].strip()
      if recname == "CRYST1":
        cryst.append(line.strip())
      elif recname == "ATOM":
        natom += 1
        index   = line[6:11]
        atmname = line[12:16]
        indicat = line[16:17]
        resname = line[17:21]
        chainid = line[21:22]
        resid   = line[22:26]
        code    = line[26:27]
        posX   = float(line[30:38])
        posY   = float(line[38:46])
        posZ   = float(line[46:54])
        occup  = line[54:60]
        Tfact  = line[60:66]
        segid  = line[72:76]
        elesym = line[76:78].split("\n")[0]
        #charge = line[78:80]
        pdb_list.append([recname,index,atmname,indicat,resname,chainid,resid,code,posX,posY,posZ,occup,Tfact,segid,elesym])
      line = f.readline()
    f.close()
    return natom, pdb_list, cryst

def write_domain(pdb_domain,cryst,f):
    nl = len(pdb_domain)
    cnt = 1
    recname = "ATOM"
    indicat = " "
    code = " "
    chainid = " "
    print ('%s' % "CG MAP",file=f)
    print ('%s' % "REMARK Generated by map_to_cg.py",file=f)
    if len(cryst) == 1:
        print ('%s' % cryst[0],file=f)
    for i in range(nl):
        atmname  = pdb_domain[i][ATMNAME]
        resname  = pdb_domain[i][RESNAME]
        chainid  = pdb_domain[i][CHAINID]
        resid    = pdb_domain[i][RESID]
        posX     = pdb_domain[i][POSX]
        posY     = pdb_domain[i][POSY]
        posZ     = pdb_domain[i][POSZ]
        index    = str(cnt)
        print ('%-6s%5s %4s%1s%-4s%1s%4s%1s   %8.3f%8.3f%8.3f' \
                    %(recname,index[0:5],atmname,indicat,resname,chainid,resid,code,posX,posY,posZ), file=f)
        cnt += 1
    print ('%s' % "END", file=f)

########################################################
# MAIN CODE STARTS HERE..... THE PDB FILE WAS READ ABOVE
# START BY LOOPING OVER EACH RESIDUE (DOMAIN)
########################################################
class map_to_cg:
    def __init__(self, infile, outfile, jsonfile, nodelwat, verbose):
        self.infile            = infile
        self.outfile           = outfile
        self.jsonfile          = jsonfile
        self.nodelwat          = nodelwat
        self.nodelwat          = nodelwat
        self.verbose           = verbose
        self.mol_list          = []
        self.cryst             = []
        self.domain_data       = []
        self.res_name_j        = {}
        self.res_map_j         = {}
        self.all_cg_pdb_domain = []
        self.wat_index         = []
        self.nwat              = 0
        self.nsod              = 0
        self.ncla              = 0
        self.cg_domain_count   = 0
        self._counter          = 0
        self.natom, self.pdb_data, self.cryst = read_pdb(infile)
        self.nres              = make_domain(self.pdb_data, self.domain_data)
        read_json(jsonfile, self.mol_list, self.res_name_j, self.res_map_j)
        self._write_init_info()
    
    def initialize(self):
        self.all_cg_pdb_domain = []
        self.cg_domain_count   = 0
        self.nwat              = 0

    def _write_init_info(self):
        print ("# of ATOMS:", self.natom)
        print ("# of RESIDUES:", self.nres)
    
    def _weight_pos(self, aa_domain, aatype, mass, com): 
        l = aa_domain[ATMNAME].index(aatype)
        x = float(aa_domain[POSX][l])
        y = float(aa_domain[POSY][l])
        z = float(aa_domain[POSZ][l])
        r = np.array([x, y, z])
        this_mass = atm_mass[aatype[0][0]]
        mass += this_mass
        com  += this_mass*r
        return mass, com

    def _com(self, aa_domain, aatypes, from_mol_list=False):
        mass = 0.0
        com  = np.zeros(3)
        if from_mol_list:
            if set(aatypes) == set(list(set(aatypes) & set(aa_domain[ATMNAME]))):
                if set(aatypes) == set(["OH2","H1","H2"]) or set(aatypes) == set(["OW","HW1","HW2"]):
                    self.nwat += 1 
                    if (self.nwat - 1) % 3 != 0: 
                        return np.zeros(3), "SKIP"
                for aatype in aatypes:
                    mass, com = self._weight_pos(aa_domain, aatype, mass, com)
            else :
                sys.exit ("ERROR:", aa_domain[RESNAME][0],"ATOM",aatypes,"IS MISSING ON MOLECULE")
        else:
            for aatype in aatypes:
                if aa_domain[ATMNAME].count(aatype):
                    mass, com = self._weight_pos(aa_domain, aatype, mass, com)
                else :
                    bad_res = aa_domain[RESID][0]
                    print ("ERROR:", aa_domain[RESNAME][0],"BACKBONE ATOM",aatypes,"IS MISSING ON RESIDUE",bad_res)
                    return np.zeros(3), "UNK"
        com /= mass
        return com, aa_domain[RESNAME][0]
    
    def _set_cg_domain(self, aa_domain, TYPE, bead_name, com):
        cg_pdb_domain = ["UNK","UNK",0,"A",0.0,0.0,0.0]
        cg_pdb_domain[CHAINID] = aa_domain[CHAINID][0]
        cg_pdb_domain[RESID]   = aa_domain[RESID][0]
        cg_pdb_domain[RESNAME] = TYPE; cg_pdb_domain[ATMNAME] = bead_name
        cg_pdb_domain[POSX] = com[0];  cg_pdb_domain[POSY] = com[1];  cg_pdb_domain[POSZ] = com[2]
        self.all_cg_pdb_domain.append([cg_pdb_domain[m] for m in range(ncomp)])
        self.cg_domain_count += 1
    
    def _bead_name(self, TYPE):
        if TYPE != "ALA":
            if TYPE != "GLY":
                return "GBM" 
            else:
                return "GBB"
        else:
            return "ABB"

    def _bead_numbering(self, bead_name, idx):
        bead_name = bead_name.strip()[:2]
        bead_name = "%s%d" % (bead_name, idx+1)
        return bead_name

    def _warn_TYPE(self, TYPE):
        print(f"WARNING: RESIDUE {TYPE} IS NOT DEFINED in RESMAP.")

    def error_TYPE(self, TYPE):
        print(f"ERROR: RESIDUE {TYPE} IS NOT DEFINED in MOL_LIST.")
        sys.exit(0)

    def finalize(self):
        fout = open_file(self.outfile)
        write_domain(self.all_cg_pdb_domain, self.cryst, fout)
        fout.close()
        if self.verbose == 1:
            print ()
            print ("Normal Termination.")

    def prot_map(self, aa_domain):
        TYPE = aa_domain[RESNAME][0]
        if TYPE != "ALA" :
            # CREATE A NEW DOMAIN FOR THIS RESIDUE BACKBONE
            # AND STORE THE COM OF THE BACKBONE ATOMS
            com, _    = self._com(aa_domain, GBB)
            bead_name = self._bead_name(TYPE)
            self._set_cg_domain(aa_domain, TYPE, bead_name, com)
            # Just finished the backbone
            # NOW TAKE CARE OF THE SIDECHAINS
            if TYPE != "GLY" :
                if TYPE in res_map: 
                    # LOOP OVER EACH SIDECHAIN BEAD HERE
                    for idx, res in enumerate(res_map[TYPE]):
                        com, bead_name  = self._com(aa_domain, res)
                        # CHECK TO SEE IF WE HAVE MULTIPLE BEAD AND NEED TO ADD INDEX TO SIDECHAIN BEAD NAME
                        if len(res_map[TYPE]) > 1:
                            bead_name = self._bead_numbering(bead_name, idx)
                        self._set_cg_domain(aa_domain, TYPE, bead_name, com)
                else:
                    self._warn_TYPE(TYPE)
        else:
            com, _    = self._com(aa_domain, ABB)
            bead_name = self._bead_name(TYPE)
            self._set_cg_domain(aa_domain, TYPE, bead_name, com)
    
    def mol_map(self, aa_domain):
        TYPE = aa_domain[RESNAME][0]
        if TYPE in self.res_map_j: 
            for idx, res in enumerate(self.res_map_j[TYPE]):
                com, ret  = self._com(aa_domain, res, True)
                if ret == "SKIP":
                    continue
                bead_name = self.res_name_j[TYPE][idx]
                self._set_cg_domain(aa_domain, TYPE, bead_name, com)
                if self._counter == 0:
                    if bead_name == "SOD":
                        self.nsod += 1
                    if bead_name == "CLA":
                        self.ncla += 1
                    if bead_name == "WAT" or bead_name == "W":
                        self.wat_index.append(self.cg_domain_count)
        else:
            self._warn_TYPE(TYPE)
    
    def delete_ex_wat(self):
        if self._counter == 0:
            del_nwat = (self.nsod*3 + self.ncla*2) // 3
            self.del_wat_index = sorted(random.sample(self.wat_index, del_nwat), reverse=True)
        for idx in self.del_wat_index:
            self.all_cg_pdb_domain.pop(idx)

    def mapping(self):
        if self.verbose == 1:
            print ("**RESIDUES IN AA PDB**")
        for i in range(self.nres):
            aa_domain = self.domain_data[i]
            TYPE = aa_domain[RESNAME][0]
            if self.verbose == 1:
                print (TYPE)
            if prot_list.count(TYPE):
                self.prot_map(aa_domain)
            elif self.mol_list.count(TYPE):
                self.mol_map(aa_domain)
            else:
                self.error_TYPE(TYPE)
        if not self.nodelwat:
            self.delete_ex_wat()
        self._counter += 1

    def run(self):
        self.mapping()
        self.finalize()

if __name__ == "__main__":
    args     = get_option()
    infile   = args.input
    outfile  = args.output
    jsonfile = args.json
    nodelwat = args.nodelwat
    verbose  = args.verbose

    mapCG = map_to_cg(infile, outfile, jsonfile, nodelwat, verbose)
    mapCG.run()
