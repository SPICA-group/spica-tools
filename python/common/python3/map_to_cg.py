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

import sys
import json

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

atm_mass = {    "H":1.008000,
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
def read_json(infile,mol_list,res_name,res_map):
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
######################################################
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

#######################################################
def read_pdb(infile,pdb_list,cryst):
    natom = 0
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
      line = f.readline()
    f.close()
    return natom

#######################################################
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
        #chainid  = pdb_domain[i][CHAINID]
        resid    = pdb_domain[i][RESID]
        posX     = pdb_domain[i][POSX]
        posY     = pdb_domain[i][POSY]
        posZ     = pdb_domain[i][POSZ]
        index    = str(cnt)
        print ('%-6s%5s %4s%1s%4s%1s%4s%1s   %8.3f%8.3f%8.3f' \
                    %(recname,index[0:5],atmname,indicat,resname,chainid,resid,code,posX,posY,posZ), file=f)
        cnt        += 1
    print ('%s' % "TER", file=f)
    print ('%s' % "END", file=f)

########################################################
# MAIN CODE STARTS HERE..... THE PDB FILE WAS READ ABOVE
# START BY LOOPING OVER EACH RESIDUE (DOMAIN)
########################################################
if __name__ == "__main__":
    args = sys.argv
    if len(args) != 3:
        print ("USAGE: map_to_cg.py <aa.pdb filename> <cg.pdb filename (output)>")
        sys.exit(0)
    infile  = args[1]
    outfile = args[2]
    jsonfile = "spica_top.json"
    mol_list          = []
    res_name_j          = {}
    res_map_j          = {}
    cryst             = []
    pdb_data          = []
    domain_data       = []
    all_cg_pdb_domain = []
    cg_pdb_domain     = ["UNK","UNK",0,"A",0.0,0.0,0.0]
    nwat          = 0
    natom             = read_pdb(infile,pdb_data,cryst)
    read_json(infile,mol_list,res_name_j,res_map_j)
    try:
        fout = open(outfile,"w")
    except:
        print ("ERROR: FILE",outfile,"CANNOT OPEN")
        sys.exit(0)
    nres = make_domain(pdb_data,domain_data)
    print ("# of ATOMS:",natom)
    print ("# of RESIDUES:",nres)
    print ()
    print ("**RESIDUES IN AA PDB**")
    for i in range(nres):
        aa_pdb_domain = domain_data[i]
        mass = 0; comx = 0; comy = 0; comz = 0
        TYPE = aa_pdb_domain[RESNAME][0]
        print (TYPE)
        if prot_list.count(TYPE):
            # if 1
            if TYPE != "ALA" :
                for j in GBB:
                    if aa_pdb_domain[ATMNAME].count(j):
                        l = aa_pdb_domain[ATMNAME].index(j)
                        x = aa_pdb_domain[POSX][l]
                        y = aa_pdb_domain[POSY][l]
                        z = aa_pdb_domain[POSZ][l]
                        this_mass = atm_mass[j[0][0]]
                        mass += this_mass
                        comx += x*this_mass
                        comy += y*this_mass
                        comz += z*this_mass
                    else :
                        bad_res = aa_pdb_domain[RESID][0]
                        print ("ERROR:",TYPE,"BACKBONE ATOM",j,"IS MISSING ON RESIDUE",bad_res)
                        sys.exit(0)
                # j
                comx/=mass; comy/=mass; comz/=mass;
# CREATE A NEW DOMAIN FOR THIS RESIDUE BACKBONE
# AND STORE THE COM OF THE BACKBONE ATOMS
                if TYPE != "GLY" :
                    #cg_pdb_domain[ATMNAME] = "GB"
                    cg_pdb_domain[ATMNAME] = "GBB"
                    cg_pdb_domain[POSX] = comx
                    cg_pdb_domain[POSY] = comy
                    cg_pdb_domain[POSZ] = comz
                else:
                    cg_pdb_domain[ATMNAME] = "GB"
                    #cg_pdb_domain[ATMNAME] = "GBM"
                    cg_pdb_domain[POSX] = comx
                    cg_pdb_domain[POSY] = comy
                    cg_pdb_domain[POSZ] = comz

                cg_pdb_domain[RESNAME] = TYPE
                cg_pdb_domain[CHAINID] = aa_pdb_domain[CHAINID][0]
                cg_pdb_domain[RESID]   = aa_pdb_domain[RESID][0]
                all_cg_pdb_domain.append([cg_pdb_domain[m] for m in range(ncomp)])
#############################
# Just finished the backbone
#############################
# NOW TAKE CARE OF THE SIDECHAINS
                if TYPE != "GLY" :
                    nmtg = 1
                    if TYPE in res_map: 
                        for j in res_map[TYPE]:
                            comx = 0; comy = 0; comz = 0; mass = 0
                            MISSATERROR = 0 
# LOOP OVER EACH SIDECHAIN BEAD HERE
                            for k in j:
                                if aa_pdb_domain[ATMNAME].count(k):
                                    l = aa_pdb_domain[ATMNAME].index(k)
                                    x = aa_pdb_domain[POSX][l] 
                                    y = aa_pdb_domain[POSY][l]
                                    z = aa_pdb_domain[POSZ][l]
                                    this_mass = atm_mass[k[0][0]]
                                    mass += this_mass
                                    comx += x*this_mass
                                    comy += y*this_mass
                                    comz += z*this_mass
                                else :
                                    print ("ERROR:",TYPE,"BACKBONE ATOM",k,"IS MISSING ON RESIDUE")
                                    MISSATERROR = 1

                            #k
                            if MISSATERROR == 0:
                                comx/=mass; comy/=mass; comz/=mass;
                                bead_name = TYPE
                            else:
                                bead_name = "UNK"
                                comx = 0; comy = 0; comz = 0
# CHECK TO SEE IF WE HAVE MULTIPLE BEAD AND NEED TO ADD INDEX TO SIDECHAIN BEAD NAME
                            if len(res_map[TYPE]) > 1:
                                bead_name = bead_name.strip()[:2]
                                bead_name = "%s%d" % (bead_name,nmtg)
                                nmtg += 1

                            if MISSATERROR == 0:
                                cg_pdb_domain[ATMNAME] = bead_name 
                                cg_pdb_domain[POSX]    = comx
                                cg_pdb_domain[POSY]    = comy
                                cg_pdb_domain[POSZ]    = comz
                            else:
                                bead_name = "UNK"
                                cg_pdb_domain[ATMNAME] = bead_name
                                cg_pdb_domain[POSX]    = comx
                                cg_pdb_domain[POSY]    = comy
                                cg_pdb_domain[POSZ]    = comz
                            all_cg_pdb_domain.append([cg_pdb_domain[m] for m in range(ncomp)])
                        #j
                    else:
                        print ("WARNING: RESIDUE",TYPE,"IS NOT DEFINED IN RESMAP")
            else: #if 1 else
                comx = 0; comy = 0; comz = 0; mass = 0
                for j in ABB:
                    if aa_pdb_domain[ATMNAME].count(j):
                        l = aa_pdb_domain[ATMNAME].index(j)
                        x = aa_pdb_domain[POSX][l] 
                        y = aa_pdb_domain[POSY][l]
                        z = aa_pdb_domain[POSZ][l]
                        this_mass = atm_mass[j[0]]
                        mass += this_mass
                        comx += x*this_mass
                        comy += y*this_mass
                        comz += z*this_mass
                    else :
                        print ("ERROR: ALANINE ATOM",j,"IS MISSING")
                        sys.exit(0)
                #cg_pdb_domain[ATMNAME] = "AB"
                cg_pdb_domain[ATMNAME] = "ABB"
                cg_pdb_domain[RESNAME] = TYPE
                cg_pdb_domain[CHAINID] = aa_pdb_domain[CHAINID][0]
                cg_pdb_domain[RESID]   = aa_pdb_domain[RESID][0]
                comx/=mass; comy/=mass; comz/=mass;
                cg_pdb_domain[POSX]    = comx
                cg_pdb_domain[POSY]    = comy
                cg_pdb_domain[POSZ]    = comz
                all_cg_pdb_domain.append([cg_pdb_domain[m] for m in range(ncomp)])
                #if 1 end
        elif mol_list.count(TYPE):
            if TYPE in res_map_j: 
                nmtg = 0
                for j in res_map_j[TYPE]:
                    comx = 0; comy = 0; comz = 0; mass = 0
                    if set(j) == set(list(set(j) & set(aa_pdb_domain[ATMNAME]))):
                        iatm = 0
                        if set(j) == set(["OH2","H1","H2"]) or set(j) == set(["OW","HW1","HW2"]):
                            nwat += 1 
                            if (nwat - 1) % 3 != 0: 
                                continue
                                for k in j:
                                    l = aa_pdb_domain[ATMNAME].index(k)
                                    x = aa_pdb_domain[POSX][l] 
                                    y = aa_pdb_domain[POSY][l]
                                    z = aa_pdb_domain[POSZ][l]
                                    this_mass = atm_mass[k[0][0]]
                                    mass += this_mass
                                    comx += x*this_mass
                                    comy += y*this_mass
                                    comz += z*this_mass
                                #k
                    else :
                        print ("ERROR:",TYPE,"ATOM",j,"IS MISSING ON MOLECULE")
                        sys.exit(0)
                    
                    comx/=mass; comy/=mass; comz/=mass;
                    bead_name = res_name_j[TYPE][nmtg]
                    cg_pdb_domain[ATMNAME] = bead_name
                    cg_pdb_domain[RESNAME] = TYPE
                    cg_pdb_domain[CHAINID] = aa_pdb_domain[CHAINID][0]
                    cg_pdb_domain[RESID]   = aa_pdb_domain[RESID][0]
                    cg_pdb_domain[POSX]    = comx
                    cg_pdb_domain[POSY]    = comy
                    cg_pdb_domain[POSZ]    = comz
                    all_cg_pdb_domain.append([cg_pdb_domain[m] for m in range(ncomp)])
                    nmtg += 1
            #j
            else:
                print ("WARNING: RESIDUE",TYPE,"IS NOT DEFINED IN RESMAP")
        else:
            print ("ERROR: RESIDUE",TYPE,"IS NOT IN MOL_LIST")
            sys.exit(0)
    #i
    write_domain(all_cg_pdb_domain,cryst,fout)
    fout.close()

    print ()
    print ("Normal Termination.")
    print ()
