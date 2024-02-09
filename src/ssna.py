#!/usr/bin/env python

import sys
from argparse import ArgumentParser


def get_option():
    argparser = ArgumentParser()
    argparser.add_argument('cgpdb', type=str,
                            help='Specify input CG PDB file name.')
    argparser.add_argument('output', type=str,
                            help='Specify output topology file name.')
    argparser.add_argument('-rna', action='store_true',
                            help='Generate a topology file for RNA (default: for DNA).')
    argparser.add_argument('-pspica', action='store_true',
                            help='Assign partial charge (0.5990) for pSPICA FF (default: 0.1118, for SPICA FF).')
    return argparser.parse_args()


def get_option_script(argv):
    argparser = ArgumentParser(usage='ssna [-h] [-pspica] [-rna] cgpdb output',
                               prog ="ssna")
    argparser.add_argument('cgpdb', type=str,
                            help='Specify input CG PDB file name.')
    argparser.add_argument('output', type=str,
                            help='Specify output topology file name.')
    argparser.add_argument('-rna', action='store_true',
                            help='Generate a topology file for RNA (default: for DNA).')
    argparser.add_argument('-pspica', action='store_true',
                            help='Assign partial charge (0.5990) for pSPICA FF (default: 0.1118, for SPICA FF).')
    return argparser.parse_args(argv)

ncomp = 7
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

###############################################################################

spec_dihed={}
spec_dihed["ADE"] = ["AD1", "AD2", "AD3", "AD4"]
spec_dihed["GUA"] = ["GU1", "GU2", "GU3", "GU4"]
#######################################################

# CHARGE HASH
charge={}
charge["PB"] = -0.1118

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


def read_pdb(cgpdb, pdb_list, cryst, ters):
    natom = 0
    try:
        f = open(cgpdb,"r")
    except:
        print("ERROR: File", cgpdb, "is not found.")
        sys.exit(1)
    line = f.readline()
    while line:
        recname = line[0:6].strip()
        if recname == "CRYST1":
            cryst.append(line.strip())
        elif recname == "ATOM":
            natom += 1
            index = line[6:11]
            atmname = line[12:16]
            indicat = line[16:17]
            resname = line[17:21]
            chainid = line[21:22]
            resid = line[22:26]
            code = line[26:27]
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
        print("ERROR: File",outfile,"cannot open.")
        sys.exit(0)
    return fout


class ssna:
    def __init__(self, cgpdb, outfile, rna, pspica):
        self.cgpdb = cgpdb
        self.outfile = outfile
        self.rna = rna
        self.pspica  = pspica
        self.nat = 0
        self.nbb = 0
        self.bbndx = []
        self.iat2ibb = []
        self.bbtype = []
        self.name = []
        self.resname = []
        self.resid = []
        self.bBackbone = []
        self.bPH1TY1 = []
        self.coord = []
        self.pdb_data = []
        self.ters = []
        cryst = []
        self.natom = read_pdb(cgpdb, self.pdb_data, cryst, self.ters)
        self.ftop = open_file(outfile)
        self._charge_mod()
        self._set_array()


    def _charge_mod(self):
        if self.pspica:
            charge_unit   = 0.5590 
            charge["PB"] = charge_unit


    def _set_array(self):
        for i in range(self.natom):
            tmp_name = self.pdb_data[i][PDB_ATMNAME].strip()
            tmp_resname = self.pdb_data[i][PDB_RESNAME].strip()
            if tmp_name in ["PB", "RB1", "RB2", "DB2"]:
                if tmp_resname in ["ADE", "GUA", "CYT", "THY", "URA"]:
                    self.name.append(tmp_name)
                else:
                    sys.exit("ERROR: Unknown resname", tmp_name)
                self.bbtype.append(tmp_resname)
                self.resname.append(tmp_resname)
                self.bbndx.append(self.nat)
                self.resid.append(int(self.pdb_data[i][PDB_RESID]))
                self.coord.append([self.pdb_data[i][PDB_POSX], self.pdb_data[i][PDB_POSY], self.pdb_data[i][PDB_POSZ]])
                self.iat2ibb.append(self.nbb)
                self.bBackbone.append(1)
                self.bPH1TY1.append(0)
                self.nat += 1; self.nbb += 1
                if tmp_name in ["RB2", "DB2"]:
                    for j in res_map[tmp_resname]:
                        self.name.append(j)
                        self.coord.append([self.pdb_data[self.nat][PDB_POSX], self.pdb_data[self.nat][PDB_POSY], self.pdb_data[self.nat][PDB_POSZ]])
                        self.resname.append(tmp_resname)
                        self.resid.append(int(self.pdb_data[i][PDB_RESID]))
                        self.iat2ibb.append(None)
                        self.bBackbone.append(0)
                        if j in ["AD3", "GU3"]:
                            self.bPH1TY1.append(0)
                        else :
                            self.bPH1TY1.append(1)
                        self.nat += 1
        tmp_resid = self.resid[0]
        j = 0
        for i in range(len(self.resid)):
            if self.resid[i] != tmp_resid:
                tmp_resid = self.resid[i]
                j += 1
            self.resid[i] = j


    def write_atom(self):
        ftop = self.ftop
        resid = self.resid
        for i in range(self.nat):
            name_i = self.name[i]
            resname_i = self.resname[i]
            I = i + 1
            if name_i in charge:
                thischrg = charge[name_i]
            else :
                thischrg = 0.0
            wtype = this_type[name_i]
            mass = bd_mass[name_i]
            print ("atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
                %(I,resname_i,name_i,wtype,mass,thischrg), file=ftop)
        print(file=ftop)


    def write_bond_BB(self):
        ftop = self.ftop
        self.total_bonds = 0
        bond_index1 = []
        bond_index2 = []
        for i in range(self.nbb - 1):
            if self.resid[self.bbndx[i]] > self.resid[self.bbndx[i+1]]:
                continue
            if self.resid[self.bbndx[i]] in self.ters:
                continue
            self.total_bonds += 1
            print("bond %5d %5d # %s-%s" \
                %(self.bbndx[i]+1, self.bbndx[i+1]+1, self.name[self.bbndx[i]], self.name[self.bbndx[i+1]]), file=ftop)
            bond_index1.append(self.bbndx[i])
            bond_index2.append(self.bbndx[i+1])
        print(file=ftop)
        self.bond_index1 = bond_index1
        self.bond_index2 = bond_index2


    def write_bond_SC(self):
        ftop = self.ftop
        for i in range(self.nbb):
            bbtype_i = self.bbtype[i]
            bbndx_i = self.bbndx[i]
            if self.name[bbndx_i] in ["PB", "RB1"]:
                continue
            if bbtype_i in spec_bond:
                for j in spec_bond[bbtype_i]:
                    bndx1 = 0
                    bndx2 = 0
                    for k in j:
                        if k == "DB2" and self.rna:
                            k = "RB2"
                        if bndx1 == 0:
                            tmp_ndx = bbndx_i
                            while self.name[tmp_ndx] != k:
                                tmp_ndx += 1
                            bndx1 = tmp_ndx
                        else:
                            tmp_ndx = bbndx_i
                            while self.name[tmp_ndx] != k:
                                tmp_ndx += 1
                            bndx2 = tmp_ndx
                    print("bond %5d %5d # %s-%s" \
                        %(bndx1+1, bndx2+1, self.name[bndx1], self.name[bndx2]), file=ftop)
                    self.total_bonds += 1
                    self.bond_index1.append(bndx1)
                    self.bond_index2.append(bndx2)
            else :
                print("bond %5d %5d # %s-%s" \
                    %(bbndx_i+1, bbndx_i+1+1, self.name[bbndx_i], self.name[bbndx_i+1]), file=ftop)
                self.total_bonds += 1
                self.bond_index1.append(bbndx_i)
                self.bond_index2.append(bbndx_i+1)
        print(file=ftop)


    def _write_angle(self,andx1,andx2,andx3):
        ftop = self.ftop
        name = self.name
        bBackbone = self.bBackbone
        bPH1TY1 = self.bPH1TY1
        resid = self.resid
        if bBackbone[andx1] + bBackbone[andx2] + bBackbone[andx3] != 0 or bPH1TY1[andx2] == 1:
            names = [name[andx1], name[andx2], name[andx3]]
            if name[andx2] in ["CY1", "CY3", "TH1", "TH3", "UR1", "UR3"]:
                return
            if set(names) == set(res_map["CYT"]) or set(names) == set(res_map["THY"]) or set(names) == set(res_map["URA"]):
                return
            if name[andx2] in ["AD1","GU1"]:
                print("angleparam %5d %5d %5d  0.0  90.0000 # %s %s %s"\
                    %(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3]), file=ftop)
            else:
                print("angle      %5d %5d %5d               # %s %s %s" \
                    %(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3]), file=ftop)
                self.angle_index.append([andx1+1, andx2+1, andx3+1])
                self.angle_names.append(names)


    def write_angle(self):
        ftop = self.ftop
        bond_index1 = self.bond_index1
        bond_index2 = self.bond_index2
        self.angle_index = []
        self.angle_names = []
        for i1 in range(self.total_bonds):
            for i2 in range(i1 + 1, self.total_bonds):
                if bond_index1[i1] == bond_index1[i2]:
                    self._write_angle(bond_index2[i1],bond_index1[i1],bond_index2[i2])
                elif bond_index1[i1] == bond_index2[i2]:
                    self._write_angle(bond_index2[i1],bond_index1[i1],bond_index1[i2])
                elif bond_index2[i1] == bond_index1[i2]:
                    self._write_angle(bond_index1[i1],bond_index2[i1],bond_index2[i2])
                elif bond_index2[i1] == bond_index2[i2]:
                    self._write_angle(bond_index1[i1],bond_index2[i1],bond_index1[i2])
        print(file=ftop)


    def write_dihedral(self):
        ftop = self.ftop
        nang = len(self.angle_index)
        improper_index = []
        improper_names = []
        for j in range(nang-1):
            for k in range(j+1, nang):
                ang_j = self.angle_index[j]
                ang_k = self.angle_index[k]
                if ang_j[1:] == ang_k[:2]:
                    ndx1 = ang_j[0]
                    ndx2 = ang_j[1]
                    ndx3 = ang_j[2]
                    ndx4 = ang_k[2]
                elif ang_j[1::] == list(reversed(ang_k[1:])):
                    ndx1 = ang_j[0]
                    ndx2 = ang_j[1]
                    ndx3 = ang_j[2]
                    ndx4 = ang_k[0]
                elif list(reversed(ang_j[:2])) == ang_k[:2]:
                    ndx1 = ang_j[2]
                    ndx2 = ang_j[1]
                    ndx3 = ang_j[0]
                    ndx4 = ang_k[2]
                elif ang_j[:2] == ang_k[1:]:
                    ndx1 = ang_j[2]
                    ndx2 = ang_j[1]
                    ndx3 = ang_j[0]
                    ndx4 = ang_k[0]
                else:
                    continue
                name1 = self.name[ndx1-1]
                name2 = self.name[ndx2-1]
                name3 = self.name[ndx3-1]
                name4 = self.name[ndx4-1]
                dihed_names = [name1, name2, name3, name4]
                if "RB1" in dihed_names:
                    if ["RB1", "PB", "DB2"] == dihed_names[1:] or ["RB1", "PB", "DB2"] == dihed_names[:3]:
                        continue
                    if ["RB1", "PB", "RB2"] == dihed_names[1:] or ["RB1", "PB", "RB2"] == dihed_names[:3]:
                        continue
                    if dihed_names.count("RB1") == 2 or dihed_names.count("DB2") == 2 or dihed_names.count("RB2") == 2:
                        print("dihedral  %5d %5d %5d %5d # %s-%s-%s-%s" \
                            %(ndx1, ndx2, ndx3, ndx4, name1, name2, name3, name4), file=ftop)
                    else:
                        improper_index.append([ndx1, ndx2, ndx3, ndx4])
                        improper_names.append(dihed_names)
        for i in range(self.nbb):
            bbtype_i = self.bbtype[i]
            bbndx_i = self.bbndx[i]
            name_i = self.name[bbndx_i]
            if bbtype_i in spec_dihed and name_i not in ["PB", "RB1"]:
                ndx1 = bbndx_i + 1
                ndx2 = ndx1 + 1
                ndx3 = ndx2 + 1
                ndx4 = ndx3 + 1
                print("dihedral  %5d %5d %5d %5d # %s-%s-%s-%s" \
                    %(ndx1+1, ndx2+1, ndx3+1, ndx4+1,
                      self.name[ndx1],self.name[ndx2],self.name[ndx3],self.name[ndx4]), file=ftop)
        print(file=ftop)
        for idx, name in zip(improper_index, improper_names):
            print ("improper  %5d %5d %5d %5d # %s-%s-%s-%s" \
                %(idx[0], idx[1], idx[2], idx[3],
                  name[0], name[1], name[2], name[3]), file=ftop)
        print(file=ftop)


    def finalize(self):
        self.ftop.close()
        if self.pspica:
            print('## ssDNA/RNA topology file has been generated for pSPICA FF ##')
            print('NOTE: Remove "-pspica" option when you use generated protein topology files for SPICA FF.')
        else:
            print('## ssDNA/RNA topology file has been generated for SPICA FF ##')
            print('NOTE: Use "-pspica" option when you use generated protein topology files for pSPICA FF.')
        print()


    def run(self):
        self.write_atom()
        self.write_bond_BB()
        self.write_bond_SC()
        self.write_angle()
        self.write_dihedral()
        self.finalize()
     

if __name__ == "__main__":
    args = get_option()
    cgpdb = args.cgpdb
    outfile = args.output
    rna = args.rna
    pspica = args.pspica
    obj = ssna(cgpdb, outfile, rna, pspica)
    obj.run()
