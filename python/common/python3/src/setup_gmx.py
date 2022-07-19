import numpy as np
import sys, textwrap
from math import sqrt
from math import pi
from pathlib import Path
from argparse import ArgumentParser

def get_angle(r1, r2, r3):
    r12   = r1 - r2
    r32   = r3 - r2
    r13_inn = np.dot(r12, r32)
    r13_mag = np.linalg.norm(r12)*np.linalg.norm(r32)
    cos13   = r13_inn/r13_mag
    if cos13 < -1.0:
       return np.pi
    elif cos13 > 1.0:
       return 0.0
    else:
       return np.arccos(cos13)

def get_dihedral(r1, r2, r3, r4):
    r12   = r1 - r2
    r32   = r3 - r2
    r43   = r4 - r3
    r23   = r2 - r3
    r13_inn = np.dot(r12, r32)
    r42_inn = np.dot(r43, r23)
    r13_mag = np.linalg.norm(r12)*np.linalg.norm(r32)
    r42_mag = np.linalg.norm(r43)*np.linalg.norm(r23)
    cos13   = r13_inn/r13_mag
    cos42   = r42_inn/r42_mag
    p1 = r12 - cos13*r32/np.linalg.norm(r32)*np.linalg.norm(r12)
    p4 = r43 - cos42*r23/np.linalg.norm(r23)*np.linalg.norm(r43)
    cosp    = np.dot(p1, p4)/(np.linalg.norm(p1)*np.linalg.norm(p4))
    p14_crs = np.cross(p1, p4)
    if np.dot(p14_crs, r32) > 0.0:
        sign = 1
    else:
        sign = -1
    if cosp < -1.0:
        return sign*np.pi
    elif cosp > 1.0:
        return 0.0
    else:
        return sign*np.arccos(cosp)
    
class Sysdat:
    nats = nbnds = nangs = nimprops = ndiheds = ntops = 0 
    total_ats = total_bnds = total_angs = total_improps = total_diheds = 0
    foundatoms = boxinfo = 0
    uniq_nats = uniq_nbnds = uniq_nangs = uniq_nimprops = uniq_ndiheds = 0
    param_bnds, param_angs  = [], []
    coordx, coordy, coordz = [], [], []
    boxx = boxy = boxz = 0.0

class Topdat:
    def __init__(self):
        self.nat = self.nbnd = self.nang = self.nimprop = self.nmol = 0
        self.bndndx1, self.bndndx2, self.bndtype = [], [], []
        self.angndx1, self.angndx2, self.angndx3, self.angtype = [], [], [], []
        self.improp_func, self.impropndx1, self.impropndx2, self.impropndx3, self.impropndx4, self.improptype   = [], [], [], [], [], []
        self.dihed_func, self.dihedndx1, self.dihedndx2, self.dihedndx3, self.dihedndx4, self.dihedtype, self.dihedn = [], [], [], [], [], [], []
        self.dihedpset, self.improppset, self.bndpset, self.angpset = [], [], [], []
        self.ind, self.parm_atomtype, self.ndihed, self.dihedeq   = [], [], [], []
        self.dihedfk, self.dihedof = [], []
        self.mass, self.charge, self.bndfk, self.bndeq, self.angfk, self.angeq, self.impropfk, self.impropeq = [], [], [], [], [], [], [], []
        self.atomname, self.atomtype, self.segid, self.resname = [], [], [], []
   
class Database:
     fbnd, bnde, fang, ange, eps, sig, angsdk = [], [], [], [], [], [], []
     nvdwtype, nbndtype, nangtype = [], [], []
     vdwtype1, vdwtype2, vdwstyle = [], [], []
     bndtype1, bndtype2 = [], []
     angtype1, angtype2, angtype3 = [], [], []

def read_pdb(fname, sysdat):
    col = 30
    with open(fname, "r") as fin:
        line = fin.readline()
        while line:
            items = line.split()
            if len(items) == 0:
                line = fin.readline()
                continue
            if items[0] == "CRYST1":
                print("FOUND BOXSIZE DATA.")
                sysdat.boxx = float(items[1])
                sysdat.boxy = float(items[2])
                sysdat.boxz = float(items[3])
            elif items[0] == "ATOM" or items[0] == "HETATM":
                if sysdat.foundatoms >= sysdat.total_ats:
                    sys.exit("ERROR: found atoms in pdb file >= total atoms in top files.")
                sysdat.coordx.append(float(line[col:col+8]))
                sysdat.coordy.append(float(line[col+8:col+16]))
                sysdat.coordz.append(float(line[col+16:col+24]))
                sysdat.foundatoms += 1
            line = fin.readline()
        if sysdat.foundatoms == 0:
            sys.exit("ERROR: DID NOT FIND ANY ATOMS IN THE PDB FILE.")
        if sysdat.boxx == 0.0:
            print("WARNING: DID NOT FIND CELL SIZE.")
            print("BOX SIZE WILL HAVE BE SET BY HAND.")

def read_coords(database, topdat, sysdat):
    with open("molecule.itp", "w") as fout:
        print("; generated by setup_gmx", file=fout)
        for idx in range(sysdat.ntops):
            print(file=fout)
            print("[ moleculetype ]", file=fout)
            print("; name nrexcl", file=fout)
            print("{}      2".format(topdat[idx].resname[0]), file=fout)
            print(file=fout)
            print("[ atoms ]", file=fout)
            print("; nr    type    resnr    residu   atom   cgnr   charge  mass", file=fout);

            for jdx in range(topdat[idx]nat):
                print("{:6d} {:6d}    1    {:6s} {:6s} {:6d} {:8.4f}".format(jdx+1,topdat[idx].type[jdx],
                        topdat[idx].resname[jdx],topdat[idx].name[jdx],jdx+1,topdat[idx].charge[jdx]), file=fout)
            if topdat[idx].nbnd > 0:
                print(file=fout)
                print("[ bonds ]", file=fout)
                print(file=fout)
                for jdx in range(topdat[idx].nbnd):
                    if topdat[idx].bndpset[jdx] == True:
                        tmpcalca = topdat[idx].bndeq[jdx]/10.0
                        tmpcalcb= topdat[idx].bndfk[jdx]*4.184*2.0*100
                        print("{:5d} {:5d}    6  {:8.4f}  {:8.4f} ; EN specified by bondparam".format(topdat[idx].bndndx1[jdx],
                                topdat[idx].bndndx2[jdx],tmpcalca,tmpcalcb), file=fout)
                    else :
                        print("{:5d} {:5d}    1".format(topdat[idx].bndndx1[jdx],topdat[idx].bndndx2[jdx]), file=fout)
            if topdat[idx].nang > 0:
                print(file=fout)
                print("[ angles ]",file=fout)
                print(file=fout)
                for jdx in range(topdat[idx].nang):
                    if topdat[idx].angpset[jdx] == False:
                        # now compare to the database
                        for kdx in range(database.nangtype):
                            if cmp_wc(database.angtype2[kdx], topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1]):
                                f1 = cmp_wc(database.angtype1[kdx], topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1])
                                f2 = cmp_wc(database.angtype3[kdx], topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1])
                                if f1 and f2:
                                    datndx = kdx
                                    break
                                f1 = cmp_wc(database.angtype3[kdx], topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1])
                                f2 = cmp_wc(database.angtype1[kdx], topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1])
                                if f1 and f2:
                                    datndx = kdx
                                    break
                        ifound = 0
                        for kdx in range(database.nvdwtype):
                            f1 = database.vdwtype1[kdx] == topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1]
                            f2 = database.vdwtype2[kdx] == topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1]
                            f3 = database.vdwtype1[kdx] == topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1]
                            f4 = database.vdwtype2[kdx] == topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1]
                            if f1 and f2:
                                ifound = 1
                                vdwtmp = kdx
                                break
                            elif f3 and f4:
                                ifound = 1
                                vdwtmp = kdx
                                break
                        eps = database.eps[vdwtmp]*4.184
                        sig = database.sig/10.0
                        tmpcalca = database.ange[datndx]
                        tmpcalcb = fang[datndx]*4.184*2.0
                        if tmpcalca < 0:
                            i1 = topdat[idx].angndx1[jdx]-1
                            i2 = topdat[idx].angndx2[jdx]-1
                            i3 = topdat[idx].angndx3[jdx]-1
                            r1 = np.array([sysdat.coordx[i1],sysdat.coordy[i1],sysdat.coordz[i1]])
                            r2 = np.array([sysdat.coordx[i2],sysdat.coordy[i2],sysdat.coordz[i2]])
                            r3 = np.array([sysdat.coordx[i3],sysdat.coordy[i3],sysdat.coordz[i3]])
                            angle_in_pdb = 180.0/pi*get_angle(r1,r2,r3)
                            print("{:5d} {:5d} {:5d}    1  {:8.4f} {:8.4f} ; {:6s} {:6s} {:6s} taken from pdb"
                                    .format(topdat[idx].angndx1[jdx],topdat[idx].angndx2[jdx],.topdat[idx].angndx3[jdx],
                                            angle_in_pdb,tmpcalcb,database.angtype1[datndx],database.angtype2[datndx],database.angtype3[datndx]), file=fout)
                        else:
                            if database.angsdk[datndx]:
                                print("{:5d} {:5d} {:5d}    1".format(topdat[idx].angndx1[jdx],topdat[idx].angndx2[jdx],topdat[idx].angndx3[jdx]), file=fout)
                            else:
                                print("{:5d} {:5d} {:5d}    5".format(topdat[idx].angndx1[jdx],topdat[idx].angndx2[jdx],topdat[idx].angndx3[jdx]), file=fout)
            if topdat[idx].ndihed > 0:
                print(file=fout)
                print("[ dihedrals ]",file=fout)
                print(file=fout)
                for jdx in range(topdat[idx].ndihed):
                    print("{:5d} {:5d} {:5d} {:5d}  1  {:<3d} {:8.4f} {:<3d} ; FROM TOP".format(
                                    topdat[idx].dihedndx1[jdx],topdat[idx].dihedndx2[jdx],topdat[idx].dihedndx3[jdx],topdat[idx].dihedndx4[jdx],
                                    topdat[idx].dihedeq[jdx],topdat[idx].dihedfk[jdx]*4.184,topdat[idx].dihedn[jdx]), file=fout)

            if topdat[idx].nimprop > 0:
                print(file=fout)
                print("[ impropers ]", file=fout)
                print(file=fout)
                for jdx in range(topdat[idx].nimprop)
                    print("{:5d} {:5d} {:5d} {:5d}  2  {:8.4f} {:8.4f} ; FROM TOP".format(
                            topdat[idx].impropndx1[jdx],topdat[idx].impropndx2[jdx],topdat[idx].impropndx3[jdx],topdat[idx].impropndx4[jdx],
                            topdat[idx].impropfk[jdx]*4.184*2.0,topdat[idx].impropeq[jdx]))
