import numpy as np
import sys, re, os
from pathlib import Path
from math import sqrt
from math import pi
from argparse import ArgumentParser

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument('input_files', type=str, nargs="+",
                            help='<topfile 1> <nmol 1> [ <topfile 2> <nmol 2> ..... <topfile n> <nmol n>] <param file> [ <coordfile> ]')
    argparser.add_argument('-prot', action="store_true",
                            help='read a pdb file to extract reference angle for protein models.')
    return argparser.parse_args()

def get_option_script(argv):
    argparser = ArgumentParser(usage='setup_gmx [-h] [-prot] input_files',
                               prog ="setup_gmx")
    argparser.add_argument('input_files', type=str, nargs="+",
                            help='<topfile 1> <nmol 1> [ <topfile 2> <nmol 2> ..... <topfile n> <nmol n>] <param file> [ <coordfile> ]')
    argparser.add_argument('-prot', action="store_true",
                            help='read a pdb file to extract reference angle for protein models.')
    return argparser.parse_args(argv)

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
    foundatoms = boxinfo = ischarged = 0
    uniq_nats = uniq_nbnds = uniq_nangs = uniq_nimprops = uniq_ndiheds = 0
    param_bnds, param_angs  = [], []
    coordx, coordy, coordz = [], [], []
    boxx = boxy = boxz = 0.0

class Topdat:
    def __init__(self):
        self.fname = None
        self.nat = self.nbnd = self.nang = self.nimprop = self.nmol = self.ngo = 0
        self.gondx1, self.gondx2, self.gofunctype, self.eps, self.sig = [], [], [], [], []
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
                print("Found boxsize data.")
                sysdat.boxx = float(items[1])
                sysdat.boxy = float(items[2])
                sysdat.boxz = float(items[3])
            elif items[0] == "ATOM" or items[0] == "HETATM":
                if sysdat.foundatoms >= sysdat.total_ats:
                    sys.exit("ERROR: Found atoms in pdb file >= total atoms in top files.")
                sysdat.coordx.append(float(line[col:col+8]))
                sysdat.coordy.append(float(line[col+8:col+16]))
                sysdat.coordz.append(float(line[col+16:col+24]))
                sysdat.foundatoms += 1
            line = fin.readline()
        if sysdat.foundatoms == 0:
            sys.exit("ERROR: Did not find any atoms in the pdb file.")
        if sysdat.boxx == 0.0:
            print("WARNING: Did not find cell size.")
            print("Box size will have be set by hand.")

def read_coords(database, topdat, sysdat):
    os.makedirs("toppar", exist_ok=True)
    for idx in range(sysdat.ntops):
        with open(f"toppar/{topdat[idx].fname}.itp", "w") as fout:
            print("; generated by cg_spica setup_gmx", file=fout)
            print(file=fout)
            print("[ moleculetype ]", file=fout)
            print("; name nrexcl", file=fout)
            if topdat[idx].resname[0] in ["PWAT", "PSOD", "PCLA"]: 
                print("{}      1".format(topdat[idx].fname), file=fout)
            else:
                print("{}      2".format(topdat[idx].fname), file=fout)
            print(file=fout)
            print("[ atoms ]", file=fout)
            print("; nr    type    resnr    residu   atom   cgnr   charge  mass", file=fout);

            for jdx in range(topdat[idx].nat):
                print("{:6d} {:>6s}    1    {:>6s} {:>6s} {:6d} {:8.4f}".format(jdx+1,topdat[idx].atomtype[jdx],
                        topdat[idx].resname[jdx],topdat[idx].atomname[jdx],jdx+1,topdat[idx].charge[jdx]), file=fout)
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
                if topdat[idx].resname[0] in ["PWAT", "PSOD", "PCLA"]:
                    print(file=fout)
                    print("[ constraints ]", file=fout)
                    print(file=fout)
                    print("{:5d} {:5d}    1  0.11".format(topdat[idx].bndndx1[0],topdat[idx].bndndx2[0]), file=fout)
            # Go model for protein backbones
            if topdat[idx].ngo > 0:
                print(file=fout)
                print("[ pairs ]", file=fout)
                print("; Go model for protein backbone", file=fout)
                print("; ai   aj    funct    c6    c12", file=fout)
                for jdx in range(topdat[idx].ngo):
                    topdat[idx].sig[jdx] /= 10
                    print("{:5d} {:5d}    1  {:14.5e} {:14.5e};  {} {:5.6f} {:5.6f}".format(topdat[idx].gondx1[jdx],topdat[idx].gondx2[jdx],
                            4.184*4.0*topdat[idx].eps[jdx]*pow(topdat[idx].sig[jdx],6),4.184*4.0*topdat[idx].eps[jdx]*pow(topdat[idx].sig[jdx],12),
                            topdat[idx].gofunctype[jdx],topdat[idx].eps[jdx],topdat[idx].sig[jdx]), file=fout)
                print(file=fout)

            if topdat[idx].nang > 0:
                print(file=fout)
                print("[ angles ]",file=fout)
                print(file=fout)
                for jdx in range(topdat[idx].nang):
                    if topdat[idx].angpset[jdx] == True:
                        if topdat[idx].angfk[jdx] == -1:
                            for kdx in range(database.nangtype):
                                if cmp_wc(database.angtype2[kdx], topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1]):
                                    f1 = cmp_wc(database.angtype1[kdx], topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1])
                                    f2 = cmp_wc(database.angtype3[kdx], topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1])
                                    if f1 and f2:
                                        topdat[idx].angfk[jdx] = database.fang[kdx]
                                        break
                                    f1 = cmp_wc(database.angtype3[kdx], topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1])
                                    f2 = cmp_wc(database.angtype1[kdx], topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1])
                                    if f1 and f2:
                                        topdat[idx].angfk[jdx] = database.fang[kdx]
                                        break
                        tmpcalca = topdat[idx].angeq[jdx]
                        tmpcalcb = topdat[idx].angfk[jdx]*4.184*2.0
                        print("{:5d} {:5d} {:5d}    1  {:8.4f} {:8.4f} ; {:>6s} {:>6s} {:>6s}".format(
                                topdat[idx].angndx1[jdx], topdat[idx].angndx2[jdx], topdat[idx].angndx3[jdx],tmpcalca, tmpcalcb,
                                topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1],topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1],
                                topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1]), file=fout)
                    else:
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
                        sig = database.sig[vdwtmp]/10.0
                        tmpcalca = database.ange[datndx]
                        tmpcalcb = database.fang[datndx]*4.184*2.0
                        if tmpcalca < 0:
                            i1 = topdat[idx].angndx1[jdx]-1
                            i2 = topdat[idx].angndx2[jdx]-1
                            i3 = topdat[idx].angndx3[jdx]-1
                            r1 = np.array([sysdat.coordx[i1],sysdat.coordy[i1],sysdat.coordz[i1]])
                            r2 = np.array([sysdat.coordx[i2],sysdat.coordy[i2],sysdat.coordz[i2]])
                            r3 = np.array([sysdat.coordx[i3],sysdat.coordy[i3],sysdat.coordz[i3]])
                            angle_in_pdb = 180.0/pi*get_angle(r1,r2,r3)
                            print("{:5d} {:5d} {:5d}    1  {:8.4f} {:8.4f} ; {:>6s} {:>6s} {:>6s} taken from pdb"
                                    .format(topdat[idx].angndx1[jdx],topdat[idx].angndx2[jdx],topdat[idx].angndx3[jdx],
                                            angle_in_pdb,tmpcalcb,database.angtype1[datndx],database.angtype2[datndx],database.angtype3[datndx]), file=fout)
                        else:
                            if database.angsdk[datndx] == 0:
                                print("{:5d} {:5d} {:5d}    1".format(topdat[idx].angndx1[jdx],topdat[idx].angndx2[jdx],topdat[idx].angndx3[jdx]), file=fout)
                            else:
                                print("{:5d} {:5d} {:5d}    4".format(topdat[idx].angndx1[jdx],topdat[idx].angndx2[jdx],topdat[idx].angndx3[jdx]), file=fout)
            if topdat[idx].ndihed > 0:
                print(file=fout)
                print("[ dihedrals ]",file=fout)
                print(file=fout)
                for jdx in range(topdat[idx].ndihed):
                    print("{:5d} {:5d} {:5d} {:5d}  1  {:<3f} {:8.4f} {:<3d} ; FROM TOP".format(
                                    topdat[idx].dihedndx1[jdx],topdat[idx].dihedndx2[jdx],topdat[idx].dihedndx3[jdx],topdat[idx].dihedndx4[jdx],
                                    topdat[idx].dihedeq[jdx],topdat[idx].dihedfk[jdx]*4.184,topdat[idx].dihedn[jdx]), file=fout)

            if topdat[idx].nimprop > 0:
                print(file=fout)
                print("[ impropers ]", file=fout)
                print(file=fout)
                for jdx in range(topdat[idx].nimprop):
                    print("{:5d} {:5d} {:5d} {:5d}  2  {:8.4f} {:8.4f} ; FROM TOP".format(
                            topdat[idx].impropndx1[jdx],topdat[idx].impropndx2[jdx],topdat[idx].impropndx3[jdx],topdat[idx].impropndx4[jdx],
                            topdat[idx].impropfk[jdx]*4.184*2.0,topdat[idx].impropeq[jdx]), file=fout)

            # Remove non-native vdw interaction between protein backbone beads forming native contact
            if True in topdat[idx].bndpset or topdat[idx].ngo > 0:
                print(file=fout)
                print("[ exclusions ]", file=fout)
                if topdat[idx].ngo > 0:
                    print(file=fout)
                    for jdx in range(topdat[idx].ngo):
                        if jdx == 0:
                            print("{:4d} {:4d} ".format(topdat[idx].gondx1[jdx],topdat[idx].gondx2[jdx]), end='', file=fout)
                        elif topdat[idx].gondx1[jdx] != topdat[idx].gondx1[jdx-1]:
                            print(file=fout)
                            print("{:4d} {:4d} ".format(topdat[idx].gondx1[jdx],topdat[idx].gondx2[jdx]), end='', file=fout)
                        else:
                            print("{:4d} ".format(topdat[idx].gondx2[jdx]), end='', file=fout)
                if True in topdat[idx].bndpset:
                    tmp_bndndx1 = 0
                    for jdx in range(topdat[idx].nbnd):
                        if topdat[idx].bndpset[jdx] == True:
                            if topdat[idx].bndndx1[jdx] != tmp_bndndx1:
                                tmp_bndndx1 = topdat[idx].bndndx1[jdx]
                                print(file=fout)
                                print("{:4d} {:4d} ".format(topdat[idx].bndndx1[jdx],topdat[idx].bndndx2[jdx]), end='', file=fout)
                            else:
                                print("{:4d} ".format(topdat[idx].bndndx2[jdx]), end='', file=fout)
            print(file=fout)

def write_psf(database, topdat, sysdat):
    with open("out.psf", "w") as fout:
        print("PSF ", file=fout)
        print(file=fout)
        print("       2 !NTITLE", file=fout)
        print("* created by setup_lammps", file=fout)
        print("* dummy", file=fout)
        print(file=fout)
        print("{:8} !NATOM".format(sysdat.total_ats), file=fout)
        atidx = molidx = 0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nmol):
                molidx += 1
                for kdx in range(topdat[idx].nat):
                    atidx += 1
                    print("{:8} {:<4}{:5} {:<4} {:<4} {:<4}  {:9.6f}  {:12.4f}".format(atidx, topdat[idx].resname[kdx], min(9999,molidx), 
                           topdat[idx].resname[kdx], topdat[idx].atomname[kdx],topdat[idx].atomtype[kdx], 
                           topdat[idx].charge[kdx], topdat[idx].mass[kdx]), file=fout)
        print(file=fout)
        print("{:8} !NBOND: bonds".format(sysdat.total_bnds), file=fout)
        bondidx = offset  = 0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nmol):
                for kdx in range(topdat[idx].nbnd):
                    bondidx += 1
                    print("{:>8}{:>8}".format(topdat[idx].bndndx1[kdx]+(jdx*topdat[idx].nat)+offset,
                                              topdat[idx].bndndx2[kdx]+(jdx*topdat[idx].nat)+offset),
                           file=fout, end="")
                    if bondidx % 4 == 0:
                        print(file=fout)
            offset += topdat[idx].nmol*topdat[idx].nat
        print(file=fout)
        print(file=fout)
        print("{:8} !NTHETA: angles".format(sysdat.total_angs),       file=fout)
        angleidx = offset = 0;
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nmol):
                for kdx in range(topdat[idx].nang):
                    angleidx += 1
                    print("{:>8}{:>8}{:>8}".format(topdat[idx].angndx1[kdx]+(jdx*topdat[idx].nat)+offset,
                                                   topdat[idx].angndx2[kdx]+(jdx*topdat[idx].nat)+offset,
                                                   topdat[idx].angndx3[kdx]+(jdx*topdat[idx].nat)+offset),
                           file=fout, end="")
                    if angleidx % 3 ==0:
                        print(file=fout)
            offset += topdat[idx].nmol*topdat[idx].nat
        print(file=fout)
        print(file=fout)
        #print("{:<7} !NPHI: dihedrals".format(sysdat.total_diheds),    file=fout)
        #print(file=fout)
        #print("{:<7} !NIMPHI: impropers".format(sysdat.total_improps), file=fout)
        #print(file=fout)
        #print("       0 !NDON: donors", file=fout)
        #print(file=fout)
        #print("       0 !NACC: acceptors", file=fout)
        #print(file=fout)

def cmp_wc(s1, s2):
    if s1[-1] == "*":
        idx = s1.find("*")
        return s1[:idx] == s2[:idx]
    elif s2[-1] == "*":
        idx = s2.find("*")
        return s1[:idx] == s2[:idx]       
    else:
        return s1 == s2

def get_unique(database, topdat, sysdat):
    uniq_nats = uniq_bnds = uniq_angs = uniq_improps = uniq_diheds = 0
    uniq_atype, uniq_mass, uniq_charge = [], [], []
    bnd_params, bnd_name1, bnd_name2 = [], [], []
    ang_params, ang_vdw = [], []
    os.makedirs("toppar", exist_ok=True)
    with open("toppar/SPICA.itp","w") as fout:
        print("; Generated by setup_gmx", file=fout)
        print(file=fout)

        # first gather the unique atom types
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nat):
                ikeep = 1
                for kdx in range(uniq_nats):
                    if topdat[idx].atomtype[jdx] == uniq_atype[kdx]:
                        ikeep = 0
                        topdat[idx].parm_atomtype.append(kdx)
                        break
                if ikeep == 1:
                    uniq_atype.append(topdat[idx].atomtype[jdx])
                    uniq_mass.append(topdat[idx].mass[jdx])
                    uniq_charge.append(topdat[idx].charge[jdx])
                    topdat[idx].parm_atomtype.append(uniq_nats)
                    uniq_nats += 1
        sysdat.uniq_nats = uniq_nats
        print("[ defaults ]", file=fout);
        print("; nbfun    comb-rule    gen-pairs    fudgeLJ fudgeQQ", file=fout);
        print("1           1             no", file=fout);
        print(file=fout);
        print("[ atomtypes ]", file=fout);
        print("; name   mass   charge   ptype   sigma   epsilon", file=fout);
        for idx in range(uniq_nats):
            print("{:>8s} {:8.4f} {:8.4f}    A    0.0    0.0".format(uniq_atype[idx], uniq_mass[idx], uniq_charge[idx]), file=fout)
        print(file=fout)

        # get pair interactions
        print("[ nonbond_params ]", file=fout);
        print("; i     j    func   C    A", file=fout);
        for idx in range(uniq_nats):
            for jdx in range(idx, uniq_nats):
                ifound=0;
                for kdx in range(database.nvdwtype):
                    if database.vdwtype1[kdx] == uniq_atype[idx] and database.vdwtype2[kdx] == uniq_atype[jdx]:
                        ifound = 1
                        vdwtmp = kdx
                        break
                    elif database.vdwtype2[kdx] == uniq_atype[idx] and database.vdwtype1[kdx] == uniq_atype[jdx]:
                        ifound = 1
                        vdwtmp = kdx
                        break
                if ifound == 0:
                    print("*********************")
                    print("WARNING: No params for VDW interaction between {} and {}".format(uniq_atype[idx],uniq_atype[jdx]))
                    print("UPDATE DATABASE!!!")
                    sys.exit(1)
                elif ifound == 1:
                    eps = database.eps[vdwtmp]
                    sig = database.sig[vdwtmp]/10.0
                    sigsq = sig*sig
                    sigcub = sig*sig*sig
                    if database.vdwstyle[vdwtmp] == "lj12_4":
	                    pf124 = 4.184*3.0*1.73205080757/2.0
	                    disp = pf124*eps*sigsq*sigsq
	                    repul = pf124*eps*sigcub*sigcub*sigcub*sigcub
                    elif database.vdwstyle[vdwtmp] == "lj9_6":
	                    pf96 = 4.184*27.0/4.0
	                    disp = pf96*eps*sigcub*sigcub
	                    repul = pf96*eps*sigcub*sigcub*sigcub
                    elif database.vdwstyle[vdwtmp] == "lj12_5":
	                    pf125 = 4.184*3.203779841
	                    disp = pf125*eps*sigcub*sigsq
	                    repul = pf125*eps*sigcub*sigcub*sigcub*sigcub
                    elif database.vdwstyle[vdwtmp] == "lj12_6":
	                    pf126 = 4.184*4.0
	                    disp = pf126*eps*sigcub*sigcub
	                    repul = pf126*eps*sigcub*sigcub*sigcub*sigcub
                    else:
                        print("ERROR: Write correct LJ type {} {} (e.g.) lj12_4".format(database.vdwtype1[vdwtmp],database.vdwtype2[vdwtmp]))
                        sys.exit(0)

                    print("{:>6s} {:>6s}    1  {:14.5e} {:14.5e} ;  {} {:5.6f} {:5.6f}".format
                            (database.vdwtype1[vdwtmp], database.vdwtype2[vdwtmp], disp, repul,database.vdwstyle[vdwtmp],eps,sig), file=fout)
        
        # get bond interactions
        print(file=fout)
        if sysdat.nbnds > 0:
            print("[ bondtypes ]", file=fout)
            print("; i     j   funct   length  force.c", file=fout)
            for idx in range(sysdat.ntops):
                for jdx in range(topdat[idx].nbnd):
                    ikeep = 1
                    datndx = -1
                    # AT THIS POINT WE WILL CHECK TO SEE IF THE PARAMS WERE GIVEN */
                    # IN THE TOP FILE....IF SO WE WILL SKIP A LOT OF THIS AND ADD */
                    # THIS AS A UNIQUE BOND....IF NOT WE GO THROUGH THE PROCEDURE */
                    if topdat[idx].bndpset[jdx] == False:
                        # now compare to the database 
                        for kdx in range(database.nbndtype):
                                f1 = cmp_wc(database.bndtype1[kdx], topdat[idx].atomtype[topdat[idx].bndndx1[jdx]-1]) 
                                f2 = cmp_wc(database.bndtype2[kdx], topdat[idx].atomtype[topdat[idx].bndndx2[jdx]-1])
                                if f1 and f2:
                                    datndx = kdx
                                    b1tmp = topdat[idx].atomtype[topdat[idx].bndndx1[jdx]-1]
                                    b2tmp = topdat[idx].atomtype[topdat[idx].bndndx2[jdx]-1]
                                    break
                                f1 = cmp_wc(database.bndtype2[kdx], topdat[idx].atomtype[topdat[idx].bndndx1[jdx]-1]) 
                                f2 = cmp_wc(database.bndtype1[kdx], topdat[idx].atomtype[topdat[idx].bndndx2[jdx]-1])
                                if f1 and f2:
                                    datndx = kdx
                                    b2tmp = topdat[idx].atomtype[topdat[idx].bndndx1[jdx]-1]
                                    b1tmp = topdat[idx].atomtype[topdat[idx].bndndx2[jdx]-1]
                                    break
                        if datndx == -1:
                                sys.exit("ERROR: Did not find bond parameters in database {} {} {} {}".format(
                                        topdat[idx].bndndx1[jdx],
                                        topdat[idx].bndndx2[jdx],
                                        topdat[idx].atomtype[topdat[idx].bndndx1[jdx]-1],
                                        topdat[idx].atomtype[topdat[idx].bndndx2[jdx]-1]))

                        # Now make sure we do not already know we have this interaction
                        for kdx in range(uniq_bnds):
                            if bnd_name1[kdx] == b1tmp and bnd_name2[kdx] == b2tmp:
                                ikeep = 0
                                topdat[idx].bndtype.append(kdx)
                                break
                            if bnd_name2[kdx] == b1tmp and bnd_name1[kdx] == b2tmp:
                                ikeep = 0
                                topdat[idx].bndtype.append(kdx)
                                break
                        # ikeep = 1 if we found a new one 
                        if ikeep == 1:
                            bnd_params.append(datndx)
                            bnd_name1.append(b1tmp)
                            bnd_name2.append(b2tmp)
                            sysdat.param_bnds.append(datndx)
                            topdat[idx].bndtype.append(uniq_bnds)
                            uniq_bnds += 1
                            tmpcalca = database.bnde[bnd_params[uniq_bnds-1]]/10.0
                            tmpcalcb = database.fbnd[bnd_params[uniq_bnds-1]]*4.184*2.0*100
                            print("{} {}".format(database.bndtype1[bnd_params[uniq_bnds-1]],database.bndtype2[bnd_params[uniq_bnds-1]]))
                            print("{:>6s} {:>6s}    1  {:8.4f}  {:8.4f}".format(b1tmp,b2tmp,tmpcalca,tmpcalcb), file=fout)
        sysdat.uniq_nbnds = uniq_bnds

        # get angle interactions
        print(file=fout)
        if sysdat.nangs > 0:
            print("[ angletypes ]", file=fout)
            print("; i     j     k    funct   angle  force.c  sigma   eps", file=fout)
            for idx in range(sysdat.ntops):
                for jdx in range(topdat[idx].nang):
                    datndx = -1
                    # AT THIS POINT WE WILL CHECK TO SEE IF THE PARAMS WERE GIVEN */
                    # IN THE TOP FILE....IF SO WE WILL SKIP A LOT OF THIS AND ADD */
                    # THIS AS A UNIQUE BOND....IF NOT WE GO THROUGH THE PROCEDURE */
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
                        if ifound == 0:
                            print("*********************");
                            print("ERROR: No params for VDW interaction between {} and {} for angle (database)".format(topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1],
                                   topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1]))
                            print("Update database.")
                        # end VDW for CG angles 
                        # No params for this interaction in the database
                        if datndx == -1:
                            print("ERROR: Did not find angle parameters in database {} {} {} ({} {} {})".format(topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1],
                                   topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1],
                                   topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1],
                                   topdat[idx].angndx1[jdx],
                                   topdat[idx].angndx2[jdx],
                                   topdat[idx].angndx3[jdx]))
                            sys.exit(1)
                        # Now make sure we do not already have this one
                        ikeep = 1
                        for kdx in range(uniq_angs):
                            if datndx == ang_params[kdx]:
                                ikeep = 0
                                topdat[idx].angtype.append(kdx)
                                break
                        if ikeep == 1:
                            eps = database.eps[vdwtmp]*4.184
                            sig = database.sig[vdwtmp]/10.0
                            ang_params.append(datndx)
                            ang_vdw.append(vdwtmp)
                            sysdat.param_angs.append(datndx)
                            topdat[idx].angtype.append(uniq_angs)
                            uniq_angs += 1
                            tmpcalca = database.ange[ang_params[uniq_angs-1]]
                            tmpcalcb = database.fang[ang_params[uniq_angs-1]]*4.184*2.0
                            if tmpcalca > 0:
                                if database.angsdk[datndx] == 0:
                                    print("{:>6s} {:>6s} {:>6s}    1  {:8.4f} {:8.4f}".format(
                                            database.angtype1[ang_params[ uniq_angs-1]],
                                            database.angtype2[ang_params[ uniq_angs-1]],
                                            database.angtype3[ang_params[ uniq_angs-1]],
                                            tmpcalca, tmpcalcb), file=fout)
                                else:
                                    print("{:>6s} {:>6s} {:>6s}    4  {:8.4f} {:8.4f} {:8.4f} {:8.4f}".format(
                                            database.angtype1[ang_params[ uniq_angs-1]],
                                            database.angtype2[ang_params[ uniq_angs-1]],
                                            database.angtype3[ang_params[ uniq_angs-1]],
                                            tmpcalca, tmpcalcb, sig, eps), file=fout)
        
        sysdat.uniq_nangs = uniq_angs
        
        # Now let handle the dihedral params the only way we do... they have to be specified in the top file 
        print(file=fout)

# Read the database file and store unique params
# Warn if you find duplicates
def read_database(fname, database):
    nvdw = nbnd = nang = 0
    with open(fname, "r") as fin:
        line = fin.readline()
        while line:
            items = line.split()
            if len(items) == 0:
                line = fin.readline()
                continue
            if items[0] == "pair":
                ikeep = 1
                vdwtype1 = items[1]
                vdwtype2 = items[2]
                vdwstyle = items[3]
                eps      = float(items[4])
                sig      = float(items[5])
                for idx in range(nvdw):
                    if   vdwtype1 == database.vdwtype1[idx] and vdwtype2 == database.vdwtype2[idx]:
                        print("WARNING: Found dup vdw param {} {}".format(vdwtype1,vdwtype2))
                        ikeep = 0
                    elif vdwtype1 == database.vdwtype2[idx] and vdwtype2 == database.vdwtype1[idx]:
                        print("WARNING: Found dup vdw param {} {}".format(vdwtype1,vdwtype2))
                        ikeep = 0
                if ikeep == 1:
                        database.vdwtype1.append(vdwtype1)
                        database.vdwtype2.append(vdwtype2)
                        database.vdwstyle.append(vdwstyle)
                        database.eps.append(eps)
                        database.sig.append(sig)
                        nvdw += 1
            if items[0] == "bond":
                ikeep = 1
                bndtype1 = items[1]
                bndtype2 = items[2]
                fbnd     = float(items[3])
                bnde     = float(items[4])
                for idx in range(nbnd):
                        if   bndtype1 == database.bndtype1[idx] and bndtype2 == database.bndtype2[idx]:
                            print("WARNING: Found dup bond param {} {}".format(bndtype1,bndtype2))
                            ikeep = 0
                        elif bndtype1 == database.bndtype2[idx] and bndtype2 == database.bndtype1[idx]:
                            print("WARNING: Found dup bond param {} {}".format(bndtype1,bndtype2))
                            ikeep = 0
                if ikeep == 1:
                    database.bndtype1.append(bndtype1)
                    database.bndtype2.append(bndtype2)
                    database.fbnd.append(fbnd)
                    database.bnde.append(bnde)
                    nbnd += 1
            if items[0] == "angle":
                if "harmonic" in line:
                    angsdk = 0
                else:
                    angsdk = 1
                ikeep = 1
                angtype1 = items[1]
                angtype2 = items[2]
                angtype3 = items[3]
                fang     = float(items[4])
                ange     = float(items[5])
                for idx in range(nang):
                    if angtype2 == database.angtype2[idx]:
                        if   angtype1 == database.angtype1[idx] and angtype3 == database.angtype3[idx]:
                                print("WARNING: Found dup angle param {} {} {}".format(angtype1,angtype2,angtype3))
                                ikeep = 0
                        elif angtype3 == database.angtype1[idx] and angtype1 == database.angtype3[idx]:
                                print("WARNING: Found dup angle param {} {} {}".format(angtype1,angtype2,angtype3))
                                ikeep = 0
                if ikeep == 1:
                    database.angtype1.append(angtype1)
                    database.angtype2.append(angtype2)
                    database.angtype3.append(angtype3)
                    database.fang.append(fang)
                    database.ange.append(ange)
                    database.angsdk.append(angsdk)
                    nang += 1
            line = fin.readline()
        database.nvdwtype = nvdw
        database.nbndtype = nbnd
        database.nangtype = nang

# Count the number of params in the database so we can allocate for storage
def count_params(fname, database):
    database.nvdwtype = 0
    database.nbndtype = 0
    database.nangtype = 0
    with open(fname, "r") as fin:
        line = fin.readline()
        while line:
            items = line.split()
            if len(items) == 0:
                line = fin.readline()
                continue
            if items[0] == "pair":
                database.nvdwtype += 1
            if items[0] == "bond":
                database.nbndtype += 1
            if items[0] == "angle":
                database.nangtype += 1
            line = fin.readline()
        
# count the number of things in the topology files so we can allocate
def count_atoms(fname, topdat, ntop):
    topdat[ntop].nat     = 0
    topdat[ntop].nbnd    = 0
    topdat[ntop].nang    = 0
    topdat[ntop].ndihed  = 0
    topdat[ntop].nimprop = 0
    with open(fname, "r") as fin:
        line = fin.readline()
        while line:
            items = line.split()
            if len(items) == 0:
                line = fin.readline()
                continue
            if items[0] == "atom":
                topdat[ntop].nat += 1
            if items[0] == "goparam":
                topdat[ntop].ngo += 1
            if items[0] == "bond"  or items[0] == "bondparam":
                topdat[ntop].nbnd += 1
            if items[0] == "angle" or items[0] == "angleparam":
                topdat[ntop].nang += 1
            if items[0] == "dihedral" or items[0] == "dihedralparam":
                topdat[ntop].ndihed  += 1
            if items[0] == "improper" or items[0] == "improperparam":
                topdat[ntop].nimprop += 1
            line = fin.readline()
        if topdat[ntop].nat == 0:
            sys.exit("ERROR: natom in {} is zero.".format(fname))

# Read the topology file and store the data
def read_top(fname, sysdat, topdat, ntop):
    log_bndprm = log_angprm = log_dihprm = log_impprm = log_charge = True
    ndx = bndx = andx = dndx = indx = lc = 0
    print("######################")
    print("##### READING {}".format(fname))
    with open(fname, "r") as fin:
        line = fin.readline()
        while line:
            lc    += 1
            items = line.split()
            if len(items) == 0:
                line = fin.readline()
                continue
            if items[0] == "atom":
                try:
                    topdat[ntop].ind.append(int(items[1]))
                    topdat[ntop].resname.append(items[2])
                    topdat[ntop].atomname.append(items[3])
                    topdat[ntop].atomtype.append(items[4])
                    topdat[ntop].mass.append(float(items[5]))
                    topdat[ntop].charge.append(float(items[6].replace("+","")))
                    topdat[ntop].segid.append(items[7])
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                if topdat[ntop].charge[ndx]*topdat[ntop].charge[ndx] > 1e-5 and log_charge:
                    print("Charge in top file {} {}".format(fname, topdat[ntop].charge[ndx]))
                    sysdat.ischarged = 1
                    log_charge = False
                ndx += 1
            # Go model for protein backbones
            if items[0] == "goparam":
                try:
                    topdat[ntop].gondx1.append(int(items[1]))
                    topdat[ntop].gondx2.append(int(items[2]))
                    topdat[ntop].gofunctype.append((items[3]))
                    topdat[ntop].eps.append((float(items[4])))
                    topdat[ntop].sig.append((float(items[5])))
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))

            if items[0] == "bond":
                try:
                    topdat[ntop].bndndx1.append(int(items[1]))
                    topdat[ntop].bndndx2.append(int(items[2]))
                    topdat[ntop].bndfk.append(None)
                    topdat[ntop].bndeq.append(None)
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].bndpset.append(False)
                bndx += 1
            if items[0] == "bondparam":
                if log_bndprm:
                    print("WARNING: Using bond parameters from the top file.")
                    log_bndprm = False
                if len(items) < 5:
                    sys.exit("ERROR: Not enough args for bondparam: must be: ndx1 ndx2 fk eq.")
                try:
                    topdat[ntop].bndndx1.append(int(items[1]))
                    topdat[ntop].bndndx2.append(int(items[2]))
                    topdat[ntop].bndfk.append(float(items[3]))
                    topdat[ntop].bndeq.append(float(items[4]))
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].bndpset.append(True)
                bndx += 1
            if items[0] == "angle":
                try:
                    topdat[ntop].angndx1.append(int(items[1]))
                    topdat[ntop].angndx2.append(int(items[2]))
                    topdat[ntop].angndx3.append(int(items[3]))
                    topdat[ntop].angfk.append(None)
                    topdat[ntop].angeq.append(None)
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].angpset.append(False)
                andx += 1
            if items[0] == "angleparam":
                if log_angprm:
                    print("WARNING: Using angle parameters from the top file.")
                    log_angprm = False
                if len(items) < 6:
                    sys.exit("ERROR: Not enough args for angleparam: must be: ndx1 ndx2 ndx3 fk eq.")
                try:
                    topdat[ntop].angndx1.append(int(items[1]))
                    topdat[ntop].angndx2.append(int(items[2]))
                    topdat[ntop].angndx3.append(int(items[3]))
                    topdat[ntop].angfk.append(float(items[4]))
                    topdat[ntop].angeq.append(float(items[5]))
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].angpset.append(True)
                andx += 1
            if items[0] == "improper":
                print("WARNING: This is not implemented. must use improperparam and assign improper parameters in the top file.")
                try:
                    topdat[ntop].impropndx1.append(int(items[1]))
                    topdat[ntop].impropndx2.append(int(items[2]))
                    topdat[ntop].impropndx3.append(int(items[3]))
                    topdat[ntop].impropndx4.append(int(items[4]))
                    topdat[ntop].impropfk.append(float(items[5]))
                    topdat[ntop].impropeq.append(float(items[6]))
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].improppset.append(False)
                indx += 1
            if items[0] == "improperparam":
                if log_impprm:
                    print("WARNING: Using improper parameters from the top file.")
                    log_impprm = False
                if len(items) < 7:
                    sys.exit("ERROR: Not enough args for improperparam: must be: ndx1 ndx2 ndx3 ndx4 fk eq.")
                try:
                    topdat[ntop].impropndx1.append(int(items[1]))
                    topdat[ntop].impropndx2.append(int(items[2]))
                    topdat[ntop].impropndx3.append(int(items[3]))
                    topdat[ntop].impropndx4.append(int(items[4]))
                    topdat[ntop].impropfk.append(float(items[5]))
                    topdat[ntop].impropeq.append(float(items[6]))
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].improppset.append(True)
                indx += 1
            if items[0] == "dihedralparam":
                if log_dihprm:
                    print("WARNING: Using dihedral parameters from the top file.")
                    log_dihprm = False
                if len(items) < 9:
                    sys.exit("ERROR: Not enough args for angleparam: must be: ndx1 ndx2 ndx3 fk n eq onefour.")
                try:
                    topdat[ntop].dihedndx1.append(int(items[1]))
                    topdat[ntop].dihedndx2.append(int(items[2]))
                    topdat[ntop].dihedndx3.append(int(items[3]))
                    topdat[ntop].dihedndx4.append(int(items[4]))
                    topdat[ntop].dihedfk.append(float(items[5]))
                    topdat[ntop].dihedn.append(int(items[6]))
                    topdat[ntop].dihedeq.append(int(float(items[7])))
                    topdat[ntop].dihedof.append(float(items[8]))
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].dihedpset.append(True)
                dndx += 1
            line = fin.readline()

def make_ndx(database,topdat,sysdat):
    ljtypes = {}
    grps = ["lj12_4", "lj9_6", "sol", "psol"]
    sols = ["W", "SOD", "CLA", "POT", "MAG", "CAL"]
    psolsp = ["WO", "SOD1", "CLA1"]
    psolsa= ["WO", "WH", "SOD1", "SOD2", "CLA1", "CLA2"]
    for kdx in range(database.nvdwtype):
        vt1 = database.vdwtype1[kdx]
        vt2 = database.vdwtype2[kdx] 
        if vt1 in sols and vt2 in sols:
            ljtypes[vt1] = "sol"
        elif vt1 in psolsa and vt2 in psolsa:
            ljtypes[vt1] = "psol"
        elif vt1 in sols or vt1 in psolsp:
            ljtypes[vt2] = database.vdwstyle[kdx]
        elif vt2 == sols or vt2 in psolsp:
            ljtypes[vt1] = database.vdwstyle[kdx]

    with open("CGindex.ndx", "w") as fpind:
        print(f"[ system ]", file=fpind)
        atindex = 0
        indcnt = 0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nmol):
                for kdx in range(topdat[idx].nat):
                    atindex += 1
                    print("{:5d} ".format(atindex), end="", file=fpind)
                    indcnt += 1
                    if indcnt%15 == 0:
                        print(file=fpind)
        print(file=fpind)
        print(file=fpind)
        for grp in grps:
            grp_flag = True
            atindex = 0
            indcnt = 0
            for idx in range(sysdat.ntops):
                for jdx in range(topdat[idx].nmol):
                    for kdx in range(topdat[idx].nat):
                        atindex += 1
                        if ljtypes[topdat[idx].atomtype[kdx]] == grp:
                            if grp_flag:
                                print(f"[ {re.sub('_', '', grp).upper()}W ]", file=fpind)
                                grp_flag = False
                            print("{:5d} ".format(atindex), end="", file=fpind)
                            indcnt += 1
                            if indcnt%15 == 0:
                                print(file=fpind)
            if grp_flag == False:
                print(file=fpind)
                print(file=fpind)

def make_top(database,topdat,sysdat):
    with open("topol.top", "w") as fout:
        print("; generated by cg_spica setup_gmax", file=fout)
        print("#include \"toppar/SPICA.itp\"", file=fout)
        for idx in range(sysdat.ntops):
          print(f"#include \"toppar/{topdat[idx].fname}.itp\"", file=fout)
        print(file=fout)
        print("[ system ]", file=fout)
        print("; Name", file=fout)
        print("CG", file=fout)
        print(file=fout)
        print("[ molecules ]", file=fout)
        print("; Compound   #mols", file=fout)
        for idx in range(sysdat.ntops):
          print("{}     {}".format(topdat[idx].fname, topdat[idx].nmol), file=fout)

# Main routine. Call and allocate                                       
# The idea is to read in the topologies and then check the database for 
# all of the required interaction params.                              

def run(args):
    inputs = args.input_files
    nargs = len(inputs)
    if args.prot:
        if nargs < 4:
            print("Dumps input files for a GROMACS run.")
            print("usage: setup_gmx -p <topfile 1> <nmol 1> [ <topfile 2> <nmol 2> ..... <topfile n> <nmol n>] <database> <pdbfile>");
            print("Takes at least four arguments (one component system): 1) Topology, 2) number of molecules, 3) parameter database, 4) PDB.")
            sys.exit(1)
        print("setup_gmx for SPICA protein model.")
        ntops = int((nargs-2)/2)
    else:
        if nargs < 3 or nargs%2  == 0:
            print("Dumps input files for a GROMACS run.")
            print("usage: setup_gmx [-p] <topfile 1> <nmol 1> [ <topfile 2> <nmol 2> ..... <topfile n> <nmol n>] <database>")
            print("Takes at least three arguments (one component system): 1) Topology, 2) number of molecules, 3) parameter database.")
            sys.exit(1)
        print("setup_gmx for SPICA.")
        ntops = int((nargs-1)/2)
    topdat   = [Topdat() for _ in range(ntops)]
    database = Database()
    sysdat   = Sysdat()
    print("Will read {} topology file(s).".format(ntops))

    # Loop through the topologies and count the number of atoms, bonds and bends
    sysdat.ntops = ntops
    sysdat.nats = 0
    sysdat.nbnds = 0
    sysdat.nangs = 0
    sysdat.total_ats = 0
    sysdat.total_bnds = 0
    sysdat.total_angs = 0
    sysdat.total_improps = 0
    sysdat.total_diheds = 0
    for idx in range(ntops):
        topfile = inputs[2*idx]
        ptop = Path(topfile)
        topdat[idx].fname = re.sub(ptop.suffix, "", ptop.name)
        topdat[idx].nmol = int(inputs[2*idx+1])
        count_atoms(topfile, topdat,idx)
        read_top(topfile, sysdat, topdat, idx)
        print("Bookkeeping:")
        print("Found: {} atoms".format(topdat[idx].nat))
        print("Found: {} bonds".format(topdat[idx].nbnd))
        print("Found: {} angles".format(topdat[idx].nang))
        print("Found: {} impropers".format(topdat[idx].nimprop))
        sysdat.nats	     += topdat[idx].nat
        sysdat.nbnds     += topdat[idx].nbnd
        sysdat.nangs     += topdat[idx].nang
        sysdat.nimprops  += topdat[idx].nimprop
        sysdat.ndiheds   += topdat[idx].ndihed
        sysdat.total_ats	 += topdat[idx].nat*topdat[idx].nmol
        sysdat.total_bnds	 += topdat[idx].nbnd*topdat[idx].nmol
        sysdat.total_angs	 += topdat[idx].nang*topdat[idx].nmol
        sysdat.total_improps += topdat[idx].nimprop*topdat[idx].nmol
        sysdat.total_diheds	 += topdat[idx].ndihed*topdat[idx].nmol
    if args.prot:
        count_params(inputs[nargs-2], database)
        read_database(inputs[nargs-2], database)
    else:
        count_params(inputs[nargs-1], database)
        read_database(inputs[nargs-1], database)
    print("Found {} unique vdw pair params".format(database.nvdwtype))
    print("Found {} unique bond params".format(database.nbndtype))
    print("Found {} unique angle params".format(database.nangtype))
    if args.prot:
        if ".pdb" in inputs[nargs-1]:
            print("Takes angles from {}".format(inputs[nargs-1]))
            read_pdb(inputs[nargs-1],sysdat)
        else:
            print("No PDB files to take angles!")
            sys.exit(1)
    get_unique(database, topdat, sysdat)
    read_coords(database, topdat, sysdat)
    make_top(database, topdat, sysdat)
    make_ndx(database, topdat, sysdat)
    write_psf(database, topdat, sysdat)

if __name__ == "__main__":
    args = get_option()
    run(args)
