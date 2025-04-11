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
    r12 = r1 - r2
    r32 = r3 - r2
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
    r12 = r1 - r2
    r32 = r3 - r2
    r43 = r4 - r3
    r23 = r2 - r3
    r13_inn = np.dot(r12, r32)
    r42_inn = np.dot(r43, r23)
    r13_mag = np.linalg.norm(r12)*np.linalg.norm(r32)
    r42_mag = np.linalg.norm(r43)*np.linalg.norm(r23)
    cos13 = r13_inn/r13_mag
    cos42 = r42_inn/r42_mag
    p1 = r12 - cos13*r32/np.linalg.norm(r32)*np.linalg.norm(r12)
    p4 = r43 - cos42*r23/np.linalg.norm(r23)*np.linalg.norm(r43)
    cosp = np.dot(p1, p4)/(np.linalg.norm(p1)*np.linalg.norm(p4))
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
    nats = nbnds = nangs = nimps = ndihs = ntops = 0 
    total_ats = total_bnds = total_angs = total_imps = total_dihs = 0
    foundatoms = boxinfo = ischarged = 0
    uniq_nats = uniq_nbnds = uniq_nangs = uniq_nimps = uniq_ndihs = 0
    param_bnds, param_angs, param_dihs, param_imps  = [], [], [], []
    coordx, coordy, coordz = [], [], []
    boxx = boxy = boxz = 0.0


class Topdat:
    def __init__(self):
        self.fname = None
        self.nat = self.nbnd = self.nang = self.nimp = self.nmol = self.ngo = 0
        self.gondx1, self.gondx2, self.gofunctype, self.eps, self.sig = [], [], [], [], []
        self.bndndx1, self.bndndx2, self.bndtype = [], [], []
        self.angndx1, self.angndx2, self.angndx3, self.angtype = [], [], [], []
        self.impndx1, self.impndx2, self.impndx3, self.impndx4 = [], [], [], []
        self.dihndx1, self.dihndx2, self.dihndx3, self.dihndx4 = [], [], [], []
        self.dih_func, self.imp_func, self.dihtype, self.imptype, self.dihn = [], [], [], [], []
        self.dihpset, self.imppset, self.bndpset, self.angpset = [], [], [], []
        self.ind, self.parm_atomtype = [], []
        self.dihfk, self.ndih, self.diheq = [], [], []
        self.mass, self.charge = [], []
        self.bndfk, self.bndeq = [], []
        self.angfk, self.angeq = [], []
        self.impfk, self.impeq = [], []
        self.atomname, self.atomtype, self.segid, self.resname = [], [], [], []
   

class Database:
    fbnd, bnde, fang, ange, eps, sig, angsdk = [], [], [], [], [], [], []
    fdih, dihn, dihe, fimp, impe = [], [], [], [], []
    nvdwtype, nbndtype, nangtype, ndihtype, nimptype = [], [], [], [], []
    vdwtype1, vdwtype2, vdwstyle = [], [], []
    bndtype1, bndtype2 = [], []
    angtype1, angtype2, angtype3 = [], [], []
    dihtype1, dihtype2, dihtype3, dihtype4 = [], [], [], []
    imptype1, imptype2, imptype3, imptype4 = [], [], [], []
    loop_pair = []


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
                print("{:6d} {:>6s}    1    {:>6s} {:>6s} {:6d} {:8.4f}".format(
                       jdx+1, 
                       topdat[idx].atomtype[jdx],
                       topdat[idx].resname[jdx], 
                       topdat[idx].atomname[jdx],
                       jdx+1,
                       topdat[idx].charge[jdx]), file=fout)

            if topdat[idx].nbnd > 0:
                print(file=fout)
                print("[ bonds ]", file=fout)
                print(file=fout)
                for jdx in range(topdat[idx].nbnd):
                    if topdat[idx].bndpset[jdx]:
                        tmpcalca = topdat[idx].bndeq[jdx]*0.1
                        tmpcalcb = topdat[idx].bndfk[jdx]*4.184*2.0*100
                        print("{:5d} {:5d}    6  {:8.4f}  {:8.4f} ; EN specified by bondparam".format(
                               topdat[idx].bndndx1[jdx],
                               topdat[idx].bndndx2[jdx],
                               tmpcalca, tmpcalcb), file=fout)
                    else :
                        print("{:5d} {:5d}    1".format(
                               topdat[idx].bndndx1[jdx],
                               topdat[idx].bndndx2[jdx]), file=fout)
                if topdat[idx].resname[0] in ["PWAT", "PSOD", "PCLA"]:
                    print(file=fout)
                    print("[ constraints ]", file=fout)
                    print(file=fout)
                    print("{:5d} {:5d}    1  0.11".format(
                           topdat[idx].bndndx1[0],
                           topdat[idx].bndndx2[0]), file=fout)

            # Go model for protein backbones
            if topdat[idx].ngo > 0:
                print(file=fout)
                print("[ pairs ]", file=fout)
                print("; Go model for protein backbone", file=fout)
                print("; ai   aj    funct    c6    c12", file=fout)
                for jdx in range(topdat[idx].ngo):
                    tmpcalcb = topdat[idx].sig[jdx]*0.1
                    print("{:5d} {:5d}    1  {:14.5e} {:14.5e};  {} {:5.6f} {:5.6f}".format(
                           topdat[idx].gondx1[jdx],
                           topdat[idx].gondx2[jdx],
                           4.184*4.0*topdat[idx].eps[jdx]*pow(tmpcalcb, 6),
                           4.184*4.0*topdat[idx].eps[jdx]*pow(tmpcalcb, 12),
                           topdat[idx].gofunctype[jdx],
                           topdat[idx].eps[jdx],
                           tmpcalcb), file=fout)
                print(file=fout)

            if topdat[idx].nang > 0:
                print(file=fout)
                print("[ angles ]",file=fout)
                print(file=fout)
                for jdx in range(topdat[idx].nang):
                    if topdat[idx].angpset[jdx]:
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
                               topdat[idx].angndx1[jdx], 
                               topdat[idx].angndx2[jdx], 
                               topdat[idx].angndx3[jdx],
                               tmpcalca, tmpcalcb,
                               topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1],
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
                            print("{:5d} {:5d} {:5d}    1  {:8.4f} {:8.4f} ; {:>6s} {:>6s} {:>6s} taken from pdb".format(
                                   topdat[idx].angndx1[jdx],
                                   topdat[idx].angndx2[jdx],
                                   topdat[idx].angndx3[jdx],
                                   angle_in_pdb, tmpcalcb, 
                                   database.angtype1[datndx],
                                   database.angtype2[datndx],
                                   database.angtype3[datndx]), file=fout)
                        else:
                            if database.angsdk[datndx]:
                                print("{:5d} {:5d} {:5d}    4".format(
                                       topdat[idx].angndx1[jdx],
                                       topdat[idx].angndx2[jdx],
                                       topdat[idx].angndx3[jdx]), file=fout)
                            else:
                                print("{:5d} {:5d} {:5d}    1".format(
                                       topdat[idx].angndx1[jdx],
                                       topdat[idx].angndx2[jdx],
                                       topdat[idx].angndx3[jdx]), file=fout)

            if topdat[idx].ndih > 0:
                print(file=fout)
                print("[ dihedrals ]",file=fout)
                print(file=fout)
                for jdx in range(topdat[idx].ndih):
                    if topdat[idx].dihpset[jdx]:
                        print("{:5d} {:5d} {:5d} {:5d}  1  {:<3f} {:8.4f} {:<3d} ; FROM TOP".format(
                               topdat[idx].dihndx1[jdx],
                               topdat[idx].dihndx2[jdx],
                               topdat[idx].dihndx3[jdx],
                               topdat[idx].dihndx4[jdx],
                               topdat[idx].diheq[jdx],
                               topdat[idx].dihfk[jdx]*4.184,
                               topdat[idx].dihn[jdx]), file=fout)
                    else:
                        print("{:5d} {:5d} {:5d} {:5d}  9".format(
                               topdat[idx].dihndx1[jdx],
                               topdat[idx].dihndx2[jdx],
                               topdat[idx].dihndx3[jdx],
                               topdat[idx].dihndx4[jdx]), file=fout)

            if topdat[idx].nimp > 0:
                print(file=fout)
                print("[ dihedrals ]", file=fout)
                print(file=fout)
                for jdx in range(topdat[idx].nimp):
                    if topdat[idx].imppset[jdx]:
                        print("{:5d} {:5d} {:5d} {:5d}  2  {:8.4f} {:8.4f} ; FROM TOP".format(
                               topdat[idx].impndx1[jdx],
                               topdat[idx].impndx2[jdx],
                               topdat[idx].impndx3[jdx],
                               topdat[idx].impndx4[jdx],
                               topdat[idx].impfk[jdx]*4.184*2.0,
                               topdat[idx].impeq[jdx]), file=fout)
                    else:
                        print("{:5d} {:5d} {:5d} {:5d}  2".format(
                               topdat[idx].impndx1[jdx],
                               topdat[idx].impndx2[jdx],
                               topdat[idx].impndx3[jdx],
                               topdat[idx].impndx4[jdx]), file=fout)

            # Remove non-native vdw interaction between protein backbone beads forming native contact
            if True in topdat[idx].bndpset or topdat[idx].ngo > 0:
                print(file=fout)
                print("[ exclusions ]", file=fout)
                if topdat[idx].ngo > 0:
                    print(file=fout)
                    for jdx in range(topdat[idx].ngo):
                        if jdx == 0:
                            print("{:4d} {:4d} ".format(
                                   topdat[idx].gondx1[jdx],
                                   topdat[idx].gondx2[jdx]), end='', file=fout)
                        elif topdat[idx].gondx1[jdx] != topdat[idx].gondx1[jdx-1]:
                            print(file=fout)
                            print("{:4d} {:4d} ".format(
                                   topdat[idx].gondx1[jdx],
                                   topdat[idx].gondx2[jdx]), end='', file=fout)
                        else:
                            print("{:4d} ".format(
                                   topdat[idx].gondx2[jdx]), end='', file=fout)
                if True in topdat[idx].bndpset:
                    tmp_bndndx1 = 0
                    for jdx in range(topdat[idx].nbnd):
                        if topdat[idx].bndpset[jdx] == True:
                            if topdat[idx].bndndx1[jdx] != tmp_bndndx1:
                                tmp_bndndx1 = topdat[idx].bndndx1[jdx]
                                print(file=fout)
                                print("{:4d} {:4d} ".format(
                                       topdat[idx].bndndx1[jdx],
                                       topdat[idx].bndndx2[jdx]), end='', file=fout)
                            else:
                                print("{:4d} ".format(
                                       topdat[idx].bndndx2[jdx]), end='', file=fout)
            print(file=fout)


def write_psf(topdat, sysdat):
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
                    print("{:8} {:<4}{:5} {:<4} {:<4} {:<4}  {:9.6f}  {:12.4f}".format(
                           atidx, topdat[idx].resname[kdx], molidx % 10000, 
                           topdat[idx].resname[kdx], topdat[idx].atomname[kdx], topdat[idx].atomtype[kdx], 
                           topdat[idx].charge[kdx], topdat[idx].mass[kdx]), file=fout)
        print(file=fout)
        print("{:8} !NBOND: bonds".format(sysdat.total_bnds), file=fout)
        bondidx = offset  = 0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nmol):
                for kdx in range(topdat[idx].nbnd):
                    bondidx += 1
                    print("{:>8}{:>8}".format(
                           topdat[idx].bndndx1[kdx]+(jdx*topdat[idx].nat)+offset,
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
                    print("{:>8}{:>8}{:>8}".format(
                           topdat[idx].angndx1[kdx]+(jdx*topdat[idx].nat)+offset,
                           topdat[idx].angndx2[kdx]+(jdx*topdat[idx].nat)+offset,
                           topdat[idx].angndx3[kdx]+(jdx*topdat[idx].nat)+offset),
                           file=fout, end="")
                    if angleidx % 3 ==0:
                        print(file=fout)
            offset += topdat[idx].nmol*topdat[idx].nat
        print(file=fout)
        print(file=fout)
        print("{:8} !NPHI: dihedrals".format(sysdat.total_dihs), file=fout)
        dihedidx = offset = 0;
        dihs_lst = []
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nmol):
                for kdx in range(topdat[idx].ndih):
                    dihedidx += 1
                    d1 = topdat[idx].dihndx1[kdx]+(jdx*topdat[idx].nat)+offset
                    d2 = topdat[idx].dihndx2[kdx]+(jdx*topdat[idx].nat)+offset
                    d3 = topdat[idx].dihndx3[kdx]+(jdx*topdat[idx].nat)+offset
                    d4 = topdat[idx].dihndx4[kdx]+(jdx*topdat[idx].nat)+offset
                    dihs_tmp = [ d1, d2, d3, d4 ]
                    if dihs_tmp in dihs_lst:
                        continue
                    dihs_lst.append(dihs_tmp)
                    print("{:>8}{:>8}{:>8}{:>8}".format(d1, d2, d3, d4),
                           file=fout, end="")
                    if dihedidx % 2 ==0:
                        print(file=fout)
            offset += topdat[idx].nmol*topdat[idx].nat
        print(file=fout)
        print(file=fout)
        print("{:8} !NIMPHI: impropers".format(sysdat.total_imps), file=fout)
        improidx = offset = 0;
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nmol):
                for kdx in range(topdat[idx].nimp):
                    improidx += 1
                    print("{:>8}{:>8}{:>8}{:>8}".format(
                           topdat[idx].impndx1[kdx]+(jdx*topdat[idx].nat)+offset,
                           topdat[idx].impndx2[kdx]+(jdx*topdat[idx].nat)+offset,
                           topdat[idx].impndx3[kdx]+(jdx*topdat[idx].nat)+offset,
                           topdat[idx].impndx4[kdx]+(jdx*topdat[idx].nat)+offset),
                           file=fout, end="")
                    if improidx % 2 ==0:
                        print(file=fout)
            offset += topdat[idx].nmol*topdat[idx].nat
        print(file=fout)
        print(file=fout)
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
    uniq_nats = uniq_bnds = uniq_angs = uniq_imps = uniq_dihs = 0
    uniq_atype, uniq_mass, uniq_charge = [], [], []
    bnd_params, bnd_name1, bnd_name2 = [], [], []
    ang_params, ang_vdw = [], []
    dih_params, imp_params = [], []
    os.makedirs("toppar", exist_ok=True)
    with open("toppar/SPICA.itp","w") as fout:
        print("; Generated by setup_gmx", file=fout)
        print(file=fout)
        # first gather the unique atom types
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nat):
                keep = True
                for kdx in range(uniq_nats):
                    if topdat[idx].atomtype[jdx] == uniq_atype[kdx]:
                        keep = False
                        topdat[idx].parm_atomtype.append(kdx)
                        break
                if keep:
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
            print("{:>8s} {:8.4f} {:8.4f}    A    0.0    0.0".format(
                   uniq_atype[idx], uniq_mass[idx], uniq_charge[idx]), file=fout)
        print(file=fout)
        # get pair interactions
        bb_sec = ['GBML','GBBL','ABBL','GBTL','ABTL','GBMS','GBBS','ABBS','GBTS','ABTS']
        print("[ nonbond_params ]", file=fout);
        print("; i     j    func   C    A", file=fout);
        for idx in range(uniq_nats):
            for jdx in range(idx, uniq_nats):
                if (uniq_atype[idx] in bb_sec and uniq_atype[jdx] in database.loop_pair) \
                    or (uniq_atype[jdx] in bb_sec and uniq_atype[idx] in database.loop_pair):
                    tmp_type1 = uniq_atype[idx]
                    tmp_type2 = uniq_atype[jdx]
                else:
                    if uniq_atype[idx] in ['GBTP','GBTN','ABTP','ABTN']:
                        tmp_type1 = uniq_atype[idx]
                    elif uniq_atype[idx][0:3] in ['GBM','GBB','GBT','ABB','ABT']:
                        tmp_type1 = uniq_atype[idx][0:3]
                    else:
                        tmp_type1 = uniq_atype[idx]
                    if uniq_atype[jdx] in ['GBTP','GBTN','ABTP','ABTN']:
                        tmp_type2 = uniq_atype[jdx]
                    elif uniq_atype[jdx][0:3] in ['GBM','GBB','GBT','ABB','ABT']:
                        tmp_type2 = uniq_atype[jdx][0:3]
                    else:
                        tmp_type2 = uniq_atype[jdx]
                found = False;
                for kdx in range(database.nvdwtype):
                    if database.vdwtype1[kdx] == tmp_type1 and database.vdwtype2[kdx] == tmp_type2:
                        found = True
                        vdwtmp = kdx
                        break
                    elif database.vdwtype2[kdx] == tmp_type1 and database.vdwtype1[kdx] == tmp_type2:
                        found = True
                        vdwtmp = kdx
                        break
                if found:
                    eps = database.eps[vdwtmp]
                    sig = database.sig[vdwtmp]*0.1
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
                        sys.exit("ERROR: Write correct LJ type {} {} (e.g.) lj12_4".format(
                                  database.vdwtype1[vdwtmp],database.vdwtype2[vdwtmp]))
                    print("{:>6s} {:>6s}    1  {:14.5e} {:14.5e} ;  {} {:5.6f} {:5.6f}".format(
                            uniq_atype[idx], uniq_atype[jdx], 
                            disp, repul, database.vdwstyle[vdwtmp], 
                            eps, sig), file=fout)
                else:
                    print("*********************")
                    print("WARNING: No params for VDW interaction between {} and {}".format(uniq_atype[idx],uniq_atype[jdx]))
                    print("Update database")
                    sys.exit(1)
        
        # get bond interactions
        print(file=fout)
        if sysdat.nbnds > 0:
            print("[ bondtypes ]", file=fout)
            print("; i     j     func     b0     kb", file=fout)
            for idx in range(sysdat.ntops):
                for jdx in range(topdat[idx].nbnd):
                    keep = True
                    datndx = -1
                    if not topdat[idx].bndpset[jdx]:
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
                                keep = False
                                topdat[idx].bndtype.append(kdx)
                                break
                            if bnd_name2[kdx] == b1tmp and bnd_name1[kdx] == b2tmp:
                                keep = False
                                topdat[idx].bndtype.append(kdx)
                                break
                        # keep = True if we found a new one 
                        if keep:
                            bnd_params.append(datndx)
                            bnd_name1.append(b1tmp)
                            bnd_name2.append(b2tmp)
                            sysdat.param_bnds.append(datndx)
                            topdat[idx].bndtype.append(uniq_bnds)
                            uniq_bnds += 1
                            tmpcalca = database.bnde[bnd_params[uniq_bnds-1]]/10.0
                            tmpcalcb = database.fbnd[bnd_params[uniq_bnds-1]]*4.184*2.0*100
                            print("{:>6s} {:>6s}    1  {:8.4f}  {:8.4f}".format(b1tmp,b2tmp,tmpcalca,tmpcalcb), file=fout)
        sysdat.uniq_nbnds = uniq_bnds

        # get angle interactions
        print(file=fout)
        if sysdat.nangs > 0:
            print("[ angletypes ]", file=fout)
            print("; i     j     k     func     theta0     ktheta     sigma     eps", file=fout)
            for idx in range(sysdat.ntops):
                for jdx in range(topdat[idx].nang):
                    datndx = -1
                    if not topdat[idx].angpset[jdx]:
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
                        found = False
                        for kdx in range(database.nvdwtype):
                            f1 = database.vdwtype1[kdx] == topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1]
                            f2 = database.vdwtype2[kdx] == topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1]
                            f3 = database.vdwtype1[kdx] == topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1]
                            f4 = database.vdwtype2[kdx] == topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1]
                            if f1 and f2:
                                found = True
                                vdwtmp = kdx
                                break
                            elif f3 and f4:
                                found = True
                                vdwtmp = kdx
                                break
                        if not found:
                            print("*********************");
                            print("ERROR: No params for VDW interaction between {} and {} for angle (database)".format(
                                   topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1],
                                   topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1]))
                            print("Update database")
                            sys.exit(1)
                        # end VDW for CG angles 
                        # No params for this interaction in the database
                        if datndx == -1:
                            sys.exit("ERROR: Did not find angle parameters in database {} {} {} ({} {} {})".format(
                                   topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1],
                                   topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1],
                                   topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1],
                                   topdat[idx].angndx1[jdx],
                                   topdat[idx].angndx2[jdx],
                                   topdat[idx].angndx3[jdx]))
                        # Now make sure we do not already have this one
                        keep = True
                        for kdx in range(uniq_angs):
                            if datndx == ang_params[kdx]:
                                keep = False
                                topdat[idx].angtype.append(kdx)
                                break
                        if keep:
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
                                if database.angsdk[datndx]:
                                    print("{:>6s} {:>6s} {:>6s}    4  {:8.4f} {:8.4f} {:8.4f} {:8.4f}".format(
                                            database.angtype1[ang_params[uniq_angs-1]],
                                            database.angtype2[ang_params[uniq_angs-1]],
                                            database.angtype3[ang_params[uniq_angs-1]],
                                            tmpcalca, tmpcalcb, sig, eps), file=fout)
                                else:
                                    print("{:>6s} {:>6s} {:>6s}    1  {:8.4f} {:8.4f}".format(
                                            database.angtype1[ang_params[uniq_angs-1]],
                                            database.angtype2[ang_params[uniq_angs-1]],
                                            database.angtype3[ang_params[uniq_angs-1]],
                                            tmpcalca, tmpcalcb), file=fout)
        sysdat.uniq_nangs = uniq_angs
        # get dihedral interactions
        print(file=fout)
        if sysdat.ndihs > 0:
            print("[ dihedraltypes ]", file=fout)
            print("; i    j    k    l    func    phi0    kphi    mult", file=fout)
            for idx in range(sysdat.ntops):
                for jdx in range(topdat[idx].ndih):
                    datndx = []
                    if not topdat[idx].dihpset[jdx]:
                        top_atype1 = topdat[idx].atomtype[topdat[idx].dihndx1[jdx]-1]
                        top_atype2 = topdat[idx].atomtype[topdat[idx].dihndx2[jdx]-1]
                        top_atype3 = topdat[idx].atomtype[topdat[idx].dihndx3[jdx]-1]
                        top_atype4 = topdat[idx].atomtype[topdat[idx].dihndx4[jdx]-1]
                        # now compare to the database
                        for kdx in range(database.ndihtype):
                            dat_atype1 = database.dihtype1[kdx]
                            dat_atype2 = database.dihtype2[kdx]
                            dat_atype3 = database.dihtype3[kdx]
                            dat_atype4 = database.dihtype4[kdx]
                            if cmp_wc(dat_atype2, top_atype2):
                                if cmp_wc(dat_atype3, top_atype3):
                                    f1 = cmp_wc(dat_atype1, top_atype1)
                                    f2 = cmp_wc(dat_atype4, top_atype4)
                                    if top_atype2 == top_atype3:
                                        f3 = cmp_wc(dat_atype1, top_atype4)
                                        f4 = cmp_wc(dat_atype4, top_atype1)
                                        if f1 and f2 or f3 and f4:
                                            datndx.append(kdx)
                                    else:
                                        if f1 and f2:
                                            datndx.append(kdx)
                            elif cmp_wc(dat_atype2, top_atype3):
                                if cmp_wc(dat_atype3, top_atype2):
                                    f1 = cmp_wc(dat_atype1, top_atype4)
                                    f2 = cmp_wc(dat_atype4, top_atype1)
                                    if f1 and f2:
                                        datndx.append(kdx)
                        # No params for this interaction in the database
                        if len(datndx) == 0:
                            sys.exit("ERROR: Did not find dihedral parameters in database {} {} {} {} ({} {} {} {})".format(
                                      top_atype1,
                                      top_atype2,
                                      top_atype3,
                                      top_atype4,
                                      topdat[idx].dihndx1[jdx],
                                      topdat[idx].dihndx2[jdx],
                                      topdat[idx].dihndx3[jdx],
                                      topdat[idx].dihndx4[jdx]))
                        # Now make sure we do not already have this one
                        keep = True
                        keeps = []
                        for ldx in datndx:
                            if database.dihe[ldx] != -1:
                                for kdx in range(uniq_dihs):
                                    if ldx == dih_params[kdx]:
                                        keep = False
                                        keeps.append(kdx)
                        if not keep:
                            topdat[idx].dihtype.append([uniq_dihs, keeps])
                        else:
                            dih_params += datndx
                            sysdat.param_dihs += datndx
                            topdat[idx].dihtype.append([uniq_dihs, [uniq_dihs + x for x in range(len(datndx))]])
                            uniq_dihs += 1
                            tmpcalca = database.dihe[dih_params[uniq_dihs-1]]
                            tmpcalcb = database.fdih[dih_params[uniq_dihs-1]]*4.184
                            tmpcalcc = database.dihn[dih_params[uniq_dihs-1]]
                            print("{:>6s} {:>6s} {:>6s} {:>6s}    9  {:8.1f} {:8.4f} {:8}".format(
                                    database.dihtype1[dih_params[uniq_dihs-1]],
                                    database.dihtype2[dih_params[uniq_dihs-1]],
                                    database.dihtype3[dih_params[uniq_dihs-1]],
                                    database.dihtype4[dih_params[uniq_dihs-1]],
                                    tmpcalca, tmpcalcb, int(tmpcalcc)), file=fout)
        sysdat.uniq_ndihs = uniq_dihs
        # get improper interactions
        print(file=fout)
        if sysdat.nimps > 0:
            print("[ dihedraltypes ]", file=fout)
            print("; i    j    k    l    func    phi0    kphi", file=fout)
            for idx in range(sysdat.ntops):
                for jdx in range(topdat[idx].nimp):
                    datndx = -1
                    if not topdat[idx].imppset[jdx]:
                        # now compare to the database
                        for kdx in range(database.nimptype):
                            if cmp_wc(database.imptype2[kdx], topdat[idx].atomtype[topdat[idx].impndx2[jdx]-1]):
                                if cmp_wc(database.imptype3[kdx], topdat[idx].atomtype[topdat[idx].impndx3[jdx]-1]):
                                    f1 = cmp_wc(database.imptype1[kdx], topdat[idx].atomtype[topdat[idx].impndx1[jdx]-1])
                                    f2 = cmp_wc(database.imptype4[kdx], topdat[idx].atomtype[topdat[idx].impndx4[jdx]-1])
                                    if f1 and f2:
                                        datndx = kdx
                                        break
                            if cmp_wc(database.imptype2[kdx], topdat[idx].atomtype[topdat[idx].impndx3[jdx]-1]):
                                if cmp_wc(database.imptype3[kdx], topdat[idx].atomtype[topdat[idx].impndx2[jdx]-1]):
                                    f1 = cmp_wc(database.imptype1[kdx], topdat[idx].atomtype[topdat[idx].impndx4[jdx]-1])
                                    f2 = cmp_wc(database.imptype4[kdx], topdat[idx].atomtype[topdat[idx].impndx1[jdx]-1])
                                    if f1 and f2:
                                        datndx = kdx
                                        break
                        # No params for this interaction in the database
                        if datndx == -1:
                            sys.exit("ERROR: Did not find improper parameters in database {} {} {} ({} {} {})".format(
                                      topdat[idx].atomtype[topdat[idx].impndx1[jdx]-1],
                                      topdat[idx].atomtype[topdat[idx].impndx2[jdx]-1],
                                      topdat[idx].atomtype[topdat[idx].impndx3[jdx]-1],
                                      topdat[idx].atomtype[topdat[idx].impndx3[jdx]-1],
                                      topdat[idx].impndx1[jdx],
                                      topdat[idx].impndx2[jdx],
                                      topdat[idx].impndx3[jdx],
                                      topdat[idx].impndx4[jdx]))
                        # Now make sure we do not already have this one
                        keep = True
                        for kdx in range(uniq_imps):
                            if datndx == imp_params[kdx]:
                                keep = False
                                topdat[idx].imptype.append(kdx)
                                break
                        if keep:
                            imp_params.append(datndx)
                            sysdat.param_imps.append(datndx)
                            topdat[idx].imptype.append(uniq_imps)
                            uniq_imps += 1
                            tmpcalca = database.impe[imp_params[uniq_imps-1]]
                            tmpcalcb = database.fimp[imp_params[uniq_imps-1]]*4.184*2.0
                            print("{:>6s} {:>6s} {:>6s} {:>6s}    2  {:8.4f} {:8.4f}".format(
                                    database.imptype1[imp_params[uniq_imps-1]],
                                    database.imptype2[imp_params[uniq_imps-1]],
                                    database.imptype3[imp_params[uniq_imps-1]],
                                    database.imptype4[imp_params[uniq_imps-1]],
                                    tmpcalca, tmpcalcb), file=fout)
        sysdat.uniq_nimps = uniq_imps
        print(file=fout)


# Read the database file and store unique params
# Warn if you find duplicates
def read_database(fname, database):
    nvdw = nbnd = nang = ndih = nimp = 0
    with open(fname, "r") as fin:
        line = fin.readline()
        while line:
            items = line.split()
            if len(items) == 0:
                line = fin.readline()
                continue
            if items[0] == "pair":
                keep = True
                vdwtype1 = items[1]
                vdwtype2 = items[2]
                vdwstyle = items[3]
                eps = float(items[4])
                sig = float(items[5])
                for idx in range(nvdw):
                    if vdwtype1 == database.vdwtype1[idx] and vdwtype2 == database.vdwtype2[idx]:
                        print("WARNING: Found dup vdw param {} {}".format(vdwtype1,vdwtype2))
                        keep = False
                    elif vdwtype1 == database.vdwtype2[idx] and vdwtype2 == database.vdwtype1[idx]:
                        print("WARNING: Found dup vdw param {} {}".format(vdwtype1,vdwtype2))
                        keep = False
                if keep:
                    database.vdwtype1.append(vdwtype1)
                    database.vdwtype2.append(vdwtype2)
                    database.vdwstyle.append(vdwstyle)
                    database.eps.append(eps)
                    database.sig.append(sig)
                    nvdw += 1
                if vdwtype1 == 'GBML':
                    database.loop_pair.append(vdwtype2)
                elif vdwtype2 == 'GBML':
                    database.loop_pair.append(vdwtype1)
            if items[0] == "bond":
                keep = True
                bndtype1 = items[1]
                bndtype2 = items[2]
                fbnd = float(items[3])
                bnde = float(items[4])
                for idx in range(nbnd):
                    if bndtype1 == database.bndtype1[idx] and bndtype2 == database.bndtype2[idx]:
                        print("WARNING: Found dup bond param {} {}".format(bndtype1,bndtype2))
                        keep = False
                    elif bndtype1 == database.bndtype2[idx] and bndtype2 == database.bndtype1[idx]:
                        print("WARNING: Found dup bond param {} {}".format(bndtype1,bndtype2))
                        keep = False
                if keep:
                    database.bndtype1.append(bndtype1)
                    database.bndtype2.append(bndtype2)
                    database.fbnd.append(fbnd)
                    database.bnde.append(bnde)
                    nbnd += 1
            if items[0] == "angle":
                if "harmonic" in line:
                    angsdk = False
                else:
                    angsdk = True
                keep = True
                angtype1 = items[1]
                angtype2 = items[2]
                angtype3 = items[3]
                fang = float(items[4])
                ange = float(items[5])
                for idx in range(nang):
                    if angtype2 == database.angtype2[idx]:
                        if angtype1 == database.angtype1[idx] and angtype3 == database.angtype3[idx]:
                            print("WARNING: Found dup angle param {} {} {}".format(angtype1,angtype2,angtype3))
                            keep = False
                        elif angtype3 == database.angtype1[idx] and angtype1 == database.angtype3[idx]:
                            print("WARNING: Found dup angle param {} {} {}".format(angtype1,angtype2,angtype3))
                            keep = False
                if keep:
                    database.angtype1.append(angtype1)
                    database.angtype2.append(angtype2)
                    database.angtype3.append(angtype3)
                    database.fang.append(fang)
                    database.ange.append(ange)
                    database.angsdk.append(angsdk)
                    nang += 1
            if items[0] == "dihedral":
                keep = True
                dihtype1 = items[1]
                dihtype2 = items[2]
                dihtype3 = items[3]
                dihtype4 = items[4]
                fdih = float(items[6])
                dihn = float(items[7])
                dihe = float(items[8])
                for idx in range(ndih):
                    if dihn == database.dihn[idx]:
                        if dihtype2 == database.dihtype2[idx] and dihtype3 == database.dihtype3[idx]:
                            if dihtype1 == database.dihtype1[idx] and dihtype4 == database.dihtype4[idx]:
                                print("WARNING: Found dup dihedral param {} {} {} {}".format(dihtype1,dihtype2,dihtype3,dihtype4))
                                keep = False
                        if dihtype2 == database.dihtype3[idx] and dihtype3 == database.dihtype2[idx]:
                            if dihtype1 == database.dihtype4[idx] and dihtype4 == database.dihtype1[idx]:
                                print("WARNING: Found dup dihedral param {} {} {} {}".format(dihtype1,dihtype2,dihtype3,dihtype4))
                                keep = False
                if keep:
                    database.dihtype1.append(dihtype1)
                    database.dihtype2.append(dihtype2)
                    database.dihtype3.append(dihtype3)
                    database.dihtype4.append(dihtype4)
                    database.fdih.append(fdih)
                    database.dihn.append(dihn)
                    database.dihe.append(dihe)
                    ndih += 1
            if items[0] == "improper":
                keep = True
                imptype1 = items[1]
                imptype2 = items[2]
                imptype3 = items[3]
                imptype4 = items[4]
                fimp     = float(items[5])
                impe     = float(items[6])
                for idx in range(nimp):
                    if imptype2 == database.imptype2[idx] and imptype3 == database.imptype3[idx]:
                        if imptype1 == database.imptype1[idx] and imptype4 == database.imptype4[idx]:
                            print("WARNING: Found dup improper param {} {} {} {}".format(imptype1,imptype2,imptype3,imptype4))
                            keep = False
                    if imptype2 == database.imptype3[idx] and imptype3 == database.imptype2[idx]:
                        if imptype1 == database.imptype4[idx] and imptype4 == database.imptype1[idx]:
                            print("WARNING: Found dup improper param {} {} {} {}".format(imptype1,imptype2,imptype3,imptype4))
                            keep = False
                if keep:
                    database.imptype1.append(imptype1)
                    database.imptype2.append(imptype2)
                    database.imptype3.append(imptype3)
                    database.imptype4.append(imptype4)
                    database.fimp.append(fimp)
                    database.impe.append(impe)
                    nimp += 1
            line = fin.readline()
        database.nvdwtype = nvdw
        database.nbndtype = nbnd
        database.nangtype = nang
        database.ndihtype = ndih
        database.nimptype = nimp


# Count the number of params in the database so we can allocate for storage
def count_params(fname, database):
    database.nvdwtype = 0
    database.nbndtype = 0
    database.nangtype = 0
    database.ndihtype = 0
    database.nimptype = 0
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
            if items[0] == "dihedral":
                database.ndihtype += 1
            if items[0] == "improper":
                database.nimptype += 1
            line = fin.readline()
        

# count the number of things in the topology files so we can allocate
def count_atoms(fname, topdat, ntop):
    topdat[ntop].nat = 0
    topdat[ntop].nbnd = 0
    topdat[ntop].nang = 0
    topdat[ntop].ndih = 0
    topdat[ntop].nimp = 0
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
                topdat[ntop].ndih  += 1
            if items[0] == "improper" or items[0] == "improperparam":
                topdat[ntop].nimp += 1
            line = fin.readline()
        if topdat[ntop].nat == 0:
            sys.exit("ERROR: # of atoms in {} is zero.".format(fname))


# Read the topology file and store the data
def read_top(fname, sysdat, topdat, ntop):
    log_bndprm = log_angprm = log_dihprm = log_impprm = log_charge = True
    ndx = bndx = andx = dndx = indx = lc = 0
    print("######################")
    print("##### Reading {}".format(fname))
    with open(fname, "r") as fin:
        line = fin.readline()
        while line:
            lc += 1
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
                    print("NOTE: Using bond parameters from the top file.")
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
                    print("NOTE: Using angle parameters from the top file.")
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
                try:
                    topdat[ntop].impndx1.append(int(items[1]))
                    topdat[ntop].impndx2.append(int(items[2]))
                    topdat[ntop].impndx3.append(int(items[3]))
                    topdat[ntop].impndx4.append(int(items[4]))
                    topdat[ntop].impfk.append(None)
                    topdat[ntop].impeq.append(None)
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].imppset.append(False)
                indx += 1
            if items[0] == "improperparam":
                if log_impprm:
                    print("NOTE: Using improper parameters from the top file.")
                    log_impprm = False
                if len(items) < 7:
                    sys.exit("ERROR: Not enough args for improperparam: must be: ndx1 ndx2 ndx3 ndx4 fk eq.")
                try:
                    topdat[ntop].impndx1.append(int(items[1]))
                    topdat[ntop].impndx2.append(int(items[2]))
                    topdat[ntop].impndx3.append(int(items[3]))
                    topdat[ntop].impndx4.append(int(items[4]))
                    topdat[ntop].impfk.append(float(items[5]))
                    topdat[ntop].impeq.append(float(items[6]))
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].imppset.append(True)
                indx += 1
            if items[0] == "dihedral":
                try:
                    topdat[ntop].dihndx1.append(int(items[1]))
                    topdat[ntop].dihndx2.append(int(items[2]))
                    topdat[ntop].dihndx3.append(int(items[3]))
                    topdat[ntop].dihndx4.append(int(items[4]))
                    topdat[ntop].dihfk.append(None)
                    topdat[ntop].dihn.append(None)
                    topdat[ntop].diheq.append(None)
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].dihpset.append(False)
                dndx += 1
            if items[0] == "dihedralparam":
                if log_dihprm:
                    print("NOTE: Using dihedral parameters from the top file.")
                    log_dihprm = False
                if len(items) < 9:
                    sys.exit("ERROR: Not enough args for angleparam: must be: ndx1 ndx2 ndx3 fk n eq onefour.")
                try:
                    topdat[ntop].dihndx1.append(int(items[1]))
                    topdat[ntop].dihndx2.append(int(items[2]))
                    topdat[ntop].dihndx3.append(int(items[3]))
                    topdat[ntop].dihndx4.append(int(items[4]))
                    topdat[ntop].dihfk.append(float(items[5]))
                    topdat[ntop].dihn.append(int(items[6]))
                    topdat[ntop].diheq.append(int(float(items[7])))
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].dihpset.append(True)
                dndx += 1
            line = fin.readline()


def make_ndx(database, topdat, sysdat):
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
        elif vt2 in sols or vt2 in psolsp:
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


def make_top(topdat, sysdat):
    with open("topol.top", "w") as fout:
        print("; generated by cg_spica setup_gmax", file=fout)
        print("#include \"toppar/SPICA.itp\"", file=fout)
        lst = [ topdat[idx].fname for idx in range(sysdat.ntops) ]
        topnames = sorted(set(lst), key=lst.index)
        for topname in topnames:
            print(f"#include \"toppar/{topname}.itp\"", file=fout)
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
        if nargs < 3 or nargs % 2  == 0:
            print("Dumps input files for a GROMACS run.")
            print("usage: setup_gmx [-p] <topfile 1> <nmol 1> [ <topfile 2> <nmol 2> ..... <topfile n> <nmol n>] <database>")
            print("Takes at least three arguments (one component system): 1) Topology, 2) number of molecules, 3) parameter database.")
            sys.exit(1)
        print("setup_gmx for SPICA.")
        ntops = int((nargs-1)/2)
    topdat = [Topdat() for _ in range(ntops)]
    database = Database()
    sysdat = Sysdat()
    print("Will read {} topology file(s).".format(ntops))
    # Loop through the topologies and count the number of atoms, bonds and bends
    sysdat.ntops = ntops
    sysdat.nats = 0
    sysdat.nbnds = 0
    sysdat.nangs = 0
    sysdat.total_ats = 0
    sysdat.total_bnds = 0
    sysdat.total_angs = 0
    sysdat.total_dihs = 0
    sysdat.total_imps = 0
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
        print("Found: {} dihedrals".format(topdat[idx].ndih))
        print("Found: {} impropers".format(topdat[idx].nimp))
        sysdat.nats	+= topdat[idx].nat
        sysdat.nbnds += topdat[idx].nbnd
        sysdat.nangs += topdat[idx].nang
        sysdat.ndihs += topdat[idx].ndih
        sysdat.nimps += topdat[idx].nimp
        sysdat.total_ats += topdat[idx].nat*topdat[idx].nmol
        sysdat.total_bnds += topdat[idx].nbnd*topdat[idx].nmol
        sysdat.total_angs += topdat[idx].nang*topdat[idx].nmol
        sysdat.total_dihs += topdat[idx].ndih*topdat[idx].nmol
        sysdat.total_imps += topdat[idx].nimp*topdat[idx].nmol
    if args.prot:
        count_params(inputs[nargs-2], database)
        read_database(inputs[nargs-2], database)
    else:
        count_params(inputs[nargs-1], database)
        read_database(inputs[nargs-1], database)
    print("Found {} unique vdw pair params".format(database.nvdwtype))
    print("Found {} unique bond params".format(database.nbndtype))
    print("Found {} unique angle params".format(database.nangtype))
    print("Found {} unique dihedral params".format(database.ndihtype))
    print("Found {} unique improper params".format(database.nimptype))
    if args.prot:
        if ".pdb" in inputs[nargs-1]:
            print("Takes angles from {}".format(inputs[nargs-1]))
            read_pdb(inputs[nargs-1],sysdat)
        else:
            sys.exit("ERROR: No pdb files to take angles!")
    get_unique(database, topdat, sysdat)
    read_coords(database, topdat, sysdat)
    make_top(topdat, sysdat)
    make_ndx(database, topdat, sysdat)
    write_psf(topdat, sysdat)

if __name__ == "__main__":
    args = get_option()
    run(args)
