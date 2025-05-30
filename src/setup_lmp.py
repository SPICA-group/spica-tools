import numpy as np
import sys, textwrap
from math import sqrt
from math import pi
from pathlib import Path
from argparse import ArgumentParser

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument('input_files', type=str, nargs="+",
                            help='<topfile 1> <nmol 1> [ <topfile 2> <nmol 2> ..... <topfile n> <nmol n>] <param file> <coordfile>')
    argparser.add_argument('-Go', action='store_true',help='Go model for protein backbone')
    return argparser.parse_args()


def get_option_script(argv):
    argparser = ArgumentParser(usage='setup_lmp [-h] [-Go] input_files',
                               prog ="setup_lmp")
    argparser.add_argument('input_files', type=str, nargs="+",
                            help='<topfile 1> <nmol 1> [ <topfile 2> <nmol 2> ..... <topfile n> <nmol n>] <param file> <coordfile>')
    argparser.add_argument('-Go', action='store_true',help='Go model for protein backbone')
    return argparser.parse_args(argv)


def get_angle(r1, r2, r3):
    r12 = r1 - r2
    r32 = r3 - r2
    r13_inn = np.dot(r12, r32)
    r13_mag = np.linalg.norm(r12)*np.linalg.norm(r32)
    cos13 = r13_inn/r13_mag
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
    nats = ngo = nbnds = nangs = nimps = ndihs = ntops = 0 
    total_ats = total_go = total_bnds = total_angs = total_imps = total_dihs = 0
    total_all_dihs = 0
    foundatoms = boxinfo = ischarged = 0
    uniq_nats = uniq_nbnds = uniq_nangs = uniq_nimps = uniq_ndihs = 0
    keep_ndihs = 0
    param_bnds, param_angs, param_dihs, param_imps  = [], [], [], []
    coordx, coordy, coordz = [], [], []
    boxx = boxy = boxz = 0.0


class Topdat:
    def __init__(self):
        self.nat = self.nbnd = self.nang = self.ndih = self.nimp = self.nmol = self.ngo = 0
        self.nalldih = 0
        self.gondx1, self.gondx2, self.gofunctype = [], [], []
        self.eps, self.sig = [], []
        self.bndndx1, self.bndndx2, self.bndtype = [], [], []
        self.angndx1, self.angndx2, self.angndx3, self.angtype = [], [], [], []
        self.dihed_func, self.improp_func = [], []
        self.dihndx1, self.dihndx2, self.dihndx3, self.dihndx4 = [], [], [], []
        self.impndx1, self.impndx2, self.impndx3, self.impndx4 = [], [], [], []
        self.imptype, self.dihtype, self.dihedn = [], [], []
        self.bndpset, self.angpset, self.dihpset, self.imppset = [], [], [], []
        self.ind, self.parm_atomtype, self.dihedeq = [], [], []
        self.dihedfk, self.dihedof = [], []
        self.mass, self.charge = [], [] 
        self.bndfk, self.bndeq = [], []
        self.angfk, self.angeq = [], []
        self.impropfk, self.impropeq = [], []
        self.atomname, self.atomtype, self.segid, self.resname = [], [], [], []
   

class Database:
    fbnd, bnde, fang, ange, eps, sig, angsdk = [], [], [], [], [], [], []
    fdih, dihn, dihd, fimp, impe = [], [], [], [], []
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
                print("Found Boxsize data")
                sysdat.boxx = float(items[1])
                sysdat.boxy = float(items[2])
                sysdat.boxz = float(items[3])
            elif items[0] == "ATOM" or items[0] == "HETATM":
                if sysdat.foundatoms >= sysdat.total_ats:
                    sys.exit(f"ERROR: Found atoms in pdb file >= total atoms in top files")
                sysdat.coordx.append(float(line[col:col+8]))
                sysdat.coordy.append(float(line[col+8:col+16]))
                sysdat.coordz.append(float(line[col+16:col+24]))
                sysdat.foundatoms += 1
            line = fin.readline()
        if sysdat.foundatoms == 0:
            sys.exit("ERROR: Did not find any atoms in the pdb file")
        if sysdat.boxx == 0.0:
            print("WARNING: Did not find cell size. Box size must be set by hand")
        
def read_coords(topdat, sysdat):
    if sysdat.ischarged:
        print("Found Charges")
    else:
        print("Did not find charges")
    if sysdat.total_ats != sysdat.foundatoms:
       print("ERROR: # of atoms read from topology and command line")
       sys.exit("Does not match # of atoms found in the coord file")

    with open("DATA.FILE", "w") as fout:
        print("LAMMPS description", file=fout)
        print(file=fout)
        print("{:<7} atoms".format(sysdat.total_ats), file=fout)
        if sysdat.nbnds > 0:
            print("{:<7} bonds".format(sysdat.total_bnds), file=fout)
        if sysdat.nangs  > 0:
            print("{:<7} angles".format(sysdat.total_angs), file=fout)
        if sysdat.total_dihs > 0 or sysdat.total_imps > 0:
            all_dihs = 0
            for idx in range(sysdat.ntops):
                all_dihs += topdat[idx].nalldih*topdat[idx].nmol
            print("{:<7} dihedrals".format(all_dihs + sysdat.total_imps), file=fout)
        #if sysdat.total_imps > 0:
        #    print("{:<7} impropers".format(sysdat.total_imps), file=fout)
        print(file=fout)
        print("{:<7} atom types".format(sysdat.uniq_nats), file=fout)
        if sysdat.nbnds > 0:
            print("{:<7} bond types".format(sysdat.uniq_nbnds), file=fout)
        if sysdat.nangs > 0:
            print("{:<7} angle types".format(sysdat.uniq_nangs), file=fout)
        if sysdat.total_dihs > 0:
            print("{:<7} dihedral types".format(sysdat.uniq_ndihs + sysdat.uniq_nimps), file=fout)
        #if sysdat.total_imps > 0:
        #    print("{:<7} improper types".format(sysdat.uniq_nimps), file=fout)
        print(file=fout)
        print("{:8.3f} {:8.3f} xlo xhi".format(-0.5*sysdat.boxx, 0.5*sysdat.boxx), file=fout)
        print("{:8.3f} {:8.3f} ylo yhi".format(-0.5*sysdat.boxy, 0.5*sysdat.boxy), file=fout)
        print("{:8.3f} {:8.3f} zlo zhi".format(-0.5*sysdat.boxz, 0.5*sysdat.boxz), file=fout)
        print(file=fout)
        print("Atoms", file=fout)
        print(file=fout)
        atidx = molidx = 0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nmol):
                molidx += 1
                for kdx in range(topdat[idx].nat):
                    atidx += 1
                    print("{:<7} {:<7} {:<7} {:7.4f} {:9.4f} {:9.4f} {:9.4f} # {}".format(
                           atidx, molidx, topdat[idx].parm_atomtype[kdx]+1, topdat[idx].charge[kdx], 
                           sysdat.coordx[atidx-1], sysdat.coordy[atidx-1], sysdat.coordz[atidx-1],
                           topdat[idx].atomtype[kdx]), file=fout)
        if sysdat.total_bnds > 0:
            print(file=fout)
            print("Bonds", file=fout)
            print(file=fout)
            bondidx = 0
            offset  = 0
            for idx in range(sysdat.ntops):
                for jdx in range(topdat[idx].nmol):
                    molidx += 1
                    for kdx in range(topdat[idx].nbnd):
                        bondidx += 1
                        print("{:<8} {:<8} {:<8} {:<8} # {:<4} {:<4}".format(
                               bondidx, topdat[idx].bndtype[kdx]+1,
                               topdat[idx].bndndx1[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].bndndx2[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].atomtype[topdat[idx].bndndx1[kdx]-1],
                               topdat[idx].atomtype[topdat[idx].bndndx2[kdx]-1]), file=fout)
                offset += topdat[idx].nmol*topdat[idx].nat
        if sysdat.total_angs > 0:
            print(file=fout);
            print("Angles", file=fout);
            print(file=fout);
            angleidx = 0
            offset   = 0
            for idx in range(sysdat.ntops):
                for jdx in range(topdat[idx].nmol):
                    molidx += 1
                    for kdx in range(topdat[idx].nang):
                        angleidx += 1
                        print("{:<8} {:<8} {:<8} {:<8} {:<8} # {:<4} {:<4} {:<4}".format(
                               angleidx, topdat[idx].angtype[kdx]+1,
                               topdat[idx].angndx1[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].angndx2[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].angndx3[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].atomtype[topdat[idx].angndx1[kdx]-1],
                               topdat[idx].atomtype[topdat[idx].angndx2[kdx]-1],
                               topdat[idx].atomtype[topdat[idx].angndx3[kdx]-1]), file=fout)
                offset += topdat[idx].nmol*topdat[idx].nat
        if sysdat.total_dihs > 0:
            print(file=fout)
            print("Dihedrals", file=fout)
            print(file=fout)
            dihedidx = 0
            offset   = 0
            for idx in range(sysdat.ntops):
                for jdx in range(topdat[idx].nmol):
                    molidx += 1
                    for kdx in range(topdat[idx].ndih):
                        for ldx in range(len(topdat[idx].dihtype[kdx][1])):
                            dihedidx += 1
                            print("{:<8} {:<8} {:<8} {:<8} {:<8} {:<8} # {:<4} {:<4} {:<4} {:<4}".format(
                                   dihedidx, topdat[idx].dihtype[kdx][1][ldx]+1,
                                   topdat[idx].dihndx1[kdx]+(jdx*topdat[idx].nat)+offset,
                                   topdat[idx].dihndx2[kdx]+(jdx*topdat[idx].nat)+offset,
                                   topdat[idx].dihndx3[kdx]+(jdx*topdat[idx].nat)+offset,
                                   topdat[idx].dihndx4[kdx]+(jdx*topdat[idx].nat)+offset,
                                   topdat[idx].atomtype[topdat[idx].dihndx1[kdx]-1],
                                   topdat[idx].atomtype[topdat[idx].dihndx2[kdx]-1], 
                                   topdat[idx].atomtype[topdat[idx].dihndx3[kdx]-1],
                                   topdat[idx].atomtype[topdat[idx].dihndx4[kdx]-1]), file=fout)
                offset += topdat[idx].nmol*topdat[idx].nat
        if sysdat.total_imps > 0:
            #print(file=fout)
            #print("Impropers", file=fout)
            #print(file=fout)
            impropidx = 0
            offset    = 0
            for idx in range(sysdat.ntops):
                for jdx in range(topdat[idx].nmol):
                    molidx += 1
                    for kdx in range(topdat[idx].nimp):
                        impropidx += 1
                        print("{:<8} {:<8} {:<8} {:<8} {:<8} {:<8} # {:<4} {:<4} {:<4} {:<4}".format(
                               impropidx + dihedidx, topdat[idx].imptype[kdx]+1 + sysdat.uniq_ndihs,
                               topdat[idx].impndx1[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].impndx2[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].impndx3[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].impndx4[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].atomtype[topdat[idx].impndx1[kdx]-1],
                               topdat[idx].atomtype[topdat[idx].impndx2[kdx]-1],
                               topdat[idx].atomtype[topdat[idx].impndx3[kdx]-1],
                               topdat[idx].atomtype[topdat[idx].impndx4[kdx]-1]), file=fout)
                offset += topdat[idx].nmol*topdat[idx].nat;

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
        tot_charge = 0.0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nmol):
                molidx += 1
                for kdx in range(topdat[idx].nat):
                    atidx += 1
                    tot_charge += topdat[idx].charge[kdx]
                    print("{:8} {:<4}{:5} {:<4} {:<4} {:<4}  {:9.6f}  {:12.4f}".format(
                           atidx, topdat[idx].resname[kdx], molidx % 10000, 
                           topdat[idx].resname[kdx], topdat[idx].atomname[kdx],topdat[idx].atomtype[kdx], 
                           topdat[idx].charge[kdx], topdat[idx].mass[kdx]), file=fout)
        if abs(tot_charge) < 1e-6:
            print("Total charge is zero")
        else:
            print("WARNING: Total charge is not zero ({})".format(tot_charge))
        print(file=fout)
        if sysdat.total_bnds > 0:
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
        if sysdat.total_angs > 0:
            print("{:8} !NTHETA: angles".format(sysdat.total_angs), file=fout)
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
        if sysdat.total_dihs > 0:
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
        if sysdat.total_imps > 0:
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

def angle_style(database, topdat, sysdat):
    sdk = False
    harmonic = False
    if sysdat.nangs > 0:
        for idx in range(sysdat.ntops):
            if sdk and harmonic:
                break
            for jdx in range(topdat[idx].nang):
                if sdk and harmonic:
                    break
                for kdx in range(database.nangtype):
                    if cmp_wc(database.angtype2[kdx], topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1]):
                        f1 = cmp_wc(database.angtype1[kdx], topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1])
                        f2 = cmp_wc(database.angtype3[kdx], topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1])
                        if f1 and f2:
                            sdk = sdk or database.angsdk[kdx]
                            harmonic = harmonic or (not database.angsdk[kdx])
                            break
                        f1 = cmp_wc(database.angtype3[kdx], topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1])
                        f2 = cmp_wc(database.angtype1[kdx], topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1])
                        if f1 and f2:
                            sdk = sdk or database.angsdk[kdx]
                            harmonic = harmonic or (not database.angsdk[kdx])
                            break
    if sdk and harmonic:
        return "hybrid"
    elif sdk:
        return "sdk"
    else:
        return "harmonic"

        
def get_unique(database, topdat, sysdat):
    uniq_nats = uniq_bnds = uniq_angs = uniq_dihs = uniq_imps = 0
    uniq_atype, uniq_mass, uniq_charge = [], [], []
    bnd_params, ang_params, ang_vdw = [], [], []
    dih_params, imp_params = [], []
    ang_style = angle_style(database, topdat, sysdat)
    keep_dihs = 0
    fout = open("PARM.FILE", "w")
    print("# Generated by setup_lmp of spica_tools", file=fout)
    print(file=fout)
    if sysdat.ischarged:
        print("pair_style      lj/sdk/coul/long        15.0", file=fout)
    else:
        print("pair_style      lj/sdk        15.0", file=fout)
    if sysdat.nbnds > 0:
        print("bond_style      harmonic", file=fout)
    if sysdat.nangs > 0:
        if ang_style == "hybrid":
            print("angle_style     hybrid sdk harmonic", file=fout)
        elif ang_style == "sdk":
            print("angle_style     sdk", file=fout)
        else:
            print("angle_style     harmonic", file=fout)
    if sysdat.total_dihs > 0 and sysdat.total_imps > 0:
            print("dihedral_style  hybrid charmm quadratic", file=fout)
    elif sysdat.total_dihs > 0 and sysdat.total_imps == 0:
        print("dihedral_style  charmm", file=fout)
    elif sysdat.total_dihs == 0 and sysdat.total_imps > 0:
        print("dihedral_style  quadratic", file=fout)
    if sysdat.ischarged:
        print("special_bonds   lj/coul 0.0 0.0 1.0 angle yes", file=fout)
    else:
        print("special_bonds   lj 0.0 0.0 1.0 angle yes", file=fout)
    print(file=fout)
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
    print("Found {} Unique Atoms".format(uniq_nats))
    sysdat.uniq_nats = uniq_nats
    for idx in range(uniq_nats):
        print("mass   {:<5}  {:6.4f} # {}".format(idx+1, uniq_mass[idx],uniq_atype[idx]), file=fout)
    print(file=fout)
    # get pair interactions
    # Go model for protein backbone
    Go_bool = [[False for i in range(uniq_nats)] for j in range(uniq_nats)]
    if sysdat.ngo > 0:
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].ngo):
                Go_parm_atomtype1 = topdat[idx].parm_atomtype[topdat[idx].gondx1[jdx]-1]
                Go_parm_atomtype2 = topdat[idx].parm_atomtype[topdat[idx].gondx2[jdx]-1]
                print("pair_coeff  {:<5} {:<5} {:<6} {:5.4f} {:5.4f} # {:<4} {:<4} Go model".format(
                        Go_parm_atomtype1+1,Go_parm_atomtype2+1,
                        topdat[idx].gofunctype[jdx],topdat[idx].eps[jdx],topdat[idx].sig[jdx],
                        topdat[idx].atomtype[topdat[idx].gondx1[jdx]-1],topdat[idx].atomtype[topdat[idx].gondx2[jdx]-1]), file=fout)
                Go_bool[Go_parm_atomtype1][Go_parm_atomtype2] = True
    bb_sec = ['GBML','GBBL','ABBL','GBTL','ABTL','GBMS','GBBS','ABBS','GBTS','ABTS']
    for idx in range(uniq_nats):
        for jdx in range(idx, uniq_nats):
            if Go_bool[idx][jdx]:
                continue
            if uniq_atype[idx][0:4] in bb_sec and uniq_atype[jdx] in database.loop_pair:
                tmp_type1 = uniq_atype[idx][0:4]
                tmp_type2 = uniq_atype[jdx]
            elif uniq_atype[jdx][0:4] in bb_sec and uniq_atype[idx] in database.loop_pair:
                tmp_type1 = uniq_atype[idx]
                tmp_type2 = uniq_atype[jdx][0:4]
            else:
                if uniq_atype[idx][0:4] in ['GBTP','GBTN','ABTP','ABTN']:
                    tmp_type1 = uniq_atype[idx][0:4]
                elif uniq_atype[idx][0:3] in ['GBM','GBB','GBT','ABB','ABT']:
                    tmp_type1 = uniq_atype[idx][0:3]
                else:
                    tmp_type1 = uniq_atype[idx]
                if uniq_atype[jdx][0:4] in ['GBTP','GBTN','ABTP','ABTN']:
                    tmp_type2 = uniq_atype[jdx][0:4]
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
                print("pair_coeff  {:<5} {:<5} {:<6} {:5.4f} {:5.4f} # {:<4} {:<4}".format(
                       idx+1, jdx+1, database.vdwstyle[vdwtmp],
                       database.eps[vdwtmp],database.sig[vdwtmp],database.vdwtype1[vdwtmp],database.vdwtype2[vdwtmp]),
                       file=fout)
            else:
                print("*********************");
                print("WARNING: No params for VDW interaction between {} and {}".format(uniq_atype[idx],uniq_atype[jdx]))
                print("Will check for wildcards");
    print(file=fout)
    # BONDS
    if sysdat.nbnds > 0:
        index0 = 0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nbnd):
                datndx = -1
                # At this point we will check to see if the params were given 
                # In the top file....if so we will skip a lot of this and add 
                # This as a unique bond....if not we go through the procedure
                if topdat[idx].bndpset[jdx] == False:
                    # now compare to the database 
                    for kdx in range(database.nbndtype):
                        f1 = cmp_wc(database.bndtype1[kdx], topdat[idx].atomtype[topdat[idx].bndndx1[jdx]-1]) 
                        f2 = cmp_wc(database.bndtype2[kdx], topdat[idx].atomtype[topdat[idx].bndndx2[jdx]-1])
                        if f1 and f2:
                            datndx = kdx
                            break
                        f1 = cmp_wc(database.bndtype2[kdx], topdat[idx].atomtype[topdat[idx].bndndx1[jdx]-1]) 
                        f2 = cmp_wc(database.bndtype1[kdx], topdat[idx].atomtype[topdat[idx].bndndx2[jdx]-1])
                        if f1 and f2:
                            datndx = kdx
                            break
                    if datndx == -1:
                        sys.exit("ERROR: Did not find bond parameters in database {} {} {} {}".format(
                                  topdat[idx].bndndx1[jdx],
                                  topdat[idx].bndndx2[jdx],
                                  topdat[idx].atomtype[topdat[idx].bndndx1[jdx]-1],
                                  topdat[idx].atomtype[topdat[idx].bndndx2[jdx]-1]))
                    # Now make sure we do not already know we have this interaction 
                    keep = True
                    if database.bnde[datndx] > 0:
                        for kdx in range(uniq_bnds):
                            if datndx == bnd_params[kdx]:
                                keep = False # found a replica
                                topdat[idx].bndtype.append(kdx)
                                break # kill the for loop
                    # keep = True if we found a new one
                    if keep:
                        bnd_params.append(datndx)
                        sysdat.param_bnds.append(datndx)
                        topdat[idx].bndtype.append(uniq_bnds)
                        uniq_bnds += 1
                        if database.bnde[datndx] > 0:
                            print("bond_coeff  {:<6} {:8.4f} {:8.4f} # {} {}".format(
                                   uniq_bnds, database.fbnd[bnd_params[uniq_bnds-1]],
                                   database.bnde[bnd_params[uniq_bnds-1]],
                                   database.bndtype1[bnd_params[uniq_bnds-1]],
                                   database.bndtype2[bnd_params[uniq_bnds-1]]), file=fout)
                        else:
                            i1 = topdat[idx].bndndx1[jdx]-1 + index0
                            i2 = topdat[idx].bndndx2[jdx]-1 + index0
                            dx = sysdat.coordx[i1] - sysdat.coordx[i2]
                            dy = sysdat.coordy[i1] - sysdat.coordy[i2]
                            dz = sysdat.coordz[i1] - sysdat.coordz[i2]
                            bond_in_pdb = sqrt(dx*dx+dy*dy+dz*dz);
                            print("bond_coeff  {:<6} {:8.4f} {:8.4f} # {} {}".format(
                                   uniq_bnds, database.fbnd[bnd_params[uniq_bnds-1]],
                                   bond_in_pdb,
                                   database.bndtype1[bnd_params[uniq_bnds-1]],
                                   database.bndtype2[bnd_params[uniq_bnds-1]]), file=fout)
                else:
                    # The params were given in the top file so lets add it to the param file
                    topdat[idx].bndtype.append(uniq_bnds)
                    bnd_params.append(-1)
                    uniq_bnds += 1
                    print("bond_coeff  {:<6} {:8.4f} {:8.4f} # {} {} FROM TOP".format(uniq_bnds,
                           topdat[idx].bndfk[jdx],topdat[idx].bndeq[jdx],
                           topdat[idx].atomtype[topdat[idx].bndndx1[jdx]-1], 
                           topdat[idx].atomtype[topdat[idx].bndndx2[jdx]-1]), file=fout)
            index0 += topdat[idx].nat*topdat[idx].nmol
        # Finished looping over top files
    sysdat.uniq_nbnds = uniq_bnds
    print(file=fout)
    # ANGLES
    if sysdat.nangs > 0:
        uniq_angs = index0 = 0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nang):
                if topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1][0:4] in ['GBTP','GBTN','ABTP','ABTN']:
                    tmp_type1 = topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1][0:4]
                elif topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1][0:3] in ['GBM','GBB','GBT','ABB','ABT']:
                    tmp_type1 = topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1][0:3]
                else:
                    tmp_type1 = topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1]
                if topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1][0:4] in ['GBTP','GBTN','ABTP','ABTN']:
                    tmp_type2 = topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1][0:4]
                elif topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1][0:3] in ['GBM','GBB','GBT','ABB','ABT']:
                    tmp_type2 = topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1][0:3]
                else:
                    tmp_type2 = topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1]
                if topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1][0:4] in ['GBTP','GBTN','ABTP','ABTN']:
                    tmp_type3 = topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1][0:4]
                elif topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1][0:3] in ['GBM','GBB','GBT','ABB','ABT']:
                    tmp_type3 = topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1][0:3]
                else:
                    tmp_type3 = topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1]
                datndx = -1
                # At this point we will check to see if the params were given
                # In the top file....if so we will skip a lot of this and add
                # This as a unique bond....if not we go through the procedure
                if topdat[idx].angpset[jdx] == -1:
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
                        f1 = database.vdwtype1[kdx] == tmp_type1
                        f2 = database.vdwtype2[kdx] == tmp_type3
                        f3 = database.vdwtype1[kdx] == tmp_type3
                        f4 = database.vdwtype2[kdx] == tmp_type1
                        if f1 and f2:
                            found = True
                            vdwtmp = kdx
                            break
                        elif f3 and f4:
                            found = True
                            vdwtmp = kdx
                            break
                    if not found:
                        print("*********************")
                        print("ERROR: No params for VDW interaction between {} and {} for angle (database)".format(
                               topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1]))
                        sys.exit("Please update database")
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
                    if database.ange[datndx] > 0:
                        for kdx in range(uniq_angs):
                            if datndx == ang_params[kdx]:
                                keep = False
                                topdat[idx].angtype.append(kdx)
                                break
                    if keep:
                        ang_params.append(datndx)
                        ang_vdw.append(vdwtmp)
                        sysdat.param_angs.append(datndx)
                        topdat[idx].angtype.append(uniq_angs)
                        uniq_angs += 1
                        if database.ange[datndx] > 0:
                            if database.angsdk[ang_params[uniq_angs-1]]:
                                if ang_style == "hybrid":
                                    print("angle_coeff {:<10}      sdk  {:8.4f} {:8.4f} {} {:8.4f} {:8.4f} # {} {} {}".format(
                                           uniq_angs,
                                           database.fang[ang_params[uniq_angs-1]],
                                           database.ange[ang_params[uniq_angs-1]],
                                           database.vdwstyle[ang_vdw[uniq_angs-1]],
                                           database.eps[ang_vdw[uniq_angs-1]],
                                           database.sig[ang_vdw[uniq_angs-1]],
                                           database.angtype1[ang_params[uniq_angs-1]],
                                           database.angtype2[ang_params[uniq_angs-1]],
                                           database.angtype3[ang_params[uniq_angs-1]]), file=fout)
                                else:
                                    print("angle_coeff {:<10} {:8.4f} {:8.4f} {} {:8.4f} {:8.4f} # {} {} {}".format(
                                           uniq_angs,
                                           database.fang[ang_params[uniq_angs-1]],
                                           database.ange[ang_params[uniq_angs-1]],
                                           database.vdwstyle[ang_vdw[uniq_angs-1]],
                                           database.eps[ang_vdw[uniq_angs-1]],
                                           database.sig[ang_vdw[uniq_angs-1]],
                                           database.angtype1[ang_params[uniq_angs-1]],
                                           database.angtype2[ang_params[uniq_angs-1]],
                                           database.angtype3[ang_params[uniq_angs-1]]), file=fout)
                            else:
                                if ang_style == "hybrid":
                                    print("angle_coeff {:<10} harmonic  {:8.4f} {:8.4f} # {} {} {}".format(
                                           uniq_angs,
                                           database.fang[ang_params[uniq_angs-1]],
                                           database.ange[ang_params[uniq_angs-1]],
                                           database.angtype1[ang_params[uniq_angs-1]],
                                           database.angtype2[ang_params[uniq_angs-1]],
                                           database.angtype3[ang_params[uniq_angs-1]]), file=fout);
                                else:
                                    print("angle_coeff {:<10} {:8.4f} {:8.4f} # {} {} {}".format(
                                           uniq_angs,
                                           database.fang[ang_params[uniq_angs-1]],
                                           database.ange[ang_params[uniq_angs-1]],
                                           database.angtype1[ang_params[uniq_angs-1]],
                                           database.angtype2[ang_params[uniq_angs-1]],
                                           database.angtype3[ang_params[uniq_angs-1]]), file=fout);
                        else:
                            i1 = topdat[idx].angndx1[jdx]-1 + index0
                            i2 = topdat[idx].angndx2[jdx]-1 + index0
                            i3 = topdat[idx].angndx3[jdx]-1 + index0
                            r1 = np.array([sysdat.coordx[i1],sysdat.coordy[i1],sysdat.coordz[i1]])
                            r2 = np.array([sysdat.coordx[i2],sysdat.coordy[i2],sysdat.coordz[i2]])
                            r3 = np.array([sysdat.coordx[i3],sysdat.coordy[i3],sysdat.coordz[i3]])
                            angle_in_pdb = 180.0/np.pi*get_angle(r1,r2,r3)
                            if database.angsdk[ang_params[uniq_angs-1]]:
                                if ang_style == "hybrid":
                                    print("angle_coeff {:<10}      sdk  {:8.4f} {:8.4f} {} {:8.4f} {:8.4f} # {} {} {}".format(
                                           uniq_angs,
                                           database.fang[ang_params[uniq_angs-1]],
                                           angle_in_pdb,
                                           database.vdwstyle[ang_vdw[ uniq_angs-1]],
                                           database.eps[ang_vdw[uniq_angs-1]],
                                           database.sig[ang_vdw[uniq_angs-1]],
                                           database.angtype1[ang_params[uniq_angs-1]],
                                           database.angtype2[ang_params[uniq_angs-1]],
                                           database.angtype3[ang_params[uniq_angs-1]]), file=fout)
                                else:
                                    print("angle_coeff {:<10} {:8.4f} {:8.4f} {} {:8.4f} {:8.4f} # {} {} {}".format(
                                           uniq_angs,
                                           database.fang[ang_params[uniq_angs-1]],
                                           angle_in_pdb,
                                           database.vdwstyle[ang_vdw[ uniq_angs-1]],
                                           database.eps[ang_vdw[uniq_angs-1]],
                                           database.sig[ang_vdw[uniq_angs-1]],
                                           database.angtype1[ang_params[uniq_angs-1]],
                                           database.angtype2[ang_params[uniq_angs-1]],
                                           database.angtype3[ang_params[uniq_angs-1]]), file=fout)
                            else:
                                if ang_style == "hybrid":
                                    print("angle_coeff {:<10} harmonic  {:8.4f} {:8.4f} # {} {} {}".format(
                                           uniq_angs,
                                           database.fang[ang_params[uniq_angs-1]],
                                           angle_in_pdb,
                                           database.angtype1[ang_params[uniq_angs-1]],
                                           database.angtype2[ang_params[uniq_angs-1]],
                                           database.angtype3[ang_params[uniq_angs-1]]), file=fout)
                                else:
                                    print("angle_coeff {:<10} {:8.4f} {:8.4f} # {} {} {}".format(
                                           uniq_angs,
                                           database.fang[ang_params[uniq_angs-1]],
                                           angle_in_pdb,
                                           database.angtype1[ang_params[uniq_angs-1]],
                                           database.angtype2[ang_params[uniq_angs-1]],
                                           database.angtype3[ang_params[uniq_angs-1]]), file=fout)
                else:
                    # This param was set in the top file
                    # Still need vdw stuff
                    # Get the VDW for the CG angles
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
                    found = False
                    for kdx in range(database.nvdwtype):
                        f1 = database.vdwtype1[kdx] == tmp_type1 
                        f2 = database.vdwtype2[kdx] == tmp_type3
                        f3 = database.vdwtype1[kdx] == tmp_type3 
                        f4 = database.vdwtype2[kdx] == tmp_type1
                        if f1 and f2:
                            found = True
                            vdwtmp = kdx
                            break
                        elif f3 and f4:
                            found = True
                            vdwtmp = kdx
                            break
                    if not found:
                        print("*********************")
                        print("ERROR: No params for VDW interaction between {} and {} for angle (topfile)".format(
                               topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1]))
                        sys.exit("Please update database")
                     #end VDW for CG angles
                    topdat[idx].angtype.append(uniq_angs)
                    uniq_angs += 1
                    ang_params.append(-1)
                    ang_vdw.append(vdwtmp)
                    if ang_style == "hybrid":
                        print("angle_coeff {:<10} harmonic  {:8.4f} {:8.4f} # {} {} {} FROM TOP".format(uniq_angs,
                        topdat[idx].angfk[jdx],
                        topdat[idx].angeq[jdx],
                        tmp_type1,tmp_type2,tmp_type3), file=fout)
                    else:
                        print("angle_coeff {:<10} {:8.4f} {:8.4f} # {} {} {} FROM TOP".format(uniq_angs,
                        topdat[idx].angfk[jdx],
                        topdat[idx].angeq[jdx],
                        tmp_type1,tmp_type2,tmp_type3), file=fout)

            index0 += topdat[idx].nat*topdat[idx].nmol;
    sysdat.uniq_nangs = uniq_angs
    print(file=fout)
    # DIHEDRAL
    if sysdat.ndihs > 0:
        uniq_dihs = index0 = 0
        for idx in range(sysdat.ntops):
            top_ndih = 0
            for jdx in range(topdat[idx].ndih):
                datndx = []
                if topdat[idx].dihpset[jdx] != 1:
                    # now compare to the database 
                    top_atype1 = topdat[idx].atomtype[topdat[idx].dihndx1[jdx]-1]
                    top_atype2 = topdat[idx].atomtype[topdat[idx].dihndx2[jdx]-1]
                    top_atype3 = topdat[idx].atomtype[topdat[idx].dihndx3[jdx]-1]
                    top_atype4 = topdat[idx].atomtype[topdat[idx].dihndx4[jdx]-1]
                    for kdx in range(database.ndihtype):
                        dat_dtype1 = database.dihtype1[kdx]
                        dat_dtype2 = database.dihtype2[kdx]
                        dat_dtype3 = database.dihtype3[kdx]
                        dat_dtype4 = database.dihtype4[kdx]
                        if cmp_wc(dat_dtype2, top_atype2):
                            if cmp_wc(dat_dtype3, top_atype3):
                                f1 = cmp_wc(dat_dtype1, top_atype1)
                                f2 = cmp_wc(dat_dtype4, top_atype4) 
                                if top_atype2 == top_atype3:
                                    f3 = cmp_wc(dat_dtype1, top_atype4)
                                    f4 = cmp_wc(dat_dtype4, top_atype1)
                                    if f1 and f2 or f3 and f4:
                                        datndx.append(kdx)
                                else:
                                    if f1 and f2:
                                        datndx.append(kdx)
                        elif cmp_wc(dat_dtype2, top_atype3):
                            if cmp_wc(dat_dtype3, top_atype2):
                                f1 = cmp_wc(dat_dtype1, top_atype4) 
                                f2 = cmp_wc(dat_dtype4, top_atype1)
                                if f1 and f2:
                                    datndx.append(kdx)
                    if len(datndx) == 0:
                            sys.exit("ERROR: Did not find dihedral parameters in database {} {} {} {} ({} {} {} {})".format(
                                    topdat[idx].dihndx1[jdx],
                                    topdat[idx].dihndx2[jdx],
                                    topdat[idx].dihndx3[jdx],
                                    topdat[idx].dihndx4[jdx],
                                    top_atype1,
                                    top_atype2,
                                    top_atype3,
                                    top_atype4))
                    # Now make sure we do not already know we have this interaction 
                    keep = True
                    keeps = []
                    for ldx in datndx:
                        if database.dihd[ldx] != -1:
                            for kdx in range(uniq_dihs):
                                if ldx ==  dih_params[kdx]:
                                    keep = False # found a replica
                                    keeps.append(kdx)
                                    keep_dihs +=1
                                    top_ndih +=1
                    if not keep: 
                        topdat[idx].dihtype.append([uniq_dihs, keeps])
                    else:
                        dih_params += datndx
                        sysdat.param_dihs += datndx
                        topdat[idx].dihtype.append([uniq_dihs, [uniq_dihs + x for x in range(len(datndx))]])
                        diht_lst0 = ["dum" for _ in range(4)]
                        wf = 0.0
                        for kdx in range(len(datndx)):
                            uniq_dihs += 1
                            top_ndih += 1
                            if database.dihd[datndx[kdx]] != -1:
                                diht1 = database.dihtype1[dih_params[uniq_dihs-1]]
                                diht2 = database.dihtype2[dih_params[uniq_dihs-1]]
                                diht3 = database.dihtype3[dih_params[uniq_dihs-1]]
                                diht4 = database.dihtype4[dih_params[uniq_dihs-1]]
                                diht_lst = [diht1, diht2, diht3, diht4]
                                #if all(x == y and type(x) == type(y) for x, y in zip(diht_lst, diht_lst0)):
                                #    wf = 0.0
                                #else:
                                #    wf = 1.0
                                if sysdat.nimps > 0:
                                    print("dihedral_coeff {:<8}    charmm {:8.4f} {:2} {:4} {:2.1f} # {} {} {} {}".format(
                                           uniq_dihs,
                                           database.fdih[dih_params[uniq_dihs-1]],
                                           database.dihn[dih_params[uniq_dihs-1]],
                                           int(database.dihd[dih_params[uniq_dihs-1]]),
                                           wf,
                                           diht1,
                                           diht2,
                                           diht3,
                                           diht4), file=fout)
                                else:
                                    print("dihedral_coeff {:<8} {:8.4f} {:2} {:4} {:2.1f} # {} {} {} {}".format(
                                           uniq_dihs,
                                           database.fdih[dih_params[uniq_dihs-1]],
                                           database.dihn[dih_params[uniq_dihs-1]],
                                           int(database.dihd[dih_params[uniq_dihs-1]]),
                                           wf,
                                           diht1,
                                           diht2,
                                           diht3,
                                           diht4), file=fout)
                                diht_lst0 = diht_lst
                            else:
                                i1 = topdat[idx].dihndx1[jdx]-1 + index0;
                                i2 = topdat[idx].dihndx2[jdx]-1 + index0;
                                i3 = topdat[idx].dihndx3[jdx]-1 + index0;
                                i4 = topdat[idx].dihndx4[jdx]-1 + index0;
                                r1 = np.array([sysdat.coordx[i1], sysdat.coordy[i1], sysdat.coordz[i1]])
                                r2 = np.array([sysdat.coordx[i2], sysdat.coordy[i2], sysdat.coordz[i2]])
                                r3 = np.array([sysdat.coordx[i3], sysdat.coordy[i3], sysdat.coordz[i3]])
                                r4 = np.array([sysdat.coordx[i4], sysdat.coordy[i4], sysdat.coordz[i4]])
                                dihedral_in_pdb = 180.0 + 180.0/pi*get_dihedral(r1,r2,r3,r4);
                                if dihedral_in_pdb > 0.0:
                                    idihedral_in_pdb = int(dihedral_in_pdb + 0.5)
                                else:
                                    idihedral_in_pdb = int(dihedral_in_pdb - 0.5)
                                print("dihedral_coeff {:<8} {:8.4f} {:3} {:5} {:2.1f} # {} {} {} {}".format(
                                       uniq_dihs,
                                       database.fdih[dih_params[uniq_dihs-1]],
                                       database.dihn[dih_params[uniq_dihs-1]],
                                       int(database.dihd[dih_params[uniq_dihs-1]]),
                                       1.0,
                                       database.dihtype1[dih_params[uniq_dihs-1]],
                                       database.dihtype2[dih_params[uniq_dihs-1]],
                                       database.dihtype3[dih_params[uniq_dihs-1]],
                                       database.dihtype4[dih_params[uniq_dihs-1]]), file=fout)
                else:
                    # The params were given in the top file so lets add it to the param file
                    topdat[idx].dihtype.append([uniq_dihs, [uniq_dihs]]) 
                    uniq_dihs += 1
                    dih_params += [-1]
                    top_ndih += 1
                    # When the dihedral phase (dihed) is set to -1, the phase parameter is taken from PDB
                    # This option will be implemeted someday.
                    # i1 = topdat[idx].dihndx1[jdx]-1 + index0;
                    # i2 = topdat[idx].dihndx2[jdx]-1 + index0;
                    # i3 = topdat[idx].dihndx3[jdx]-1 + index0;
                    # i4 = topdat[idx].dihndx4[jdx]-1 + index0;
                    # r1 = np.array([sysdat.coordx[i1], sysdat.coordy[i1], sysdat.coordz[i1]])
                    # r2 = np.array([sysdat.coordx[i2], sysdat.coordy[i2], sysdat.coordz[i2]])
                    # r3 = np.array([sysdat.coordx[i3], sysdat.coordy[i3], sysdat.coordz[i3]])
                    # r4 = np.array([sysdat.coordx[i4], sysdat.coordy[i4], sysdat.coordz[i4]])
                    # dihedral_in_pdb = 180.0 + 180.0/pi*get_dihedral(r1,r2,r3,r4);
                    # if dihedral_in_pdb > 0.0:
                    #     idihedral_in_pdb = int(dihedral_in_pdb + 0.5)
                    # else:
                    #     idihedral_in_pdb = int(dihedral_in_pdb - 0.5)
                    if sysdat.nimps > 0:
                        print("dihedral_coeff {:<8}    charmm {:8.4f} {:2} {:4} {:2.1f} # {} {} {} {} FROM TOP".format(
                               uniq_dihs,
                               topdat[idx].dihedfk[jdx],
                               topdat[idx].dihedn[jdx],
                               int(topdat[idx].dihedeq[jdx]),
                               topdat[idx].dihedof[jdx],
                               topdat[idx].atomtype[topdat[idx].dihndx1[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].dihndx2[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].dihndx3[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].dihndx4[jdx]-1]), file=fout)
                    else:
                        print("dihedral_coeff {:<8} {:8.4f} {:3} {:5} {:2.1f} # {} {} {} {} FROM TOP".format(
                               uniq_dihs,
                               topdat[idx].dihedfk[jdx],
                               topdat[idx].dihedn[jdx],
                               int(topdat[idx].dihedeq[jdx]),
                               topdat[idx].dihedof[jdx],
                               topdat[idx].atomtype[topdat[idx].dihndx1[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].dihndx2[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].dihndx3[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].dihndx4[jdx]-1]), file=fout)
            index0 += topdat[idx].nat*topdat[idx].nmol
            topdat[idx].nalldih = top_ndih
    # Finished looping over top files
    sysdat.uniq_ndihs = uniq_dihs
    sysdat.keep_ndihs = keep_dihs
    # IMPROPER
    if sysdat.nimps > 0:
        uniq_imps = index0 = 0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nimp):
                datndx = -1
                if topdat[idx].imppset[jdx] != 1:
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
                    if datndx == -1:
                        sys.exit("ERROR: Did not find improper parameters in database {} {} {} {} ({} {} {} {})".format(
                                topdat[idx].impndx1[jdx],
                                topdat[idx].impndx2[jdx],
                                topdat[idx].impndx3[jdx],
                                topdat[idx].impndx4[jdx],
                                topdat[idx].atomtype[topdat[idx].impndx1[jdx]-1],
                                topdat[idx].atomtype[topdat[idx].impndx2[jdx]-1],
                                topdat[idx].atomtype[topdat[idx].impndx3[jdx]-1],
                                topdat[idx].atomtype[topdat[idx].impndx4[jdx]-1]))
                    # Now make sure we do not already know we have this interaction 
                    keep = True
                    for kdx in range(uniq_imps):
                        if datndx == imp_params[kdx]:
                            keep = False # found a replica
                            topdat[idx].imptype.append(kdx)
                            break # kill the for loop
                    if keep:
                        imp_params.append(datndx)
                        sysdat.param_imps.append(datndx)
                        topdat[idx].imptype.append(uniq_imps)
                        uniq_imps += 1
                        if sysdat.ndihs > 0:
                            print("dihedral_coeff {:<8} quadratic {:8.4f} {:8.4f} # {} {} {} {}".format(
                                   uniq_imps + uniq_dihs,
                                   database.fimp[imp_params[uniq_imps-1]],
                                   database.impe[imp_params[uniq_imps-1]],
                                   database.imptype1[imp_params[uniq_imps-1]],
                                   database.imptype2[imp_params[uniq_imps-1]],
                                   database.imptype3[imp_params[uniq_imps-1]],
                                   database.imptype4[imp_params[uniq_imps-1]]), file=fout)
                        else:
                            print("dihedral_coeff {:<8} {:8.4f} {:8.4f} # {} {} {} {}".format(
                                   uniq_imps + uniq_dihs,
                                   database.fimp[imp_params[uniq_imps-1]],
                                   database.impe[imp_params[uniq_imps-1]],
                                   database.imptype1[imp_params[uniq_imps-1]],
                                   database.imptype2[imp_params[uniq_imps-1]],
                                   database.imptype3[imp_params[uniq_imps-1]],
                                   database.imptype4[imp_params[uniq_imps-1]]), file=fout)
                else:
                    # The params were given in the top file so lets add it to the param file
                    topdat[idx].imptype.append(uniq_imps)
                    uniq_imps += 1
                    imp_params.append(-1)
                    if sysdat.ndihs > 0:
                        print("dihedral_coeff {:<8} quadratic {:8.4f} {:8.4f} # {} {} {} {}".format(
                               uniq_imps + uniq_dihs,
                               topdat[idx].impropfk[jdx],
                               topdat[idx].impropeq[jdx],
                               topdat[idx].atomtype[topdat[idx].impndx1[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].impndx2[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].impndx3[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].impndx4[jdx]-1]), file=fout)
                    else:
                        print("dihedral_coeff {:<8} {:8.4f} {:8.4f} # {} {} {} {} FROM TOP".format(
                               uniq_imps + uniq_dihs,
                               topdat[idx].impropfk[jdx],
                               topdat[idx].impropeq[jdx],
                               topdat[idx].atomtype[topdat[idx].impndx1[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].impndx2[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].impndx3[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].impndx4[jdx]-1]), file=fout)
    sysdat.uniq_nimps = uniq_imps
    fout.close()

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
                        print("WARNING: Found dup vdw pair param {} {}".format(vdwtype1,vdwtype2))
                        keep = False
                    elif vdwtype1 == database.vdwtype2[idx] and vdwtype2 == database.vdwtype1[idx]:
                        print("WARNING: Found dup vdw pair param {} {}".format(vdwtype1,vdwtype2))
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
                dihn = int(float(items[7]))
                dihd = float(items[8])
                for idx in range(ndih):
                    if dihtype2 == database.dihtype2[idx] and dihtype3 == database.dihtype3[idx]:
                        if dihtype1 == database.dihtype1[idx] and dihtype4 == database.dihtype4[idx]:
                            if dihn == database.dihn[idx] and dihd == database.dihd[idx]:
                                print("WARNING: Found dup dihedral param {} {} {} {}".format(dihtype1,dihtype2,dihtype3, dihtype4))
                                keep = False
                    if dihtype2 == database.dihtype3[idx] and dihtype3 == database.dihtype2[idx]:
                        if dihtype1 == database.dihtype4[idx] and dihtype4 == database.dihtype1[idx]:
                            if dihn == database.dihn[idx] and dihd == database.dihd[idx]:
                                print("WARNING: Found dup dihedral param {} {} {} {}".format(dihtype1,dihtype2,dihtype3, dihtype4))
                                keep = False
                if keep:
                    database.dihtype1.append(dihtype1)
                    database.dihtype2.append(dihtype2)
                    database.dihtype3.append(dihtype3)
                    database.dihtype4.append(dihtype4)
                    database.fdih.append(fdih)
                    database.dihn.append(dihn)
                    database.dihd.append(dihd)
                    ndih += 1
            if items[0] == "improper":
                keep = True
                imptype1 = items[1]
                imptype2 = items[2]
                imptype3 = items[3]
                imptype4 = items[4]
                fimp = float(items[5])
                impe = float(items[6])
                for idx in range(nimp):
                    if imptype2 == database.imptype2[idx] and imptype3 == database.imptype3[idx]:
                        if imptype1 == database.imptype1[idx] and imptype4 == database.imptype4[idx]:
                            print("WARNING: Found dup improper param {} {} {} {}".format(imptype1,imptype2,imptype3, imptype4))
                            keep = False
                    if imptype2 == database.imptype3[idx] and imptype3 == database.imptype2[idx]:
                        if imptype1 == database.imptype4[idx] and imptype4 == database.imptype1[idx]:
                            print("WARNING: Found dup improper param {} {} {} {}".format(imptype1,imptype2,imptype3, imptype4))
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
        
# count the number of things in the topology files so we can allocate
def count_atoms(fname, topdat, ntop):
    topdat[ntop].nat = 0
    topdat[ntop].ngo = 0
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
            sys.exit("ERROR: natom in {} is zero.".format(fname))

# Read the topology file and store the data
def read_top(fname, sysdat, topdat, ntop):
    log_bndprm = log_angprm = log_dihprm = log_impprm = log_charge = True
    ndx = bndx = andx = dndx = indx = lc = 0
    print("######################")
    print("##### Reading {}".format(fname))
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
                topdat[ntop].angpset.append(-1)
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
                topdat[ntop].angpset.append(1)
                andx += 1
            if items[0] == "dihedral":
                try:
                    topdat[ntop].dihndx1.append(int(items[1]))
                    topdat[ntop].dihndx2.append(int(items[2]))
                    topdat[ntop].dihndx3.append(int(items[3]))
                    topdat[ntop].dihndx4.append(int(items[4]))
                    topdat[ntop].dihedfk.append(None)
                    topdat[ntop].dihedn.append(None)
                    topdat[ntop].dihedeq.append(None)
                    topdat[ntop].dihedof.append(None)
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].dihpset.append(-1)
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
                    topdat[ntop].dihedfk.append(float(items[5]))
                    topdat[ntop].dihedn.append(int(items[6]))
                    topdat[ntop].dihedeq.append(int(float(items[7])))
                    topdat[ntop].dihedof.append(float(items[8]))
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].dihpset.append(1)
                dndx += 1
            if items[0] == "improper":
                try:
                    topdat[ntop].impndx1.append(int(items[1]))
                    topdat[ntop].impndx2.append(int(items[2]))
                    topdat[ntop].impndx3.append(int(items[3]))
                    topdat[ntop].impndx4.append(int(items[4]))
                    topdat[ntop].impropfk.append(None)
                    topdat[ntop].impropeq.append(None)
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].imppset.append(-1)
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
                    topdat[ntop].impropfk.append(float(items[5]))
                    topdat[ntop].impropeq.append(float(items[6]))
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].imppset.append(1)
                indx += 1
            line = fin.readline()

def read_top_Go(fname, sysdat, topdat, ntop, ndup, bbind):
    log_bndprm = log_angprm = log_dihprm = log_impprm = log_charge = True
    ndx = bndx = andx = dndx = indx = 0
    print("######################")
    print("##### Reading {}".format(fname))
    # Count Number of GBM GBB GBT ABB ABT GBTP GBTN ABTP ABTN 
    nbb = {'GBM':0,'GBB':0,'GBT':0,'ABB':0,'ABT':0,'GBML':0,'GBBL':0,'ABBL':0,'GBTL':0,'ABTL':0,
           'GBMS':0,'GBBS':0,'ABBS':0,'GBTS':0,'ABTS':0,'GBTP':0,'GBTN':0,'ABTP':0,'ABTN':0}
    with open(fname, "r") as fin:
        lines = fin.readlines()
        for line in lines:
            items = line.split()
            if len(items) == 0:
                continue
            if items[0] == "atom":
                if items[4] in nbb:
                    nbb[items[4]] += 1
        for i in range(len(lines)):
            items = lines[i].split()
            if len(items) == 0:
                continue
            if items[0] == "atom":
                try:
                    topdat[ntop].ind.append(int(items[1]))
                    topdat[ntop].resname.append(items[2])
                    topdat[ntop].atomname.append(items[3])
                    if items[4] in bbind:
                        bbind[items[4]] += 1
                        topdat[ntop].atomtype.append(items[4]+str(bbind[items[4]]))
                    else:
                        topdat[ntop].atomtype.append(items[4])
                    topdat[ntop].mass.append(float(items[5]))
                    topdat[ntop].charge.append(float(items[6].replace("+","")))
                    topdat[ntop].segid.append(items[7])
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, i+1))
                for i in range(1,ndup):
                    topdat[ntop+i].ind.append(int(items[1]))
                    topdat[ntop+i].resname.append(items[2])
                    topdat[ntop+i].atomname.append(items[3])
                    if items[4] in nbb:
                        topdat[ntop+i].atomtype.append(items[4]+str(bbind[items[4]]+i*nbb[items[4]]))
                    else:
                        topdat[ntop+i].atomtype.append(items[4])
                    topdat[ntop+i].mass.append(float(items[5]))
                    topdat[ntop+i].charge.append(float(items[6].replace("+","")))
                    topdat[ntop+i].segid.append(items[7])
                if topdat[ntop].charge[ndx]*topdat[ntop].charge[ndx] > 1e-5 and log_charge:
                    print("Charge in top file {} {}".format(fname, topdat[ntop].charge[ndx]))
                    sysdat.ischarged = 1
                    log_charge = False
                ndx += 1
            if items[0] == "goparam":
                try:
                    topdat[ntop].gondx1.append(int(items[1]))
                    topdat[ntop].gondx2.append(int(items[2]))
                    topdat[ntop].gofunctype.append((items[3]))
                    topdat[ntop].eps.append((float(items[4])))
                    topdat[ntop].sig.append((float(items[5])))
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, i+1))
                for i in range(1,ndup):
                    topdat[ntop+i].gondx1.append(int(items[1]))
                    topdat[ntop+i].gondx2.append(int(items[2]))
                    topdat[ntop+i].gofunctype.append((items[3]))
                    topdat[ntop+i].eps.append((float(items[4])))
                    topdat[ntop+i].sig.append((float(items[5])))

            if items[0] == "bond":
                try:
                    topdat[ntop].bndndx1.append(int(items[1]))
                    topdat[ntop].bndndx2.append(int(items[2]))
                    topdat[ntop].bndfk.append(None)
                    topdat[ntop].bndeq.append(None)
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, i+1))
                topdat[ntop].bndpset.append(False)
                for i in range(1,ndup):
                    topdat[ntop+i].bndndx1.append(int(items[1]))
                    topdat[ntop+i].bndndx2.append(int(items[2]))
                    topdat[ntop+i].bndfk.append(None)
                    topdat[ntop+i].bndeq.append(None)
                    topdat[ntop+i].bndpset.append(False)
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
                    sys.exit("ERROR: File {}, line {}".format(fname, i+1))
                topdat[ntop].bndpset.append(True)
                for i in range(1,ndup):
                    topdat[ntop+i].bndndx1.append(int(items[1]))
                    topdat[ntop+i].bndndx2.append(int(items[2]))
                    topdat[ntop+i].bndfk.append(float(items[3]))
                    topdat[ntop+i].bndeq.append(float(items[4]))
                    topdat[ntop+i].bndpset.append(True)
                bndx += 1
            if items[0] == "angle":
                try:
                    topdat[ntop].angndx1.append(int(items[1]))
                    topdat[ntop].angndx2.append(int(items[2]))
                    topdat[ntop].angndx3.append(int(items[3]))
                    topdat[ntop].angfk.append(None)
                    topdat[ntop].angeq.append(None)
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, i+1))
                topdat[ntop].angpset.append(-1)
                for i in range(1,ndup):
                    topdat[ntop+i].angndx1.append(int(items[1]))
                    topdat[ntop+i].angndx2.append(int(items[2]))
                    topdat[ntop+i].angndx3.append(int(items[3]))
                    topdat[ntop+i].angfk.append(None)
                    topdat[ntop+i].angeq.append(None)
                    topdat[ntop+i].angpset.append(-1)
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
                    sys.exit("ERROR: File {}, line {}".format(fname, i+1))
                topdat[ntop].angpset.append(1)
                for i in range(1,ndup):
                    topdat[ntop+i].angndx1.append(int(items[1]))
                    topdat[ntop+i].angndx2.append(int(items[2]))
                    topdat[ntop+i].angndx3.append(int(items[3]))
                    topdat[ntop+i].angfk.append(float(items[4]))
                    topdat[ntop+i].angeq.append(float(items[5]))
                    topdat[ntop+i].angpset.append(1)
                andx += 1
            if items[0] == "dihedral":
                try:
                    topdat[ntop].dihndx1.append(int(items[1]))
                    topdat[ntop].dihndx2.append(int(items[2]))
                    topdat[ntop].dihndx3.append(int(items[3]))
                    topdat[ntop].dihndx4.append(int(items[4]))
                    topdat[ntop].dihedfk.append(None)
                    topdat[ntop].dihedn.append(None)
                    topdat[ntop].dihedeq.append(None)
                    topdat[ntop].dihedof.append(None)
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, lc))
                topdat[ntop].dihpset.append(-1)
                for i in range(1,ndup):
                    topdat[ntop+i].dihndx1.append(int(items[1]))
                    topdat[ntop+i].dihndx2.append(int(items[2]))
                    topdat[ntop+i].dihndx3.append(int(items[3]))
                    topdat[ntop+i].dihndx4.append(int(items[4]))
                    topdat[ntop+i].dihedfk.append(None)
                    topdat[ntop+i].dihedn.append(None)
                    topdat[ntop+i].dihedeq.append(None)
                    topdat[ntop+i].dihedof.append(None)
                    topdat[ntop+i].dihpset.append(-1)
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
                    topdat[ntop].dihedfk.append(float(items[5]))
                    topdat[ntop].dihedn.append(int(items[6]))
                    topdat[ntop].dihedeq.append(int(float(items[7])))
                    topdat[ntop].dihedof.append(float(items[8]))
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, i+1))
                topdat[ntop].dihpset.append(1)
                for i in range(1,ndup):
                    topdat[ntop+i].dihndx1.append(int(items[1]))
                    topdat[ntop+i].dihndx2.append(int(items[2]))
                    topdat[ntop+i].dihndx3.append(int(items[3]))
                    topdat[ntop+i].dihndx4.append(int(items[4]))
                    topdat[ntop+i].dihedfk.append(float(items[5]))
                    topdat[ntop+i].dihedn.append(int(items[6]))
                    topdat[ntop+i].dihedeq.append(int(float(items[7])))
                    topdat[ntop+i].dihedof.append(float(items[8]))
                    topdat[ntop+i].dihpset.append(1)
                dndx += 1
            if items[0] == "improper":
                try:
                    topdat[ntop].impndx1.append(int(items[1]))
                    topdat[ntop].impndx2.append(int(items[2]))
                    topdat[ntop].impndx3.append(int(items[3]))
                    topdat[ntop].impndx4.append(int(items[4]))
                    topdat[ntop].impropfk.append(None)
                    topdat[ntop].impropeq.append(None)
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, i+1))
                topdat[ntop].imppset.append(-1)
                for i in range(1,ndup):
                    topdat[ntop+i].impndx1.append(int(items[1]))
                    topdat[ntop+i].impndx2.append(int(items[2]))
                    topdat[ntop+i].impndx3.append(int(items[3]))
                    topdat[ntop+i].impndx4.append(int(items[4]))
                    topdat[ntop+i].impropfk.append(None)
                    topdat[ntop+i].impropeq.append(None)
                    topdat[ntop+i].imppset.append(-1)
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
                    topdat[ntop].impropfk.append(float(items[5]))
                    topdat[ntop].impropeq.append(float(items[6]))
                except:
                    sys.exit("ERROR: File {}, line {}".format(fname, i+1))
                topdat[ntop].imppset.append(1)
                for i in range(1,ndup):
                    topdat[ntop+i].impndx1.append(int(items[1]))
                    topdat[ntop+i].impndx2.append(int(items[2]))
                    topdat[ntop+i].impndx3.append(int(items[3]))
                    topdat[ntop+i].impndx4.append(int(items[4]))
                    topdat[ntop+i].impropfk.append(float(items[5]))
                    topdat[ntop+i].impropeq.append(float(items[6]))
                    topdat[ntop+i].imppset.append(1)
                indx += 1
    for bb in nbb:
        bbind[bb] += nbb[bb]*(ndup-1)


# Main routine. Call and allocate                                       
# The idea is to read in the topologies and then check the database for 
# all of the required interaction params.                               
def run(args):
    inputs = args.input_files
    Go = args.Go
    nargs = len(inputs)
    if nargs < 4:
       print("usage: setup_lmp [-Go] <topfile 1> <nmol 1> [ <topfile 2> <nmol 2> ..... <topfile n> <nmol n>] <paramfile> <coordfile>");
       print("Prints out input files for a lammps run. Takes a pdb file as the coordfile");
       sys.exit(1)
    ntops = int((nargs - 2)/2)
    topdat = [Topdat() for _ in range(ntops)]
    database = Database() 
    sysdat = Sysdat()
    if Go:
        bbind = {'GBM':0,'GBB':0,'GBT':0,'ABB':0,'ABT':0,'GBML':0,'GBBL':0,'ABBL':0,'GBTL':0,'ABTL':0,
                 'GBMS':0,'GBBS':0,'ABBS':0,'GBTS':0,'ABTS':0,'GBTP':0,'GBTN':0,'ABTP':0,'ABTN':0}
        tmp_ntops = int((nargs-2)/2)
        print("Will read {} topology file(s).".format(tmp_ntops))
        print()
        rdtp = 0
        ndup = []
        protein = []
        for i in range(tmp_ntops):
            f = open(inputs[2*i], 'r')
            line = f.readline() 
            while line:
                if line.split()[0] == "atom" and line.split()[3] in ['GBT','ABT','GBTS','ABTS','GBTL','ABTL','GBTP','ABTP']:
                    protein.append(True)
                    break
                elif line.split()[0] == "atom":
                    protein.append(False)
                    break
                else:
                    line = f.readline()
            f.close()
            count_atoms(inputs[2*i], topdat, rdtp)
            if protein[i]:
                ndup.append(int(inputs[2*i+1]))
                topdat[rdtp].nmol = 1
                for j in range(1,ndup[i]):
                    topdat[rdtp+j].nmol = 1
                    topdat[rdtp+j].nat = topdat[rdtp].nat 
                    topdat[rdtp+j].ngo = topdat[rdtp].ngo 
                    topdat[rdtp+j].nbnd = topdat[rdtp].nbnd 
                    topdat[rdtp+j].nang = topdat[rdtp].nang 
                    topdat[rdtp+j].nimp = topdat[rdtp].nimp 
                    topdat[rdtp+j].ndih = topdat[rdtp].ndih 
                rdtp += ndup[i]
            else:
                ndup.append(1)
                topdat[rdtp].nmol = int(inputs[2*i+1])
                rdtp += 1
        rdtp = 0
        for idx in range(tmp_ntops):
            print("Topfile {}".format(inputs[2*idx]))
            print("Found: {} Atoms".format(topdat[rdtp].nat))
            print("Found: {} Bonds".format(topdat[rdtp].nbnd))
            print("Found: {} Angles".format(topdat[rdtp].nang))
            print("Found: {} Dihedrals".format(topdat[rdtp].ndih))
            print("Found: {} Impropers".format(topdat[rdtp].nimp))
            print()
            sysdat.nats  += topdat[rdtp].nat
            sysdat.ngo   += topdat[rdtp].ngo
            sysdat.nbnds += topdat[rdtp].nbnd
            sysdat.nangs += topdat[rdtp].nang
            sysdat.nimps += topdat[rdtp].nimp
            sysdat.ndihs += topdat[rdtp].ndih
            sysdat.total_ats  += topdat[rdtp].nat*topdat[rdtp].nmol*ndup[idx]
            sysdat.total_go   += topdat[rdtp].ngo*topdat[rdtp].nmol*ndup[idx]
            sysdat.total_bnds += topdat[rdtp].nbnd*topdat[rdtp].nmol*ndup[idx]
            sysdat.total_angs += topdat[rdtp].nang*topdat[rdtp].nmol*ndup[idx]
            sysdat.total_imps += topdat[rdtp].nimp*topdat[rdtp].nmol*ndup[idx]
            sysdat.total_dihs += topdat[rdtp].ndih*topdat[rdtp].nmol*ndup[idx]
            rdtp += ndup[idx]
        print("Totals:")
        print("Found: {} Atoms".format(sysdat.total_dihs))
        print("Found: {} Bonds".format(sysdat.total_dihs))
        print("Found: {} Angles".format(sysdat.total_dihs))
        print("Found: {} Dihderals".format(sysdat.total_dihs))
        print("Found: {} Impropers".format(sysdat.total_imps))
        rdtp = 0
        for idx in range(tmp_ntops):
            if protein[idx]:
                read_top_Go(inputs[2*idx], sysdat, topdat, rdtp, ndup[idx], bbind)
                rdtp += ndup[idx]
            else:
                read_top(inputs[2*idx], sysdat, topdat, rdtp)
                rdtp += ndup[idx]
        sysdat.ntops = rdtp
    else:
        print("Will read {} topology file(s).".format(ntops))
        print()
        rdtp = 0
        # Loop through the topologies and count the number of atoms, bonds and bends
        while rdtp < ntops:
            topdat[rdtp].nmol = int(inputs[(2*rdtp) + 1])
            count_atoms(inputs[2*rdtp], topdat, rdtp)
            rdtp += 1
        sysdat.ntops = ntops
        for idx in range(ntops):
            print("Topfile {}".format(inputs[2*idx]))
            print("Found: {} Atoms".format(topdat[idx].nat))
            print("Found: {} Bonds".format(topdat[idx].nbnd))
            print("Found: {} Angles".format(topdat[idx].nang))
            print("Found: {} Dihedrals".format(topdat[idx].ndih))
            print("Found: {} Impropers".format(topdat[idx].nimp))
            print()
            sysdat.nats  += topdat[idx].nat
            sysdat.nbnds += topdat[idx].nbnd
            sysdat.nangs += topdat[idx].nang
            sysdat.ndihs += topdat[idx].ndih
            sysdat.nimps += topdat[idx].nimp
            sysdat.total_ats  += topdat[idx].nat*topdat[idx].nmol
            sysdat.total_bnds += topdat[idx].nbnd*topdat[idx].nmol
            sysdat.total_angs += topdat[idx].nang*topdat[idx].nmol
            sysdat.total_dihs += topdat[idx].ndih*topdat[idx].nmol
            sysdat.total_imps += topdat[idx].nimp*topdat[idx].nmol
        print("Totals:")
        print("Found: {} Atoms".format(sysdat.total_ats))
        print("Found: {} Bonds".format(sysdat.total_bnds))
        print("Found: {} Angles".format(sysdat.total_angs))
        print("Found: {} Dihderals".format(sysdat.total_dihs))
        print("Found: {} Impropers".format(sysdat.total_imps))
        rdtp = 0
        while rdtp < ntops:
            topdat[rdtp].nmol = int(inputs[(2*rdtp + 1)])
            read_top(inputs[2*rdtp], sysdat, topdat, rdtp)
            rdtp += 1
    read_database(inputs[nargs-2], database)
    print("###########################")
    print("######  Database Summary ######")
    print("Found {} Unique Vdw pair params".format(database.nvdwtype))
    print("Found {} Unique Bond params".format(database.nbndtype))
    print("Found {} Unique Angle params".format(database.nangtype))
    print("Found {} Unique Dihdral params".format(database.ndihtype))
    print("Found {} Unique Improper params".format(database.nimptype))
    # read a input file formatted in PDB
    read_pdb(inputs[nargs-1], sysdat)
    # find parameters and write PARM.FILE included in LAMMPS inputs
    get_unique(database, topdat, sysdat) 
    # write DATA.FILE including topology, configuration, and box size info for LAMMPS simulations
    read_coords(topdat, sysdat)
    # write a file formatted in PSF for visualization
    write_psf(topdat, sysdat)


if __name__ == "__main__":
    args = get_option()
    run(args)
