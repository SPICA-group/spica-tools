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
        
def read_coords(fname, database, topdat, sysdat):
    global ischarged
    if ischarged:
        print("FOUND CHARGES.")
    else:
        print("DID NOT FIND CHARGES.")
    if sysdat.total_ats != sysdat.foundatoms:
       print("ERROR: # OF ATOMS READ FROM TOPOLOGY AND COMMAND LINE.")
       sys.exit("DOES NOT MATCH # OF ATOMS FOUND IN THE COORD FILE.")

    with open("DATA.FILE", "w") as fout:
        print("LAMMPS description", file=fout)
        print(file=fout)
        print("{:<7} atoms".format(sysdat.total_ats), file=fout)
        if sysdat.nbnds > 0:
            print("{:<7} bonds".format(sysdat.total_bnds), file=fout)
        if sysdat.nangs  > 0:
            print("{:<7} angles".format(sysdat.total_angs), file=fout)
        if sysdat.total_diheds > 0:
            print("{:<7} dihedrals".format(sysdat.total_diheds), file=fout)
        if sysdat.total_improps > 0:
            print("{:<7} impropers".format(sysdat.total_improps), file=fout)
        print(file=fout)
        print("{:<7} atom types".format(sysdat.uniq_nats), file=fout)
        if sysdat.nbnds > 0:
                print("{:<7} bond types".format(sysdat.uniq_nbnds), file=fout)
        if sysdat.nangs > 0:
                print("{:<7} angle types".format(sysdat.uniq_nangs), file=fout)
        if sysdat.total_diheds > 0:
                print("{:<7} dihedral types".format(sysdat.uniq_ndiheds), file=fout)
        if sysdat.total_improps > 0:
                print("{:<7} improper types".format(sysdat.uniq_nimprops), file=fout)
        print(file=fout)
        print("{} {} xlo xhi".format(-0.5*sysdat.boxx, 0.5*sysdat.boxx), file=fout)
        print("{} {} ylo yhi".format(-0.5*sysdat.boxy, 0.5*sysdat.boxy), file=fout)
        print("{} {} zlo zhi".format(-0.5*sysdat.boxz, 0.5*sysdat.boxz), file=fout)
        print(file=fout)
        print("Atoms", file=fout)
        print(file=fout)
        atidx = molidx = 0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nmol):
                molidx += 1
                for kdx in range(topdat[idx].nat):
                    atidx += 1
                    print("{} {} {} {:7.4f} {:7.4f} {:7.4f} {:7.4f} # {}".format(atidx, molidx, topdat[idx].parm_atomtype[kdx]+1,
                           topdat[idx].charge[kdx], sysdat.coordx[atidx-1], sysdat.coordy[atidx-1], sysdat.coordz[atidx-1],
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
                        print("{} {} {} {} # {} {}".format(bondidx,topdat[idx].bndtype[kdx]+1,
                               topdat[idx].bndndx1[kdx]+(jdx*topdat[idx].nat)+offset,topdat[idx].bndndx2[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].atomtype[topdat[idx].bndndx1[kdx]-1],topdat[idx].atomtype[topdat[idx].bndndx2[kdx]-1]), file=fout)
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
                        print("{} {} {} {} {} # {} {} {}".format(angleidx,topdat[idx].angtype[kdx]+1,
                               topdat[idx].angndx1[kdx]+(jdx*topdat[idx].nat)+offset,topdat[idx].angndx2[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].angndx3[kdx]+(jdx*topdat[idx].nat)+offset,topdat[idx].atomtype[topdat[idx].angndx1[kdx]-1],
                               topdat[idx].atomtype[topdat[idx].angndx2[kdx]-1],topdat[idx].atomtype[topdat[idx].angndx3[kdx]-1]), file=fout)
                offset += topdat[idx].nmol*topdat[idx].nat
        if sysdat.total_diheds > 0:
            print(file=fout)
            print("Dihedrals", file=fout)
            print(file=fout)
            dihedidx = 0
            offset   = 0
            for idx in range(sysdat.ntops):
                for jdx in range(topdat[idx].nmol):
                    molidx += 1
                    for kdx in range(topdat[idx].ndihed):
                        dihedidx += 1
                        print("{} {} {} {} {} {} # {} {} {} {}".format(dihedidx,topdat[idx].dihedtype[kdx]+1,
                               topdat[idx].dihedndx1[kdx]+(jdx*topdat[idx].nat)+offset,topdat[idx].dihedndx2[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].dihedndx3[kdx]+(jdx*topdat[idx].nat)+offset,topdat[idx].dihedndx4[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].atomtype[topdat[idx].dihedndx1[kdx]-1],topdat[idx].atomtype[topdat[idx].dihedndx2[kdx]-1], 
                               topdat[idx].atomtype[topdat[idx].dihedndx3[kdx]-1],topdat[idx].atomtype[topdat[idx].dihedndx4[kdx]-1]), file=fout)
                offset += topdat[idx].nmol*topdat[idx].nat
        if sysdat.total_improps > 0:
            print(file=fout)
            print("Impropers", file=fout)
            print(file=fout)
            impropidx = 0
            offset    = 0
            for idx in range(sysdat.ntops):
                for jdx in range(topdat[idx].nmol):
                    molidx += 1
                    for kdx in range(topdat[idx].nimprop):
                        impropidx += 1
                        print("{} {} {} {} {} {} # {} {} {} {}".format(impropidx,topdat[idx].improptype[kdx]+1,
                               topdat[idx].impropndx1[kdx]+(jdx*topdat[idx].nat)+offset,topdat[idx].impropndx2[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].impropndx3[kdx]+(jdx*topdat[idx].nat)+offset,topdat[idx].impropndx4[kdx]+(jdx*topdat[idx].nat)+offset,
                               topdat[idx].atomtype[topdat[idx].impropndx1[kdx]-1],topdat[idx].atomtype[topdat[idx].impropndx2[kdx]-1],
                               topdat[idx].atomtype[topdat[idx].impropndx3[kdx]-1],topdat[idx].atomtype[topdat[idx].impropndx4[kdx]-1]), file=fout)
                offset += topdat[idx].nmol*topdat[idx].nat;

def write_psf(fname, database, topdat, sysdat):
    global ischarged
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
        print("{:<7} !NPHI: dihedrals".format(sysdat.total_diheds),    file=fout)
        print("{:<7} !NIMPHI: impropers".format(sysdat.total_improps), file=fout)
        print("       0 !NDON: donors", file=fout)
        print(file=fout)
        print("       0 !NACC: acceptors", file=fout)
        print(file=fout)
        
def cmp_wc(s1, s2):
    if s1[-1] == "*" or s2[-1] == "*":
        return s1[:-1] == s2[:-1]       
    else:
        return s1 == s2
        
def get_unique(database, topdat, sysdat):
    global ischarged
    uniq_nats = uniq_bnds = uniq_angs = uniq_diheds = uniq_improps = 0
    uniq_atype, uniq_mass, uniq_charge, bnd_params, ang_params, ang_vdw = [], [], [], [], [], []
    fout = open("PARM.FILE", "w")
    print("# Generated by setup_lammps", file=fout)
    print(file=fout)
    if ischarged:
        print("pair_style      lj/sdk/coul/long        15.0", file=fout)
    else:
        print("pair_style      lj/sdk        15.0", file=fout)
    if sysdat.nbnds > 0:
        print("bond_style      harmonic", file=fout)
    if sysdat.nangs > 0:
        print("angle_style     hybrid sdk harmonic", file=fout)
    if sysdat.total_diheds > 0:
        print("dihedral_style  charmm", file=fout)
    if sysdat.total_improps > 0:
        print("improper_style  harmonic", file=fout)
    print("special_bonds   lj/coul 0.0 0.0 1.0", file=fout)
    print(file=fout)
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
    print("FOUND {} UNIQUE ATOMS.".format(uniq_nats))
    sysdat.uniq_nats = uniq_nats
    for idx in range(uniq_nats):
        print("mass   {:<5}  {:6.4f} # {}".format(idx+1, uniq_mass[idx],uniq_atype[idx]), file=fout)
    print(file=fout)
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
                    print("*********************");
                    print("WARNING: No params for VDW interaction between {} and {}".format(uniq_atype[idx],uniq_atype[jdx]))
                    print("WILL CHECK FOR WILDCARDS.");
            elif ifound == 1:
                    print("pair_coeff  {:<5} {:<5} {:<6} {:5.4f} {:5.4f} # {:<4} {:<4}".format(idx+1,jdx+1,database.vdwstyle[vdwtmp],
                           database.eps[vdwtmp],database.sig[vdwtmp],database.vdwtype1[vdwtmp],database.vdwtype2[vdwtmp]),
                           file=fout)
    print(file=fout)
    # BONDS
    if sysdat.nbnds > 0:
        index0 = 0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nbnd):
                datndx = -1
                # AT THIS POINT WE WILL CHECK TO SEE IF THE PARAMS WERE GIVEN 
                # IN THE TOP FILE....IF SO WE WILL SKIP A LOT OF THIS AND ADD 
                # THIS AS A UNIQUE BOND....IF NOT WE GO THROUGH THE PROCEDURE
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
                            sys.exit("ERROR: DID NOT FIND BOND PARAMETERS IN DATABASE {} {} {} {}".format(
                                    topdat[idx].bndndx1[jdx],
                                    topdat[idx].bndndx2[jdx],
                                    topdat[idx].atomtype[topdat[idx].bndndx1[jdx]-1],
                                    topdat[idx].atomtype[topdat[idx].bndndx2[jdx]-1]))
                    # Now make sure we do not already know we have this interaction 
                    ikeep = 1
                    if database.bnde[datndx] > 0:
                        for kdx in range(uniq_bnds):
                            if datndx == bnd_params[kdx]:
                                ikeep = 0 # found a replica
                                topdat[idx].bndtype.append(kdx)
                                break # kill the for loop
                    # ikeep = 1 if we found a new one
                    if ikeep == 1:
                        bnd_params.append(datndx)
                        sysdat.param_bnds.append(datndx)
                        topdat[idx].bndtype.append(uniq_bnds)
                        uniq_bnds += 1
                        #print("#####  {}".format(datndx))
                        if database.bnde[datndx] > 0:
                            print("bond_coeff  {:<6} {:8.4f} {:8.4f} # {} {}".format(uniq_bnds,database.fbnd[bnd_params[uniq_bnds-1]],
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
                            print("bond_coeff  {:<6} {:8.4f} {:8.4f} # {} {}".format(uniq_bnds,database.fbnd[bnd_params[uniq_bnds-1]],
                                   bond_in_pdb,
                                   database.bndtype1[bnd_params[uniq_bnds-1]],
                                   database.bndtype2[bnd_params[uniq_bnds-1]]), file=fout)
                else:
                    # THE PARAMS WERE GIVEN IN THE TOP FILE SO LETS ADD IT TO THE PARAM FILE */
                    topdat[idx].bndtype.append(uniq_bnds)
                    bnd_params.append(-1)
                    uniq_bnds += 1
                    print("bond_coeff  {:<6} {:8.4f} {:8.4f} # {} {} FROM TOP".format(uniq_bnds,
                           topdat[idx].bndfk[jdx],topdat[idx].bndeq[jdx],
                           topdat[idx].atomtype[topdat[idx].bndndx1[jdx]-1], 
                           topdat[idx].atomtype[topdat[idx].bndndx2[jdx]-1]), file=fout)
            index0 += topdat[idx].nat*topdat[idx].nmol
        # FINISHED LOOPING OVER TOP FILES */
    sysdat.uniq_nbnds = uniq_bnds
    print(file=fout)
    # ANGLES
    if sysdat.nangs > 0:
        uniq_angs = index0 = 0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nang):
                datndx = -1
                # AT THIS POINT WE WILL CHECK TO SEE IF THE PARAMS WERE GIVEN */
                # IN THE TOP FILE....IF SO WE WILL SKIP A LOT OF THIS AND ADD */
                # THIS AS A UNIQUE BOND....IF NOT WE GO THROUGH THE PROCEDURE */
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
                        print("UPDATE DATABASE.")
                    # end VDW for CG angles 
                    # No params for this interaction in the database
                    if datndx == -1:
                        print("ERROR: DID NOT FIND ANGLE PARAMETERS IN DATABASE {} {} {} ({} {} {})".format(topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1],
                               topdat[idx].angndx1[jdx],
                               topdat[idx].angndx2[jdx],
                               topdat[idx].angndx3[jdx]))
                        sys.exit(1)
                    # Now make sure we do not already have this one
                    ikeep = 1
                    if database.ange[datndx] > 0:
                        for kdx in range(uniq_angs):
                            if datndx == ang_params[kdx]:
                                ikeep = 0
                                topdat[idx].angtype.append(kdx)
                                break
                    if ikeep == 1:
                        ang_params.append(datndx)
                        ang_vdw.append(vdwtmp)
                        sysdat.param_angs.append(datndx)
                        topdat[idx].angtype.append(uniq_angs)
                        uniq_angs += 1
                        if database.ange[datndx] > 0:
                            if database.angsdk[ang_params[uniq_angs-1]]:
                                print("angle_coeff {:<10}      sdk  {:8.4f} {:8.4f} {} {:8.4f} {:8.4f} # {} {} {}".format(uniq_angs,
                                       database.fang[ang_params[uniq_angs-1]],
                                       database.ange[ang_params[uniq_angs-1]],
                                       database.vdwstyle[ang_vdw[uniq_angs-1]],
                                       database.eps[ang_vdw[uniq_angs-1]],
                                       database.sig[ang_vdw[uniq_angs-1]],
                                       database.angtype1[ang_params[uniq_angs-1]],
                                       database.angtype2[ang_params[uniq_angs-1]],
                                       database.angtype3[ang_params[uniq_angs-1]]), file=fout)
                            else:
                                print("angle_coeff {:<10} harmonic  {:8.4f} {:8.4f} # {} {} {}".format(uniq_angs,
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
                                print("angle_coeff {:<10}      sdk  {:8.4f} {:8.4f} {} {:8.4f} {:8.4f} # {} {} {}".format(uniq_angs,
                                       database.fang[ang_params[uniq_angs-1]],
                                       angle_in_pdb,
                                       database.vdwstyle[ang_vdw[ uniq_angs-1]],
                                       database.eps[ang_vdw[uniq_angs-1]],
                                       database.sig[ang_vdw[uniq_angs-1]],
                                       database.angtype1[ang_params[uniq_angs-1]],
                                       database.angtype2[ang_params[uniq_angs-1]],
                                       database.angtype3[ang_params[uniq_angs-1]]), file=fout)
                            else:
                                print("angle_coeff {:<10} harmonic  {:8.4f} {:8.4f} # {} {} {}".format(uniq_angs,
                                       database.fang[ang_params[uniq_angs-1]],
                                       angle_in_pdb,
                                       database.angtype1[ang_params[uniq_angs-1]],
                                       database.angtype2[ang_params[uniq_angs-1]],
                                       database.angtype3[ang_params[uniq_angs-1]]), file=fout)
                else:
                    # THIS PARAM WAS SET IN THE TOP FILE
                    # still need vdw stuff
                    # RHD Get the VDW for the CG angles
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
                        print("*********************")
                        print("ERROR: No params for VDW interaction between {} and {} for angle (topfile)".format(topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1],
                               topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1]))
                        print("UPDATE DATABASE.")
                        sys.exit(1)
                    # end VDW for CG angles
                    topdat[idx].angtype[jdx] = uniq_angs
                    uniq_angs += 1
                    if database.angsdk[ang_params[uniq_angs-1]]:
                        print("angle_coeff {:<10}      sdk  {:8.4f} {:8.4f} {} {:8.4f} {:8.4f} # {} {} {} FROM TOP".format(uniq_angs,
                        topdat[idx].angfk[jdx],
                        topdat[idx].angeq[jdx],
                        database.vdwstyle[vdwtmp],
                        database.eps[vdwtmp],
                        database.sig[vdwtmp],
                        topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1],
                        topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1],
                        topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1]), file=fout)
                    else:
                        print("angle_coeff {:<10} harmonic  {:8.4f} {:8.4f} # {} {} {} FROM TOP".format(uniq_angs,
                        topdat[idx].angfk[jdx],
                        topdat[idx].angeq[jdx],
                        topdat[idx].atomtype[topdat[idx].angndx1[jdx]-1],
                        topdat[idx].atomtype[topdat[idx].angndx2[jdx]-1],
                        topdat[idx].atomtype[topdat[idx].angndx3[jdx]-1]), file=fout)

            index0 += topdat[idx].nat*topdat[idx].nmol;
    sysdat.uniq_nangs = uniq_angs
    print(file=fout)
    # NOW LET HANDLE THE DIHEDRAL PARAMS THE ONLY WAY WE DO...THEY HAVE TO BE SPECIFIED IN THE TOP FILE */
    print("DIHEDS TEST {}".format(sysdat.total_diheds))
    if sysdat.total_diheds > 0:
        uniq_diheds = index0 = 0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].ndihed):
                if topdat[idx].dihedpset[jdx] == 1:
                    topdat[idx].dihedtype.append(uniq_diheds) 
                    uniq_diheds += 1
                    i1 = topdat[idx].dihedndx1[jdx]-1 + index0;
                    i2 = topdat[idx].dihedndx2[jdx]-1 + index0;
                    i3 = topdat[idx].dihedndx3[jdx]-1 + index0;
                    i4 = topdat[idx].dihedndx4[jdx]-1 + index0;
                    r1 = np.array([sysdat.coordx[i1], sysdat.coordy[i1], sysdat.coordz[i1]])
                    r2 = np.array([sysdat.coordx[i2], sysdat.coordy[i2], sysdat.coordz[i2]])
                    r3 = np.array([sysdat.coordx[i3], sysdat.coordy[i3], sysdat.coordz[i3]])
                    r4 = np.array([sysdat.coordx[i4], sysdat.coordy[i4], sysdat.coordz[i4]])
                    dihedral_in_pdb = 180.0 + 180.0/pi*get_dihedral(r1,r2,r3,r4);
                    if dihedral_in_pdb > 0.0:
                        idihedral_in_pdb = int(dihedral_in_pdb + 0.5)
                    else:
                        idihedral_in_pdb = int(dihedral_in_pdb - 0.5)
                    print("dihedral_coeff {:<10} {:8.4f} {:<3} {:<3} {:2.1f} # {} {} {} {}".format(uniq_diheds,
                           topdat[idx].dihedfk[jdx],
                           topdat[idx].dihedn[jdx],
                           topdat[idx].dihedeq[jdx],
                           topdat[idx].dihedof[jdx],
                           topdat[idx].atomtype[topdat[idx].dihedndx1[jdx]-1],
                           topdat[idx].atomtype[topdat[idx].dihedndx2[jdx]-1],
                           topdat[idx].atomtype[topdat[idx].dihedndx3[jdx]-1],
                           topdat[idx].atomtype[topdat[idx].dihedndx4[jdx]-1]), file=fout)

            index0 += topdat[idx].nat*topdat[idx].nmol;
    sysdat.uniq_ndiheds=uniq_diheds;
    print(file=fout)
    # NOW LET HANDLE THE IMPROPS PARAMS THE ONLY WAY WE DO...THEY HAVE TO BE SPECIFIED IN THE TOP FILE */
    print("IMPROPS TEST {}".format(sysdat.total_improps))
    if sysdat.total_improps > 0:
        uniq_improps = 0
        for idx in range(sysdat.ntops):
            for jdx in range(topdat[idx].nimprop):
                if topdat[idx].improppset[jdx] == 1:
                    topdat[idx].improptype[jdx] = uniq_improps
                    uniq_improps += 1
                    print("improper_coeff {:<10} {:8.4f} {:8.4f} # {} {} {} {} FROM TOP".format(uniq_improps,
                           topdat[idx].impropfk[jdx],
                           topdat[idx].impropeq[jdx],
                           topdat[idx].atomtype[topdat[idx].impropndx1[jdx]-1],
                           topdat[idx].atomtype[topdat[idx].impropndx2[jdx]-1],
                           topdat[idx].atomtype[topdat[idx].impropndx3[jdx]-1],
                           topdat[idx].atomtype[topdat[idx].impropndx4[jdx]-1]), file=fout)
    sysdat.uniq_nimprops = uniq_improps
    fout.close()

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
                        print("WARNING: FOUND DUP VDW PARAM {} {}".format(vdwtype1,vdwtype2))
                        ikeep = 0
                    elif vdwtype1 == database.vdwtype2[idx] and vdwtype2 == database.vdwtype1[idx]:
                        print("WARNING: FOUND DUP VDW PARAM {} {}".format(vdwtype1,vdwtype2))
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
                            print("WARNING: FOUND DUP BOND PARAM {} {}".format(bndtype1,bndtype2))
                            ikeep = 0
                        elif bndtype1 == database.bndtype2[idx] and bndtype2 == database.bndtype1[idx]:
                            print("WARNING: FOUND DUP BOND PARAM {} {}".format(bndtype1,bndtype2))
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
                                print("WARNING: FOUND DUP ANGLE PARAM {} {} {}".format(angtype1,angtype2,angtype3))
                                ikeep = 0
                        elif angtype3 == database.angtype1[idx] and angtype1 == database.angtype3[idx]:
                                print("WARNING: FOUND DUP ANGLE PARAM {} {} {}".format(angtype1,angtype2,angtype3))
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
def read_top(fname, topdat, ntop):
    global ischarged
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
                    sys.exit("ERROR at FILE {}, line {}".format(fname, lc))
                if topdat[ntop].charge[ndx]*topdat[ntop].charge[ndx] > 1e-5:
                    print("CHARGE IN TOP FILE {} {}".format(fname, topdat[ntop].charge[ndx]))
                    ischarged = 1
                ndx += 1
            if items[0] == "bond":
                try:
                    topdat[ntop].bndndx1.append(int(items[1]))
                    topdat[ntop].bndndx2.append(int(items[2]))
                    topdat[ntop].bndfk.append(None)
                    topdat[ntop].bndeq.append(None)
                except:
                    sys.exit("ERROR at FILE {}, line {}".format(fname, lc))
                topdat[ntop].bndpset.append(False)
                bndx += 1
            if items[0] == "bondparam":
                print("WARNING: Using bond parameters from the top file.")
                if len(items) < 5:
                    sys.exit("ERROR: Not enough args for bondparam: must be: ndx1 ndx2 fk eq.")
                try:
                    topdat[ntop].bndndx1.append(int(items[1]))
                    topdat[ntop].bndndx2.append(int(items[2]))
                    topdat[ntop].bndfk.append(float(items[3]))
                    topdat[ntop].bndeq.append(float(items[4]))
                except:
                    sys.exit("ERROR at FILE {}, line {}".format(fname, lc))
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
                    sys.exit("ERROR at FILE {}, line {}".format(fname, lc))
                topdat[ntop].angpset.append(-1)
                andx += 1
            if items[0] == "angleparam":
                print("WARNING: Using angle parameters from the top file.")
                if len(items) < 6:
                    sys.exit("ERROR: Not enough args for angleparam: must be: ndx1 ndx2 ndx3 fk eq.")
                try:
                    topdat[ntop].angndx1.append(int(items[1]))
                    topdat[ntop].angndx2.append(int(items[2]))
                    topdat[ntop].angndx3.append(int(items[3]))
                    topdat[ntop].angfk.append(float(items[4]))
                    topdat[ntop].angeq.append(float(items[5]))
                except:
                    sys.exit("ERROR at FILE {}, line {}".format(fname, lc))
                topdat[ntop].angpset.append(1)
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
                    sys.exit("ERROR at FILE {}, line {}".format(fname, lc))
                topdat[ntop].improppset.append(-1)
                indx += 1
            if items[0] == "improperparam":
                print("WARNING: Using improper parameters from the top file.")
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
                    sys.exit("ERROR at FILE {}, line {}".format(fname, lc))
                topdat[ntop].dihedpset.append(1)
                indx += 1
            if items[0] == "dihedralparam":
                print("WARNING: Using dihedral parameters from the top file.")
                if len(items) < 9:
                    sys.exit("ERROR: Not enough args for angleparam: must be: ndx1 ndx2 ndx3 fk n eq onefour.")
                try:
                    topdat[ntop].dihedndx1.append(int(items[1]))
                    topdat[ntop].dihedndx2.append(int(items[2]))
                    topdat[ntop].dihedndx3.append(int(items[3]))
                    topdat[ntop].dihedndx4.append(int(items[4]))
                    topdat[ntop].dihedfk.append(float(items[5]))
                    topdat[ntop].dihedn.append(int(items[6]))
                    topdat[ntop].dihedeq.append(float(items[7]))
                    topdat[ntop].dihedof.append(float(items[8]))
                except:
                    sys.exit("ERROR at FILE {}, line {}".format(fname, lc))
                topdat[ntop].dihedpset.append(1)
                dndx += 1
            line = fin.readline()

# Main routine. Call and allocate                                       
# The idea is to read in the topologies and then check the database for 
# all of the required interaction params.                               
def run():
    nargs = len(sys.argv)
    if nargs < 5:
       print("usage: setup_lmp <topfile 1> <nmol 1> [ <topfile 2> <nmol 2> ..... <topfile n> <nmol n>] <paramfile> <coordfile>");
       print("Prints out input files for a lammps run. Takes a pdb file as the coordfile");
       sys.exit(1)
    topdat    = [Topdat() for _ in range(1000)]
    database  = Database() 
    sysdat    = Sysdat()
    ischarged = 0
    ntops = int((nargs - 3)/2)
    print("WILL READ {} TOPOLOGY FILE(S).".format(ntops))
    print()
    rdtp = 0
    # Loop through the topologies and count the number of atoms, bonds and bends
    while rdtp < ntops:
        topdat[rdtp].nmol = int(sys.argv[(2*rdtp) + 2])
        count_atoms(sys.argv[(2*rdtp)+1], topdat, rdtp)
        rdtp += 1
    sysdat.ntops = ntops
    for idx in range(ntops):
        print("TOPFILE {}".format(sys.argv[(2*idx)+1]))
        print("FOUND: {} atoms".format(topdat[idx].nat))
        print("FOUND: {} bonds".format(topdat[idx].nbnd))
        print("FOUND: {} angles".format(topdat[idx].nang))                                                                             
        print("FOUND: {} impropers".format(topdat[idx].nimprop))
        print("FOUND: {} dihedrals".format(topdat[idx].ndihed))
        print()
        sysdat.nats     += topdat[idx].nat
        sysdat.nbnds    += topdat[idx].nbnd
        sysdat.nangs    += topdat[idx].nang
        sysdat.nimprops += topdat[idx].nimprop
        sysdat.ndiheds  += topdat[idx].ndihed
        
        sysdat.total_ats     += topdat[idx].nat*topdat[idx].nmol
        sysdat.total_bnds    += topdat[idx].nbnd*topdat[idx].nmol
        sysdat.total_angs    += topdat[idx].nang*topdat[idx].nmol
        sysdat.total_improps += topdat[idx].nimprop*topdat[idx].nmol
        sysdat.total_diheds  += topdat[idx].ndihed*topdat[idx].nmol
        
    print("TOTALS:")
    print("FOUND: {} impropers".format(sysdat.total_improps))
    print("FOUND: {} dihderals".format(sysdat.total_diheds))
    rdtp = 0
    while rdtp < ntops:
            topdat[rdtp].nmol = int(sys.argv[(2*rdtp + 2)])
            read_top(sys.argv[(2*rdtp) + 1], topdat, rdtp)
            rdtp += 1
    count_params(sys.argv[nargs-2], database)
    read_database(sys.argv[nargs-2], database)
    print("###########################")
    print("####  DATABASE SUMMARY ####")
    print("FOUND {} UNIQUE VDW PAIR PARAMS".format(database.nvdwtype))
    print("FOUND {} UNIQUE BOND PARAMS".format(database.nbndtype))
    print("FOUND {} UNIQUE ANGLE PARAMS".format(database.nangtype))
    read_pdb(sys.argv[nargs-1], sysdat)
    # write PARM.FILE
    get_unique(database, topdat, sysdat) 
    # write DATA.FILE
    read_coords(sys.argv[nargs-1], database, topdat, sysdat)
    write_psf(sys.argv[nargs-1], database, topdat, sysdat)

if __name__ == "__main__":
    run()
