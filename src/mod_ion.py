#!/usr/bin/env python

import sys, json, random, math
from pathlib import Path
import numpy as np
from argparse import ArgumentParser

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument('input', type=str,
                            help='Specify input CG PDB file name.')
    argparser.add_argument('output', type=str,
                            help='Specify output CG PDB file name.')
    argparser.add_argument('-conc', type=float,
                            default=0.15,
                            help='NaCl concentration [M] (default: 0.15).')
    argparser.add_argument('-soluq', type=int,
                            default=0,
                            help='total charge of solute molecules [e] (default: 0).')
    argparser.add_argument('-pspica', action='store_true',
                            help='run for pSPICA FF (default: for SPICA FF).')
    argparser.add_argument('-verbose', type=int,
                            default=1,
                            help='activate verbose logging, 0 : off, 1 : on (default: 1).')
    return argparser.parse_args()

def get_option_script(argv):
    argparser = ArgumentParser(usage='modion [-h] [-conc CONC] [-soluq SOLUQ] [-verbose VERBOSE] input output',
                               prog ="modion")
    argparser.add_argument('input', type=str,
                            help='Specify input CG PDB file name.')
    argparser.add_argument('output', type=str,
                            help='Specify output CG PDB file name.')
    argparser.add_argument('-conc', type=float,
                            default=0.15,
                            help='NaCl concentration [M] (defualt 0.15).')
    argparser.add_argument('-soluq', type=int,
                            default=0,
                            help='total charge of solute molecules [e] (default: 0).')
    argparser.add_argument('-pspica', action='store_true',
                            help='run for pSPICA FF (default: for SPICA FF).')
    argparser.add_argument('-verbose', type=int,
                            default=1,
                            help='activate verbose logging, 0 : off, 1 : on (default: 1).')
    return argparser.parse_args(argv)

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

def read_pdb(infile):
    pdb_list = []
    cryst    = []
    ters     = []
    natom    = 0
    try:
        f = open(infile,"r")
    except:
        sys.exit("ERROR: FILE",infile,"IS NOT FOUND")
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
            resid   = line[22:27]
            code    = " "
            #resid   = line[22:26]
            #code    = line[26:27]
            posX   = float(line[30:38])
            posY   = float(line[38:46])
            posZ   = float(line[46:54])
            occup  = line[54:60]
            Tfact  = line[60:66]
            segid  = line[72:76]
            elesym = line[76:78].split("\n")[0]
            #charge = line[78:80]
            pdb_list.append([recname,index,atmname,indicat,resname,chainid,resid,code,posX,posY,posZ,occup,Tfact,segid,elesym])
        elif recname == "TER":
            ters.append(resid)
        line = f.readline()
    f.close()
    return natom, pdb_list, cryst, ters

def rem(lst, xs):
    for x in xs:
        lst.remove(x)
    return lst

class mod_ion:
    def __init__(self, infile, outfile, conc, soluq, pspica, verbose):
        self.infile            = infile
        self.outfile           = outfile
        self.conc              = conc
        self.soluq             = soluq
        self.pspica            = pspica
        self.verbose           = verbose
        self.cryst             = []
        self.nwat              = 0
        self.nsod              = 0
        self.ncla              = 0
    
    def _collect_solvent(self):
        solvent_idx = []
        solvent_name = []
        solute_idx = []
        for idx, data in enumerate(self.pdb_data):
            aname = data[PDB_ATMNAME].strip()
            if (aname in ["W"] and not self.pspica) or (aname in ["W1"] and self.pspica):
                self.nwat += 1
                solvent_idx.append(idx)
                solvent_name.append("wat")
            elif (aname in ["SOD"] and not self.pspica) or (aname in ["SOD1"] and self.pspica):
                self.nsod += 1
                solvent_idx.append(idx)
                solvent_name.append("sod")
            elif (aname in ["CLA"] and not self.pspica) or (aname in ["CLA1"] and self.pspica):
                self.ncla += 1
                solvent_idx.append(idx)
                solvent_name.append("cla")
            elif aname in ["W2", "SOD2", "CLA2"] and self.pspica:
                pass
            else:
                solute_idx.append(idx)
        self.solvent_index = np.array(solvent_idx)
        self.solvent_name = np.array(solvent_name)
        self.solute_index = solute_idx
                
    def _write_init_info(self):
        self._collect_solvent()
        print("# of CG WAT:", self.nwat)
        print("# of CG SOD:", self.nsod)
        print("# of CG CLA:", self.ncla)

    def _convert_atomic(self):
        self.nwat_aa = 3.0*(self.nwat + self.nsod) + 2.0*self.ncla
        self.nsod_aa = self.nsod
        self.ncla_aa = self.ncla
    
    def _calc_nmol(self):
        conc = self.conc
        rho  = 0.997062 # bulk water density at 298.15 [K]
        mwat = 18.01528 # water molecular weight
        msod = 22.9898
        mcla = 35.453
        ntot0 = self.nwat + self.nsod + self.ncla
        if self.nsod > self.ncla:
            nnacl = self.ncla
            ncion = self.nsod - self.ncla
        else:
            nnacl = self.nsod
            ncion = -(self.ncla - self.nsod)
        conc_prev = nnacl*rho/(self.nwat_aa*mwat)*10.0**3
        self.conc_prev = conc_prev
        print(f"Current conc.: {round(conc_prev, 4)} [M]")
        print(f"Total solvent CG particles: {ntot0}")
        print()
        nall_aa = self.nwat_aa + self.nsod_aa + self.ncla_aa
        nact_aa = nall_aa - self.soluq
        cnt = 0
        ntot = 1e10
        while ntot > ntot0:
            nact_aa_0 = nact_aa - cnt
            nnacl_ret = 0
            for _ in range(10):
                nwat_aa_wk = nact_aa_0 - nnacl_ret
                nnacl_ret = round(nwat_aa_wk*mwat/(rho*10.0**3)*conc)
            conc_new = nnacl_ret*rho/(nwat_aa_wk*mwat)*10.0**3 
            if self.soluq > 0:
                nwat_cg = int((nwat_aa_wk - 5*nnacl_ret - 2*self.soluq)/3.0)
                ntot = nwat_cg + 2*nnacl_ret + self.soluq
            else:
                nwat_cg = int((nwat_aa_wk - 5*nnacl_ret + 3*self.soluq)/3.0)
                ntot = nwat_cg + 2*nnacl_ret - self.soluq
            cnt += 1
        self.conc_new = conc_new
        print("# of CG WAT:", nwat_cg)
        self.new_nwat = nwat_cg
        if self.soluq > 0:
            self.new_nsod = nnacl_ret
            self.new_ncla = nnacl_ret + self.soluq
        else:
            self.new_nsod = nnacl_ret - self.soluq
            self.new_ncla = nnacl_ret 
        print("# of CG SOD:", self.new_nsod)
        print("# of CG CLA:", self.new_ncla)
        print(f"New conc. {round(conc_new, 4)} [M]")
        print(f"Total solvent CG particles: {ntot}")
        print()


    def _mod_solvent(self):
        dsod = self.new_nsod - self.nsod
        dcla = self.new_ncla - self.ncla
        sol_index = {}
        sod2wat_index = []
        cla2wat_index = []
        for sol in ["wat", "sod", "cla"]:
            sol_index[sol] = list(self.solvent_index[self.solvent_name == sol])
        if dsod > 0:
            tmp_index = sorted(random.sample(sol_index["wat"], 
                               dsod), reverse=True)
            sol_index["sod"] += tmp_index
            rem(sol_index["wat"], tmp_index)
        else: 
            sod2wat_index = sorted(random.sample(sol_index["sod"], 
                               -dsod), reverse=True)
            rem(sol_index["sod"], sod2wat_index)
        if dcla > 0:
            tmp_index = sorted(random.sample(sol_index["wat"], 
                               dcla), reverse=True)
            sol_index["cla"] += tmp_index
            rem(sol_index["wat"], tmp_index)
        else: 
            cla2wat_index = sorted(random.sample(sol_index["cla"], 
                               -dcla), reverse=True)
            rem(sol_index["cla"], cla2wat_index)
        sol_index["wat"] += sod2wat_index
        sol_index["wat"] += cla2wat_index
        self.sol_index = sol_index
    
    def _write_pdb(self):
        nl = len(self.pdb_data)
        cnt = 1 
        recname = "ATOM"
        indicat = " " 
        code = " " 
        chainid = " " 
        ter_id = None
        ter_flag = False
        f = open(self.outfile, "w")
        print(f"REMARK Generated by mod_ions.py, {round(self.conc_prev,4)} -> {round(self.conc_new,4)} [M]", file=f)
        print(f"REMARK #W {self.nwat} -> {self.new_nwat}",file=f)
        print(f"REMARK #S {self.nsod} -> {self.new_nsod}",file=f)
        print(f"REMARK #C {self.ncla} -> {self.new_ncla}",file=f)
        if len(self.cryst) == 1:
            print('%s' % self.cryst[0],file=f)
        for i in self.solute_index:
            atmname  = self.pdb_data[i][PDB_ATMNAME]
            resname  = self.pdb_data[i][PDB_RESNAME]
            chainid  = self.pdb_data[i][PDB_CHAINID]
            resid    = self.pdb_data[i][PDB_RESID]
            posX     = self.pdb_data[i][PDB_POSX]
            posY     = self.pdb_data[i][PDB_POSY]
            posZ     = self.pdb_data[i][PDB_POSZ]
            occup    = self.pdb_data[i][PDB_OCCUP]
            if len(occup.split()) == 0:
                occup = 0.0
            else:
                occup = float(occup)
            tfact    = self.pdb_data[i][PDB_TFACT]
            if len(tfact.split()) == 0:
                tfact = 0.0
            else:
                tfact = float(occup)
            segid    = self.pdb_data[i][PDB_SEGID].strip()
            index    = str(cnt)
            if ter_flag and ter_id != resid:
               print("TER", file=f)
               ter_flag = False
            print('%-6s%5s %4s%1s%-4s%1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%10s' \
                        %(recname,index[0:5],atmname,indicat,resname,chainid,resid,posX,posY,posZ,occup,tfact,segid), file=f)
            if resid in self.ters and not ter_flag:
                ter_id = resid
                ter_flag = True
            cnt += 1
        if self.pspica:
            np = 2
        else:
            np = 1
        for k in sorted(self.sol_index["wat"]):
            for j in range(np):  
                i = k + j
                if self.pspica:
                    atmname  = "W" + str(j+1)
                    resname  = "PWAT"
                else:
                    atmname  = "W"
                    resname  = "WAT"
                chainid  = self.pdb_data[i][PDB_CHAINID]
                resid    = self.pdb_data[i][PDB_RESID]
                posX     = self.pdb_data[i][PDB_POSX]
                posY     = self.pdb_data[i][PDB_POSY]
                posZ     = self.pdb_data[i][PDB_POSZ]
                occup    = self.pdb_data[i][PDB_OCCUP]
                if len(occup.split()) == 0:
                    occup = 0.0
                else:
                    occup = float(occup)
                tfact    = self.pdb_data[i][PDB_TFACT]
                if len(tfact.split()) == 0:
                    tfact = 0.0
                else:
                    tfact = float(occup)
                segid    = self.pdb_data[i][PDB_SEGID].strip()
                index    = str(cnt)
                print('%-6s%5s %4s%1s%-4s%1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%10s' \
                            %(recname,index[0:5],atmname,indicat,resname,chainid,resid,posX,posY,posZ,occup,tfact,segid), file=f)
                cnt += 1
        for k in sorted(self.sol_index["sod"]):
            for j in range(np):  
                i = k + j
                if self.pspica:
                    atmname  = "SOD" + str(j+1)
                    resname  = "PSOD"
                else:
                    atmname  = "SOD"
                    resname  = "SOD"
                chainid  = self.pdb_data[i][PDB_CHAINID]
                resid    = self.pdb_data[i][PDB_RESID]
                posX     = self.pdb_data[i][PDB_POSX]
                posY     = self.pdb_data[i][PDB_POSY]
                posZ     = self.pdb_data[i][PDB_POSZ]
                occup    = self.pdb_data[i][PDB_OCCUP]
                if len(occup.split()) == 0:
                    occup = 0.0
                else:
                    occup = float(occup)
                tfact    = self.pdb_data[i][PDB_TFACT]
                if len(tfact.split()) == 0:
                    tfact = 0.0
                else:
                    tfact = float(occup)
                segid    = self.pdb_data[i][PDB_SEGID].strip()
                index    = str(cnt)
                print('%-6s%5s %4s%1s%-4s%1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%10s' \
                            %(recname,index[0:5],atmname,indicat,resname,chainid,resid,posX,posY,posZ,occup,tfact,segid), file=f)
                cnt += 1
        for k in sorted(self.sol_index["cla"]):
            for j in range(np):  
                i = k + j
                if self.pspica:
                    atmname  = "CLA" + str(j+1)
                    resname  = "PCLA"
                else:
                    atmname  = "CLA"
                    resname  = "CLA"
                chainid  = self.pdb_data[i][PDB_CHAINID]
                resid    = self.pdb_data[i][PDB_RESID]
                posX     = self.pdb_data[i][PDB_POSX]
                posY     = self.pdb_data[i][PDB_POSY]
                posZ     = self.pdb_data[i][PDB_POSZ]
                occup    = self.pdb_data[i][PDB_OCCUP]
                if len(occup.split()) == 0:
                    occup = 0.0
                else:
                    occup = float(occup)
                tfact    = self.pdb_data[i][PDB_TFACT]
                if len(tfact.split()) == 0:
                    tfact = 0.0
                else:
                    tfact = float(occup)
                segid    = self.pdb_data[i][PDB_SEGID].strip()
                index    = str(cnt)
                print('%-6s%5s %4s%1s%-4s%1s%5s   %8.3f%8.3f%8.3f%6.2f%6.2f%10s' \
                            %(recname,index[0:5],atmname,indicat,resname,chainid,resid,posX,posY,posZ,occup,tfact,segid), file=f)
                cnt += 1
        print('%s' % "END", file=f)
        f.close()

    def run(self):
        self.natom, self.pdb_data, self.cryst, self.ters = read_pdb(self.infile)
        self._write_init_info()
        self._convert_atomic()
        self._calc_nmol()
        self._mod_solvent()
        self._write_pdb()

if __name__ == "__main__":
    args     = get_option()
    infile   = args.input
    outfile  = args.output
    conc     = args.conc
    soluq    = args.soluq
    pspica   = args.pspica
    verbose  = args.verbose

    modion = mod_ion(infile, outfile, conc, soluq, pspica, verbose)
    modion.run()
