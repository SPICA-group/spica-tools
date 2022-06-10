import MDAnalysis as mda
import numpy as np
import map_to_cg as m2c
from tqdm import tqdm
from pathlib import Path
from argparse import ArgumentParser
import warnings
warnings.simplefilter("ignore")

def get_option():
    json   = Path(__file__).parents[0] / "spica_top.json"
    outpdb = "fin.cg.pdb"
    argparser = ArgumentParser()
    argparser.add_argument('inputAAPDB', type=str,
                            help='Specify input AA PDB file name.')
    argparser.add_argument('inputAATRAJ', type=str,
                            help='Specify input AA TRAJ file name.')
    argparser.add_argument('outputCGTRAJ', type=str,
                            help='Specify output CG TRAJ file name.')
    argparser.add_argument('-outpdb', type=str,
                            default=outpdb,
                            help='output PDB file name (default: fin.cg.pdb).')
    argparser.add_argument('-json', type=str,
                            default=json,
                            help='input json file name (default: spica_top.json).')
    argparser.add_argument('-begin', type=int,
                            default=0,
                            help='begining frame of mappin (default: 0).')
    argparser.add_argument('-last', type=int,
                            default=-1,
                            help='last frame of mapping (default: -1).')
    argparser.add_argument('-nodelwat', action='store_true',
                            help='not delete excess water due to CG ion mapping (default: off).')
    return argparser.parse_args()

def get_option_script(argv):
    json   = Path(__file__).parents[0] / "spica_top.json"
    outpdb = "fin.cg.pdb"
    argparser = ArgumentParser(usage='maptraj [-h] [-outpdb OUTPDB] [-json JSON] [-begin BEGIN] [-last LAST] [-nodelwat] inputPDB inputTRAJ outputTRAJ',
                               prog ="maptraj")
    argparser.add_argument('inputAAPDB', type=str,
                            help='Specify input AA PDB file name.')
    argparser.add_argument('inputAATRAJ', type=str,
                            help='Specify input AA TRAJ file name.')
    argparser.add_argument('outputCGTRAJ', type=str,
                            help='Specify output CG TRAJ file name.')
    argparser.add_argument('-outpdb', type=str,
                            default=outpdb,
                            help='output PDB file name (default: fin.cg.pdb).')
    argparser.add_argument('-json', type=str,
                            default=json,
                            help='input json file name (default: spica_top.json).')
    argparser.add_argument('-begin', type=int,
                            default=0,
                            help='begining frame of mappin (default: 0).')
    argparser.add_argument('-last', type=int,
                            default=-1,
                            help='last frame of mapping (default: -1).')
    argparser.add_argument('-nodelwat', action='store_true',
                            help='not delete excess water due to CG ion mapping (default: off).')
    return argparser.parse_args(argv)

class map_traj:
    def __init__(self, inPDB, inTRAJ, outTRAJ, outPDB, jsonfile, nodelwat, beg, last):
        self.outTRAJ = outTRAJ
        self.beg     = beg
        self.last    = last

        self.cg = m2c.map_to_cg(inPDB, outPDB, jsonfile, nodelwat, 0)
        self.cg.mapping()
        cg_info = np.array(self.cg.all_cg_pdb_domain)
        self.aa_pdb0 = np.array(self.cg.pdb_data)
        self.ncg, self.resindices = self._set_residx(cg_info)
        self.u_aa = mda.Universe(inPDB, inTRAJ)
        self.u_cg = mda.Universe.empty(self.ncg, trajectory=True)
        #u_cg = mda.Universe.empty(ncg, n_residues=ncg,
        #                          atom_resindex=resindices,
        #                          residue_segindex=[0]*ncg,
        #                          trajectory=True) 
        #u_cg.add_TopologyAttr('name',    cg_info[:,0])
        #u_cg.add_TopologyAttr('resname', cg_info[:,1])

    def _set_residx(self, cg_info):
        resid_col = 2
        cnt = 0; x0 = None
        resindices = []
        for x in cg_info[:, resid_col]:
            if x != x0:
                cnt += 1
            resindices.append(cnt)
            x0 = x
        return len(resindices), resindices
    
    def run(self):
        pdb_sta = 8
        pdb_end = 11
        crd_col = -3
        cg      = self.cg
        aa_pdb0 = self.aa_pdb0
        u_cg    = self.u_cg
        u_aa    = self.u_aa
        with mda.Writer(self.outTRAJ, self.ncg) as w:
            for ts in tqdm(u_aa.trajectory[self.beg:self.last]):
                aa_pdb0[:, pdb_sta:pdb_end] = u_aa.atoms.positions
                cg.pdb_data     = aa_pdb0
                cg.domain_data  = []
                m2c.make_domain(cg.pdb_data, cg.domain_data)
                cg.initialize()
                cg.mapping()
                u_cg.atoms.positions = np.array(cg.all_cg_pdb_domain)[:, crd_col:]
                u_cg.dimensions = u_aa.dimensions
                w.write(u_cg.atoms)
        cg.finalize()

if __name__ == "__main__":
    args     = get_option()
    inPDB    = args.inputPDB
    inTRAJ   = args.inputTRAJ
    outTRAJ  = args.outputTRAJ
    outPDB   = args.outpdb
    jsonfile = args.json
    nodelwat = args.nodelwat
    beg      = args.begin
    last     = args.last

    mt = map_traj(inPDB, inTRAJ, outTRAJ, outPDB, jsonfile, nodelwat, beg, last)
    mt.run()
