import sys, itertools, re
from argparse import ArgumentParser

def get_option():
    conf = "conf.pdb"
    ndxf = "CGindex.ndx"
    outf = "npt.mdp"
    argparser = ArgumentParser()
    argparser.add_argument('-pdb', type=str,
                            default=conf,
                            help='input PDB file name (default: conf.pdb).')
    argparser.add_argument('-ndx', type=str,
                            default=ndxf,
                            help='input index file name (default: CGindex.ndx).')
    argparser.add_argument('-output', type=str,
                            default=outf,
                            help='specify GROMACS input file name (default: npt.mdp).')
    argparser.add_argument('-T', type=float,
                            default=310.0,
                            help='system temperature [K] (default: 310).')
    argparser.add_argument('-P', type=float,
                            default=1.0,
                            help='system pressure [bar] (default: 1.0).')
    argparser.add_argument('-pspica', action='store_true',
                            help='generate GROMACS input for pSPICA FF (default: for SPICA FF).')
    argparser.add_argument('-em', action='store_true',
                            help='generate EM input for pSPICA FF (default: off).')
    return argparser.parse_args()

def get_option_script(argv):
    conf = "conf.pdb"
    ndxf = "CGindex.ndx"
    outf = "npt.mdp"
    argparser = ArgumentParser(usage='gen_gmxin [-h] [-conf CONF] [-ndxf NDXF] [-output OUTPUT] [-pspica] [-em]',
                               prog ="gen_gmxin")
    argparser.add_argument('-pdb', type=str,
                            default=conf,
                            help='input PDB file name (default: conf.pdb).')
    argparser.add_argument('-ndx', type=str,
                            default=ndxf,
                            help='input index file name (default: CGindex.ndx).')
    argparser.add_argument('-output', type=str,
                            default=outf,
                            help='specify GROMACS input file name (default: npt.mdp).')
    argparser.add_argument('-T', type=float,
                            default=310.0,
                            help='system temperature [K] (default: 310).')
    argparser.add_argument('-P', type=float,
                            default=1.0,
                            help='system pressure [bar] (default: 1.0).')
    argparser.add_argument('-pspica', action='store_true',
                            help='generate GROMACS input for pSPICA FF (default: for SPICA FF).')
    argparser.add_argument('-em', action='store_true',
                            help='generate EM input for pSPICA FF (default: off).')
    return argparser.parse_args(argv)

class gen_gmx_inp:
    def __init__(self, conf, ndxf, outf, temp, press, pspica, em):
        self.conf = conf    
        self.ndxf = ndxf    
        self.outf = outf
        self.temp = temp
        self.press = press
        self.pspica = pspica
        self.em = em

    def gen_table(self, n, m, tname):
       rcut = 2.5
       dr = 0.002
       nbin = int(rcut/dr) + 1
       with open(tname, "w") as f:
           print(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, file=f)
           for i in range(1, nbin):
               r = dr*i
               r0 = round(r,5)
               r1 = round(-1./r**m, 10)
               r2 = round(-m/r**(m+1), 10)
               r3 = round(1./r**n, 10)
               r4 = round(n/r**(n+1), 10)
               print(f"{r0:.5f} 0.0 0.0 {r1:.10e} {r2:.10e} {r3:.10e} {r4:.10e}", file=f)

    def energy_grps(self):
        ref_grp = ["LJ124W", "LJ96W", "SOLW", "PSOLW"]
        grps = []
        with open(self.ndxf, "r") as f:
            for line in f.readlines():
                if line.strip().startswith("["):
                    line = re.sub("\[","", line)
                    line = re.sub("\]","", line)
                    grp = line.strip()
                    if grp in ref_grp:
                        grps.append(grp)       
        return grps
 
    def gen_enegrps(self):
        grps = self.energy_grps()
        s = f"{' '.join(grps)}"
        return s

    def gen_grptable(self):
        grps = self.energy_grps()
        lgrps = []
        self.gen_table(9, 6, "table.xvg")
        for x, y in list(itertools.combinations_with_replacement(grps, 2)):
            if (x == "PSOLW" or y == "PSOLW") or (x == "SOLW" or y == "SOLW"):
                lgrps.append(sorted([x, y]))
        for x, y in lgrps:
            if x == "LJ124W":
                self.gen_table(12, 4, f"table_{x}_{y}.xvg")
            if x == "LJ96W":
                self.gen_table(9, 6, f"table_{x}_{y}.xvg")
            if x == "SOLW":
                self.gen_table(12, 4, f"table_{x}_{y}.xvg")
            if x == "PSOLW":
                self.gen_table(12, 5, f"table_{x}_{y}.xvg")
        ss = [ f"{x} {y}" for x, y in lgrps ]      
        s = f"{' '.join(ss)}"
        return s
  
    def read_box(self):
        with open(self.conf, "r") as f:
            for line in f.readlines():
                items = line.split()
                if "CRYST1" in items:
                    boxl = [float(items[1]), float(items[2]), float(items[3])]
                    return boxl
            sys.exit("ERROR: Box length information is not found.") 
    
    def run(self):
        f = open(self.outf, "w")
        if self.pspica:
            print("##### Generate GROMACS input file for pSPICA #####")
        else:
            print("##### Generate GROMACS input file for SPICA #####")
        print()
        if self.pspica:
            print("; GROMACS input for pSPICA FF", file=f)
        else:
            print("; GROMACS input for SPICA FF", file=f)
        print(file=f)
        print(f"title                = NPT simulation", file=f)
        print(f"pbc                  = xyz", file=f)
        print(f"integrator           = md", file=f)
        print(f"dt                   = 0.01", file=f)
        print(f"nsteps               = 10000000", file=f)
        print(f"nstxtcout            = 10000", file=f)
        print(f"nstxout              = 0", file=f)
        print(f"nstvout              = 0", file=f)
        print(f"nstfout              = 0", file=f)
        print(f"nstlog               = 1000", file=f)
        print(f"nstenergy            = 100", file=f)
        print(f"comm_mode            = linear", file=f)
        print(file=f)
        print(f"vdw-type             = user", file=f)
        print(f"energygrps           = {self.gen_enegrps()}", file=f)
        print(f"energygrp_table      = {self.gen_grptable()}", file=f)
        print(f"cutoff-scheme        = Group", file=f)
        print(file=f)
        print(f"ns_type              = grid", file=f)
        print(f"nstlist              = 5", file=f)
        print(f"rlist                = 1.5", file=f)
        print(f"rcoulomb             = 1.5", file=f)
        print(f"rvdw                 = 1.5", file=f)
        print(file=f)
        print(f"coulombtype          = PME", file=f)
        print(f"pme_order            = 4", file=f)
        if self.pspica:
            print(f"fourierspacing       = 0.2", file=f)
        else:
            print(f"fourierspacing       = 0.5", file=f)
        print(file=f)
        print(f"Tcoupl               = Nose-Hoover", file=f)
        print(f"tau_t                = 1.0", file=f)
        print(f"tc-grps              = system", file=f)
        print(f"ref_t                = {self.temp}", file=f)
        print(file=f)
        print(f"Pcoupl               = Parrinello-Rahman", file=f)
        boxl = self.read_box()
        if len(set(boxl)) == 1:
            cpl = "isotropic"
            print(f"pcoupltype           = {cpl}", file=f)
            print(f"tau_p                = 5.0", file=f)
            print(f"ref_p                = {self.press}", file=f)
            print(f"compressibility      = 4.5e-5", file=f)
        elif boxl[0] == boxl[1] and boxl[0] != boxl[2]:
            cpl = "semiisotropic"
            print(f"pcoupltype           = {cpl}", file=f)
            print(f"tau_p                = 5.0", file=f)
            print(f"ref_p                = {self.press} {self.press}", file=f)
            print(f"compressibility      = 4.5e-5  4.5e-5", file=f)
        else:
            cpl = "anisotropic"
            print(f"pcoupltype           = {cpl}", file=f)
            print(f"tau_p                = 5.0", file=f)
            print(f"ref_p                = {self.press} {self.press} {self.press} {self.press} {self.press} {self.press}", file=f)
            print(f"compressibility      = 4.5e-5  4.5e-5  4.5e-5  0  0  0", file=f)
        print(f"refcoord_scaling     = com", file=f)
        print(file=f)
        if self.pspica:
            print(f"constraint_algorithm = LINCS", file=f)
        print(file=f)
        f.close()
        if self.em:
            self.gen_em()
        print(f"-Pressure coupling is predicted to be '{cpl}' from the simulation box symmetry")
        print()
        if self.pspica:
            print(f"##### GROMACS input file '{self.outf}' for pSPICA has been generated! #####")
        else:
            print(f"##### GROMACS input file '{self.outf}' for SPICA  has been generated! #####")
        print()
        if self.pspica:
            print(f"NOTE: '{self.outf}' is just to run 100 ns GROMACS MD with NPT using pSPICA.")
        else:
            print(f"NOTE: '{self.outf}' is just to run 100 ns GROMACS MD with NPT using SPICA.")
        print(f"NOTE: Please check and modify the input for your simulation purpose.")
        print(f"NOTE: Please add '-table table.xvg' option when excuting 'gmx mdrun'.")

    def gen_em(self):
        f = open("minim.mdp", "w")
        if self.pspica:
            print("; GROMACS input for pSPICA FF", file=f)
        else:
            print("; GROMACS input for SPICA FF", file=f)
        print(file=f)
        print(f"title                = Energy minimization", file=f)
        print(f"pbc                  = xyz", file=f)
        print(f"integrator           = steep", file=f)
        print(f"emtol                = 1000", file=f)
        print(f"nsteps               = 50000", file=f)
        print(file=f)
        print(f"vdw-type             = user", file=f)
        print(f"energygrps           = {self.gen_enegrps()}", file=f)
        print(f"energygrp_table      = {self.gen_grptable()}", file=f)
        print(f"cutoff-scheme        = Group", file=f)
        print(file=f)
        print(f"ns_type              = grid", file=f)
        print(f"nstlist              = 5", file=f)
        print(f"rlist                = 1.5", file=f)
        print(f"rcoulomb             = 1.5", file=f)
        print(f"rvdw                 = 1.5", file=f)
        print(file=f)
        print(f"coulombtype          = PME", file=f)
        print(f"pme_order            = 4", file=f)
        if self.pspica:
            print(f"fourierspacing       = 0.2", file=f)
        else:
            print(f"fourierspacing       = 0.5", file=f)
        f.close()

if __name__ == "__main__":
    args  = get_option()
    conf = args.pdb
    ndxf = args.ndx
    outf = args.output
    temp = args.T
    press = args.P
    pspica = args.pspica
    em = args.em
    obj = gen_lmp_inp(conf, ndxf, outf, temp, press, pspica, em)
    obj.run()
