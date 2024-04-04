import sys
from argparse import ArgumentParser

def get_option():
    data = "DATA.FILE"
    parm = "PARM.FILE"
    outf = "in.lammps"
    argparser = ArgumentParser()
    argparser.add_argument('-data', type=str,
                            default=data,
                            help='input LAMMPS data file name (default: DATA.FILE).')
    argparser.add_argument('-parm', type=str,
                            default=parm,
                            help='input included LAMMPS parameter file name (default: PARM.FILE).')
    argparser.add_argument('-output', type=str,
                            default=outf,
                            help='specify LAMMPS input file name (default: in.lammps).')
    argparser.add_argument('-T', type=float,
                            default=310.,
                            help='system temperature [K] (default: 310.0).')
    argparser.add_argument('-P', type=float,
                            default=1.,
                            help='system pressure [atm] (default: 1.0).')
    argparser.add_argument('-pspica', action='store_true',
                            help='generate LAMMPS input for pSPICA FF (default: for SPICA FF).')
    return argparser.parse_args()

def get_option_script(argv):
    data = "DATA.FILE"
    parm = "PARM.FILE"
    outf = "in.lammps"
    argparser = ArgumentParser(usage='gen_lmpin [-h] [-data DATA] [-parm PARM] [-output OUTPUT] [-pspica]',
                               prog ="gen_lmpin")
    argparser.add_argument('-data', type=str,
                            default=data,
                            help='input LAMMPS data file name (default: DATA.FILE).')
    argparser.add_argument('-parm', type=str,
                            default=parm,
                            help='input included LAMMPS parameter file name (default: PARM.FILE).')
    argparser.add_argument('-output', type=str,
                            default=outf,
                            help='specify LAMMPS input file name (default: in.lammps).')
    argparser.add_argument('-T', type=float,
                            default=310.,
                            help='system temperature [K] (default: 310.0).')
    argparser.add_argument('-P', type=float,
                            default=1.,
                            help='system pressure [atm] (default: 1.0).')
    argparser.add_argument('-pspica', action='store_true',
                            help='generate LAMMPS input for pSPICA FF (default: for SPICA FF).')
    return argparser.parse_args(argv)

class gen_lmp_inp:
    def __init__(self, datafile, parmfile, outf, temp, press, pspica):
        self.datafile = datafile
        self.parmfile = parmfile
        self.outf = outf
        self.temp = temp
        self.press = press
        self.pspica = pspica

    def tune_mesh(self, boxl):
        mesh = []
        if self.pspica:
            mesh_init = list(map(lambda x: int(x//2), boxl))
        else:
            mesh_init = list(map(lambda x: int(x//5), boxl))
        for l in mesh_init:
            while True:
                if self.check_mesh(l):
                    mesh.append(l)
                    break
                else:
                    l += 1
        return mesh

    def pfactor(self, N):
        N_wk = N
        ret = {}
        i = 2
        while i*i <= N:
            while N_wk % i == 0:
                N_wk = N_wk // i
                if i in ret:
                    ret[i] += 1
                else:
                    ret[i] = 1
            i += 1
        if N_wk != 1:
            ret[N_wk] = 1
        return ret

    def check_mesh(self, N):
        lst = {2, 3, 5}
        ret = self.pfactor(N)
        if set(ret.keys()) <= lst:
            return True
        else:
            return False

    def read_polarbond(self):
        polarbond = []
        with open(self.parmfile, "r") as f:
            for line in f.readlines():
                items = line.split()
                if "bond_coeff" in items:
                    b1 = items[-1] == "WO" and items[-2] == "WH"
                    b2 = items[-1] == "WH" and items[-2] == "WO"
                    if b1 or b2:
                        polarbond.append(int(items[1]))
                    b1 = items[-1] == "SOD1" and items[-2] == "SOD2"
                    b2 = items[-1] == "SOD2" and items[-2] == "SOD1"
                    if b1 or b2:
                        polarbond.append(int(items[1]))
                    b1 = items[-1] == "CLA1" and items[-2] == "CLA2"
                    b2 = items[-1] == "CLA2" and items[-2] == "CLA1"
                    if b1 or b2:
                        polarbond.append(int(items[1]))
        s = ""
        if len(polarbond) > 0:
            for x in polarbond:
                s += f"{x} "
            return s
        else:
            return s

    def read_box(self):
        xlo = xhi = ylo = yhi = zlo = zhi = 0.0
        boxl = [0, 0, 0]
        with open(self.datafile, "r") as f:
            for line in f.readlines():
                items = line.split()
                if "xlo" in items and "xhi" in items:
                    xlo = float(items[0]); xhi = float(items[1])
                    boxl[0] = xhi - xlo
                if "ylo" in items and "yhi" in items:
                    ylo = float(items[0]); yhi = float(items[1])
                    boxl[1] = yhi - ylo
                if "zlo" in items and "zhi" in items:
                    zlo = float(items[0]); zhi = float(items[1])
                    boxl[2] = zhi - zlo
                if "Atoms" in items:
                    for x in boxl:
                        if x <= 0.:
                            sys.exit("ERROR: Box lengths include zero or negative values.")
                        return boxl
            sys.exit("ERROR: Box length information is not found.") 
    
    def check_coul(self):
        with open(self.parmfile, "r") as f:
            for line in f.readlines():
                items = line.split()
                if "pair_style" in items:
                    if "coul" in items[1]:
                        return True
                    else:
                        return False
                if "pair_coeff" in items:
                    return False
            return False 
    
    def run(self):
        f = open(self.outf, "w")
        if self.pspica:
            print("##### Generate LAMMPS input file for pSPICA #####")
            print("### LAMMPS input for pSPICA FF ###", file=f)
            gs = 2
        else:
            print("##### Generate LAMMPS input file for SPICA #####")
            print("### LAMMPS input for SPICA FF ###", file=f)
            gs = 5
        print()
        print("# Define variables", file=f)
        print(f"variable       t     equal {round(self.temp, 2)}", file=f)
        print(f"variable       p     equal {round(self.press, 3)}", file=f)
        print(f"variable       ip    index none", file=f)
        print(f"variable       op    index run_00", file=f)
        print(f"variable       min   index no", file=f)
        print(f"variable       dt    index 10", file=f)
        print(f"variable       nstep index 100000", file=f)
        print(f"variable       nlog  index 1000", file=f)
        print(f"variable       ndump index 10000", file=f)
        print(f"variable       nrest index 50000", file=f)
        print(file=f)
        print("# Unit and data file style", file=f)
        print(f"units          real", file=f)
        print(f"atom_style     full", file=f)
        print(file=f)
        print("# Set Newton's third law for interactions", file=f)
        print(f'if "$(is_active(package,gpu))" then &', file=f)
        print(f'"newton         off" &', file=f)
        print(f'else &', file=f)
        print(f'"newton         on off"', file=f)
        print(file=f)
        print("# Load data and parameter files", file=f)
        print('if "${ip} == none" then &', file=f)
        print(f'"read_data      {self.datafile}" &', file=f)
        print(f'else &', file=f)
        print('"read_restart   ${ip}.rest"', file=f)

        print(f"include        {self.parmfile}", file=f)
        print(file=f)
        boxl = self.read_box()
        if self.check_coul():
            print("# Coulomb interaction", file=f)
            print(f"variable       gs equal {gs}", file=f)
            print("variable       nx equal round(lx/${gs})", file=f)
            print("variable       ny equal round(ly/${gs})", file=f)
            print("variable       nz equal round(lz/${gs})", file=f)
            if self.pspica:
                print("kspace_style   pppm 1.0e-5", file=f)
            else:
                print("kspace_style   pppm/cg 1.0e-5", file=f)
            #mesh = self.tune_mesh(boxl)
            if self.pspica:
                print("kspace_modify  mesh ${nx} ${ny} ${nz}", file=f)
            else:
                print("kspace_modify  mesh ${nx} ${ny} ${nz} order 3", file=f)
        print(file=f)
        print("# Neighbor list setup", file=f)
        print(f"neighbor       2.0 bin", file=f)
        print(f"neigh_modify   delay 4", file=f)
        print(file=f)
        print("# Minimization", file=f)
        print('if "${min} == yes" then &', file=f)
        print('"min_style      cg" &', file=f)
        print('"minimize       1.0e-4 1.0e-6 100 1000" &', file=f)
        print('"velocity       all create $t 12345 dist gaussian rot yes"', file=f)
        print(file=f)
        print("# Zero timestep", file=f)
        print(f"reset_timestep 0", file=f)
        print(file=f)
        print("# Timestep for integration [fs]", file=f)
        print("timestep       ${dt}", file=f)
        print(file=f)
        print("# Specify fix for time integration", file=f)
        if self.pspica:
            polarbond = self.read_polarbond()
            if len(polarbond) > 0:
                print(f"fix            0 all shake 1.0e-5 10 0 b {polarbond}", file=f)
        print(f"fix            1 all momentum 1 linear 1 1 1", file=f)
        if len(set(boxl)) == 1:
            print(f"fix            2 all npt temp $t $t 500. iso $p $p 5000.", file=f)
            cpl = "isotropic"
        elif boxl[0] == boxl[1] and boxl[0] != boxl[2]:
            print(f"fix            2 all npt temp $t $t 500. aniso $p $p 5000. couple xy", file=f)
            cpl = "semiisotropic"
        else:
            print(f"fix            2 all npt temp $t $t 500. aniso $p $p 5000.", file=f)
            cpl = "anisotropic"
        print(file=f)
        print("# Output log information", file=f)
        print("thermo         ${nlog}", file=f)
        print(f"thermo_style   custom cpu step temp pxx pyy pzz lx ly lz pe density", file=f)
        print(file=f)
        print("# Output trajectory options", file=f)
        print("dump           1 all xtc ${ndump} ${op}.xtc", file=f)
        print(f"dump_modify    1 unwrap yes", file=f)
        print(file=f)
        print("# Number of steps to generate restart file", file=f)
        print("restart        ${nrest} ${op}.rest ${op}.rest", file=f)
        print(file=f)
        print("# Number of steps to integrate and Perform MD", file=f)
        print("run            ${nstep}", file=f)
        print(file=f)
        print("# Output the final snapshot", file=f)
        print("write_dump     all xtc ${op}.fin.xtc modify unwrap yes", file=f)
        print(file=f)
        f.close()
        print(f"-Pressure coupling is predicted to be '{cpl}' from the simulation box symmetry")
        print()
        if self.pspica:
            print(f"##### LAMMPS input file '{self.outf}' for pSPICA has been generated! #####")
            print(f"NOTE: '{self.outf}' is just to run 1 ns LAMMPS MD with NPT using pSPICA.")
        else:
            print(f"##### LAMMPS input file '{self.outf}' for SPICA has been generated! #####")
            print(f"NOTE: '{self.outf}' is just to run 1 ns LAMMPS MD with NPT using SPICA.")
        print("NOTE: Please check and modify the input for your simulation purpose.")
        print("Command-line Example:")
        print(f"$ lmp_mpi -in {self.outf} -v min yes -v op npt_00 -v nstep 50000 -sf gpu -pk gpu 2")
        print(f"$ lmp_mpi -in {self.outf} -v ip npt_00 -v op npt_01 -v nstep 100000 -v -nlog 5000")

if __name__ == "__main__":
    args  = get_option()
    datafile = args.data
    parmfile = args.parm
    outf = args.output
    temp = args.T     
    press = args.P     
    pspica = args.pspica
    obj = gen_lmp_inp(datafile, parmfile, outf, temp, press, pspica)
    obj.run()
