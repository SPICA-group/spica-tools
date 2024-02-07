import sys, json
from pathlib import Path
from argparse import ArgumentParser

def get_option():
    json = Path(__file__).parents[0] / "spica_top.json"
    argparser = ArgumentParser()
    argparser.add_argument('resname', type=str,
                            help='Specify a resname defined in json.')
    argparser.add_argument('-json', type=str,
                            default=json,
                            help='input json file name (default: spica_top.json).')
    argparser.add_argument('-pspica', action="store_true",
                            help='assign pSPICA partial charges.')
    argparser.add_argument('-onec', action="store_true",
                            help='apply a unit charge of 1.0.')
    argparser.add_argument('-dupang', action="store_true",
                            help='not delete duplicated angle indices.')
    return argparser.parse_args()

def get_option_script(argv):
    json = Path(__file__).parents[0] / "spica_top.json"
    argparser = ArgumentParser(usage='json2top [-h] [-json JSON] [-pspica] [-onec] [-dupang] resname',
                               prog ="json2top")
    argparser.add_argument('resname', type=str,
                            help='Specify a resname defined in json.')
    argparser.add_argument('-json', type=str,
                            default=json,
                            help='input json file name (default: spica_top.json).')
    argparser.add_argument('-pspica', action="store_true",
                            help='assign pSPICA partial charges.')
    argparser.add_argument('-onec', action="store_true",
                            help='apply a unit charge of 1.0.')
    argparser.add_argument('-dupang', action="store_true",
                            help='not delete duplicated angle indices.')
    return argparser.parse_args(argv)

class json_to_top:
    def __init__(self, jsonfile, pspica, onec, dupang):
        self.res_type = {}
        self.res_name = {}
        self.res_ch   = {}
        self.res_ms   = {}
        self.bnd_list = {}
        self.jsonfile = jsonfile
        self.pspica   = pspica
        self.onec     = onec
        self.dupang   = dupang
        self.read_json(jsonfile)

    def read_json(self, jsonfile):
        try:
            f = open(jsonfile,"r")
        except:
            sys.exit (f"ERROR: '{jsonfile}' cannot be opened.")
        else:
            jsn = json.load(f)
            f.close()
            for ir in jsn["topo"].keys():
                self.res_type[ir] = jsn["topo"][ir]["type"]
                self.res_name[ir] = jsn["topo"][ir]["name"]
                self.res_ch[ir]   = jsn["topo"][ir]["charge"]
                self.res_ms[ir]   = jsn["topo"][ir]["mass"]
                self.bnd_list[ir] = jsn["topo"][ir]["bonds"]
    
    def name2idx(self, ir):
        d    = {}
        bndx = []
        if ir not in self.bnd_list.keys():
            sys.exit (f"ERROR: Resname '{ir}' is not found in '{self.jsonfile}'.")
        if self.bnd_list[ir][0] == "none":
            print(f"{ir} consists of a particle.")
            return True
        for idx, name in enumerate(self.res_name[ir]):
            d[name] = idx
        self.natom = idx
        for idx in range(len(self.bnd_list[ir])):
            bndx.append([d[x]+1 for x in self.bnd_list[ir][idx]])
        self.bndx = bndx
        return False

    def mod_charge(self, res_ch, res_type):
        unit_spica  = 0.1118
        unit_pspica = 0.5590
        if res_type == "WO":
            if self.onec:
                return -0.8360
            else:
                return -0.4673
        elif res_type == "WH":
            if self.onec:
                return  0.8360
            else:
                return  0.4673
        elif res_type == "SOD1":
            if self.onec:
                return  1.3360
            else:
                return  0.7468
        elif res_type == "SOD2":
            if self.onec:
                return -0.3360
            else:
                return -0.1878
        elif res_type == "CLA1":
            if self.onec:
                return -1.3360
            else:
                return -0.7468
        elif res_type == "CLA2":
            if self.onec:
                return  0.3360
            else:
                return  0.1878
        else:
            if res_ch == 0.0:
                return res_ch
            else:
               rch = round(res_ch, 4)
               q1 = rch == 0
               q2 = rch % unit_spica == 0
               q3 = rch % unit_pspica == 0
               if q1:
                   if self.pspica:
                       return rch * unit_pspica
                   elif self.onec:
                       return rch
                   else:
                       return res_ch * unit_spica
               if q2:
                   if self.pspica:
                       return rch/unit_spica * unit_pspica
                   elif self.onec:
                       return rch/unit_spica
                   else:
                       return res_ch
               if q3:
                   if self.pspica:
                       return rch
                   elif self.onec:
                       return rch/unit_pspica
                   else:
                       return res_ch/unit_pspica * unit_spica
    
    def write_atom(self, ir, f):
        res_type = self.res_type[ir]
        res_name = self.res_name[ir]
        res_ms   = self.res_ms[ir]
        res_ch   = self.res_ch[ir]
        for i, (n, t, m, c) in enumerate(zip(res_name, res_type, res_ms, res_ch)):
            print ("atom %5d %5s %5s %5s  %8.4f  %8.4f  U" % (i+1, ir, n, t, m, self.mod_charge(c, t)), file=f)
        print(file=f)
    
    def write_bond(self, ir, f):
        bndx = self.bndx
        name = self.res_name[ir]
        for i in range(len(bndx)):
            bnd1 = bndx[i][0]
            bnd2 = bndx[i][1]
            print ("bond %5d %5d # %s %s" % (bnd1, bnd2, name[bnd1-1], name[bnd2-1]), file=f)
        print(file=f)
    
    def write_angle(self, ir, f):
        bndx = self.bndx
        name = self.res_name[ir]
        andxs = []
        sets  = []
        for i in range(len(bndx)-1):
            for j in range(i+1, len(bndx)):
                bA = False
                if bndx[i][1] == bndx[j][0]: # i j : j k
                    andx1 = bndx[i][0]
                    andx2 = bndx[i][1]
                    andx3 = bndx[j][1]
                    bA = True
                elif bndx[i][0] == bndx[j][0]: # j i : j k
                    andx1 = bndx[i][1]
                    andx2 = bndx[i][0]
                    andx3 = bndx[j][1]
                    bA = True
                elif bndx[i][1] == bndx[j][1]: # i j : k j
                    andx1 = bndx[i][0]
                    andx2 = bndx[i][1]
                    andx3 = bndx[j][0]
                    bA = True
                elif bndx[i][0] == bndx[j][1]: # j i : k j
                    andx1 = bndx[i][1]
                    andx2 = bndx[i][0]
                    andx3 = bndx[j][0]
                    bA = True
                if bA:
                    andxs.append([andx1, andx2, andx3])
                    sets.append(set([andx1, andx2, andx3]))
        for i, andx in enumerate(andxs):     
            if self.dupang:
                andx1 = andx[0]
                andx2 = andx[1]
                andx3 = andx[2]
                print ("angle %5d %5d %5d # %s %s %s" \
                    % (andx1, andx2, andx3, name[andx1-1], name[andx2-1], name[andx3-1]), file=f)
                continue
            if sets.count(sets[i]) == 1:
                andx1 = andx[0]
                andx2 = andx[1]
                andx3 = andx[2]
                print ("angle %5d %5d %5d # %s %s %s" \
                    % (andx1, andx2, andx3, name[andx1-1], name[andx2-1], name[andx3-1]), file=f)
        print(file=f)
    
    def finalize(self, outfile):
        print(f"{outfile} was generated.")
        print(f"Json file: {self.jsonfile}.")
        print("NOTE: Please make sure that the bond and angle indices are correct")
        print("      in the generated top file, especially for molecules including")
        print("      ring groups.")
    
    def run(self, res, outfile):
        with open(outfile, "w") as f:
            nobnd = self.name2idx(res)
            self.write_atom(res, f)
            if nobnd:
                self.finalize(outfile)
                return
            self.write_bond(res, f)
            if len(self.bndx) == 1:
                self.finalize(outfile)
                return
            self.write_angle(res, f)
            self.finalize(outfile)
            return

if __name__ == "__main__":
    args   = get_option()
    res    = args.resname
    top    = f"{res}.top"
    jfile  = args.json
    pspica = args.pspica
    onec   = args.onec
    da     = args.dupang
    jt = json_to_top(jfile, pspica, onec, da)
    jt.run(res, top)
