#!/usr/bin/env python
import sys

charged_list = ["WO","WH","NC","NH","NC4","PH","PHE","CO","SO3","SO4","SOD1","SOD2","CLA1","CLA2","AR2","LY2","ASP","GLU","GBTC","GBTP","GBTN"]
atom_list = []
pair_s   ='pair_style      hybrid/overlay  lj/sdk 15.0 coul/long 15.0 table linear 2999'
bond_s   ='bond_style      harmonic'
angle_s  ='angle_style     hybrid sdk harmonic'
special_b='special_bonds   lj 0.0 0.0 1.0'
WO_WO    ='table 12_5ww.t TABLE 15.0 # WO   WO'
WO_SOD   ='table 12_5ws.t TABLE 15.0 # WO  SOD1'
WO_CLA   ='table 12_5wc.t TABLE 15.0 # WO  CLA1'
SOD_SOD  ='table 12_5ss.t TABLE 15.0 # SOD1 SOD1'
CLA_CLA  ='table 12_5cc.t TABLE 15.0 # CLA1 CLA1'
SOD_CLA  ='table 12_5sc.t TABLE 15.0 # SOD1 CLA1'

args = sys.argv
if len(args) != 3:
	print "USAGE: conver_PARM.py input_filename output_filename"
	sys.exit(0)

inp = args[1]
out = args[2]

f1 = open(inp,"r")
fp = open(out,"w")
print "#atom#"
for line in f1.readlines():
        line = line.strip()
        if line.find("mass") == 0:
            items = line.split()
            print items[4]
            atom_list.append(items[4])
            print >> fp, line
            continue
        elif line.find("pair_style") == 0: 
            print >> fp, pair_s
            continue
        elif line.find("bond_style") == 0: 
            print >> fp, bond_s
            continue
        elif line.find("angle_style") == 0: 
            print >> fp, angle_s
            continue
        elif line.find("special_bonds") == 0: 
            print >> fp, special_b
            continue
        elif line.find("pair_coeff") == 0: 
            items = line.split()
            ni = len(items)
            if atom_list[int(items[1]) - 1]=="WO" and atom_list[int(items[2])- 1]=="WO":
                print >> fp, "%-12s %-8d %-8d %s" % (items[0],int(items[1]),int(items[2]),WO_WO)
            elif atom_list[int(items[1]) - 1]=="SOD1" and atom_list[int(items[2])- 1]=="WO":
                print >> fp, "%-12s %-8d %-8d %s" % (items[0],int(items[1]),int(items[2]),WO_SOD)
            elif atom_list[int(items[1]) - 1]=="WO" and atom_list[int(items[2])- 1]=="SOD1":
                print >> fp, "%-12s %-8d %-8d %s" % (items[0],int(items[1]),int(items[2]),WO_SOD)
            elif atom_list[int(items[1]) - 1]=="CLA1" and atom_list[int(items[2])- 1]=="WO":
                print >> fp, "%-12s %-8d %-8d %s" % (items[0],int(items[1]),int(items[2]),WO_CLA)
            elif atom_list[int(items[1]) - 1]=="WO" and atom_list[int(items[2])- 1]=="CLA1":
                print >> fp, "%-12s %-8d %-8d %s" % (items[0],int(items[1]),int(items[2]),WO_CLA)
            elif atom_list[int(items[1]) - 1]=="SOD1" and atom_list[int(items[2])- 1]=="SOD1":
                print >> fp, "%-12s %-8d %-8d %s" % (items[0],int(items[1]),int(items[2]),SOD_SOD)
            elif atom_list[int(items[1]) - 1]=="CLA1" and atom_list[int(items[2])- 1]=="CLA1":
                print >> fp, "%-12s %-8d %-8d %s" % (items[0],int(items[1]),int(items[2]),CLA_CLA)
            elif atom_list[int(items[1]) - 1]=="SOD1" and atom_list[int(items[2])- 1]=="CLA1":
                print >> fp, "%-12s %-8d %-8d %s" % (items[0],int(items[1]),int(items[2]),SOD_CLA)
            elif atom_list[int(items[1]) - 1]=="CLA1" and atom_list[int(items[2])- 1]=="SOD1":
                print >> fp, "%-12s %-8d %-8d %s" % (items[0],int(items[1]),int(items[2]),SOD_CLA)
            else:
                fp.write("%-12s %-8d %-8d lj/sdk %-6s %8.4f %8.4f %2s %5s %5s\n" % \
                (items[0],int(items[1]),int(items[2]),items[3],float(items[4]),float(items[5]),items[6],items[7],items[8]))
            if atom_list[int(items[1]) - 1 ] in charged_list and atom_list[int(items[2]) - 1 ] in charged_list :
                fp.write("%-12s %-8d %-8d coul/long %2s %5s %5s\n" % \
                (items[0],int(items[1]),int(items[2]),items[6],items[7],items[8]))
            continue
        print >> fp, line

print "Normal termination."
