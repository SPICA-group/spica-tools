#!/usr/bin/env python
# original perl script was made
# BY RUSSELL DEVANE 2/2008
###############################################################################
# READS CG PDB FILE AND DUMPS CONNECTIVITY BASED ON MODEL
# WILL GENERATE ALL BACKBONE-BACKBONE BONDS 
# SIDECHAIN BACKBONE AND SIDECHAIN SIDECAIN BONDS BASED ON MAPPINGS BELOW
# WILL DUMP ANGLES IF WANTED BASED ON ALL BONDS
# WILL DUMP SIDECHAIN IMPROPERS IF WANTED
# WILL DUMP BACKBONE TORSIONS IF WANTED
################################################################################
#
# USAGE: gen_elastic_network.py <cg.pdb filename> < cg.top filename (output)>
#
#################################################################################
# CAN HANDLE MULTIPLE CHAINS
# DOES REVERSE ANGLES NOW
#
#################################################################################
# DOES REQUIRE THE CHEMISTRY MODULE
#################################################################################

#################################################################################
# By Shuhei Kawamoto 3/2012
#   Deleted angle for ring of side chain.
#   Added terminal backbone particle GBT and ABT.
#   Fixed bugs.
#################################################################################
# and converted to PYTHON from PERL
# BY YUSUKE MIYAZAKI 2020/3

import sys
import math

# maximum bond length 
MAXdr = 9.0
# force constant for elastic network model 
kENM  = 1.195  

ncomp = 7
ncomp_pdb = 15
(	PDB_RECNAME,
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
(	ATMNAME,
	RESNAME,
	RESID,
	CHAINID,
	POSX, 
	POSY,
	POSZ,
) = range(ncomp)

res_map={}
res_map["ALA"] = ["ALA"];
res_map["VAL"] = ["VAL"];
res_map["LEU"] = ["LEU"];
res_map["ILE"] = ["ILE"];
res_map["PRO"] = ["PRO"];
res_map["MET"] = ["MET"];
res_map["SER"] = ["SER"];
res_map["THR"] = ["THR"];
res_map["CYS"] = ["CYS"];
res_map["ASN"] = ["ASN"];
res_map["GLN"] = ["GLN"];
res_map["ASP"] = ["ASP"];
res_map["GLU"] = ["GLU"];
res_map["ARG"] = ["AR1", "AR2"];
res_map["LYS"] = ["LY1", "LY2"];
res_map["PHE"] = ["PH1", "PH2", "PH3", "PH4"];
res_map["TYR"] = ["TY1", "TY2", "TY3", "TY4"];
res_map["TRP"] = ["TR1", "TR2", "TR3"];
res_map["HIS"] = ["HI1", "HI2", "HI3"];
res_map["HSP"] = ["HI1", "HI2", "HI3"];
res_map["HSD"] = ["HI1", "HI2", "HI3"];
res_map["HSE"] = ["HI1", "HI2", "HI3"];

this_type = {}
this_type["GBB"] = "GBB";
this_type["GBM"] = "GBM";
this_type["ABB"] = "ABB";
this_type["VAL"] = "VAL";
this_type["LEU"] = "LEU";
this_type["ILE"] = "ILE";
this_type["PRO"] = "PRO";
this_type["MET"] = "MET";
this_type["PH1"] = "XYR";
this_type["PH2"] = "BER";
this_type["PH3"] = "BER";
this_type["PH4"] = "BER";
this_type["TR1"] = "TR1";
this_type["TR2"] = "TR2";
this_type["TR3"] = "TR3";
this_type["SER"] = "SER";
this_type["THR"] = "THR";
this_type["CYS"] = "CYS";
this_type["ASN"] = "ASN";
this_type["GLN"] = "GLN";
this_type["TY1"] = "XYR";
this_type["TY2"] = "BER";
this_type["TY3"] = "PHR";
this_type["TY4"] = "BER";
this_type["ASP"] = "ASP";
this_type["GLU"] = "GLU";
this_type["LY1"] = "LY1";
this_type["LY2"] = "LY2";
this_type["AR1"] = "AR1";
this_type["AR2"] = "AR2";
this_type["HI1"] = "HI1";
this_type["HI2"] = "HI2";
this_type["HI3"] = "HI3";
this_type["GBT"] = "GBT";
this_type["ABT"] = "ABT";

###############################################################################
#  IDENTIFY UNIQUE SIDECHAIN BONDS...I.E. RINGS
#  OR MULTIPLE SITES..
#  IF NOTHING HERE IT ASSUMES STRAIGHTFORWARD BONDING
######################################################
spec_bond={}

lysbnd1 = ["GBM", "LY1"]
lysbnd2 = ["LY1", "LY2"]
spec_bond["LYS"] = [lysbnd1, lysbnd2]

argbnd1 = ["GBM", "AR1"]
argbnd2 = ["AR1", "AR2"]
spec_bond["ARG"] = [argbnd1, argbnd2]

hisbnd1 = ["GBM", "HI1"]
hisbnd2 = ["HI1", "HI2"]
hisbnd3 = ["HI1", "HI3"]
hisbnd4 = ["HI2", "HI3"]
spec_bond["HIS"] = [hisbnd1, hisbnd2, hisbnd3, hisbnd4];
spec_bond["HSP"] = [hisbnd1, hisbnd2, hisbnd3, hisbnd4];
spec_bond["HSD"] = [hisbnd1, hisbnd2, hisbnd3, hisbnd4];
spec_bond["HSE"] = [hisbnd1, hisbnd2, hisbnd3, hisbnd4];

trpbnd1 = ["GBM", "TR1"]
trpbnd2 = ["TR1", "TR2"]
trpbnd3 = ["TR1", "TR3"]
trpbnd4 = ["TR2", "TR3"]
spec_bond["TRP"] = [trpbnd1, trpbnd2, trpbnd3, trpbnd4];

phebnd1 = ["GBM", "PH1"];
phebnd2 = ["PH1", "PH2"];
phebnd3 = ["PH2", "PH3"];
phebnd4 = ["PH3", "PH4"];
phebnd5 = ["PH4", "PH1"];
phebnd6 = ["PH2", "PH4"];
spec_bond["PHE"] = [phebnd1, phebnd2, phebnd3, phebnd4, phebnd5];

tyrbnd1 = ["GBM", "TY1"];
tyrbnd2 = ["TY1", "TY2"];
tyrbnd3 = ["TY2", "TY3"];
tyrbnd4 = ["TY3", "TY4"];
tyrbnd5 = ["TY4", "TY1"];
tyrbnd6 = ["TY2", "TY4"];
spec_bond["TYR"] = [tyrbnd1, tyrbnd2, tyrbnd3, tyrbnd4, tyrbnd5];


###############################################################################
#  IDENTIFY RINGS DIHEDRALS
#  THEY CAN EASILY BE SPECIFIED HERE
######################################################
spec_dihed={}

# THIS IS NOT FULLY IMPLEMENTED.  THE CODE JUST ASSUMES
# THE FOUR CG SITE ARE PART OF THE DIHEDRAL...FINE FOR NOW
# YOU STILL NEED THESE LINE!!!!
spec_dihed["PHE"] = ["PH1", "PH2", "PH3", "PH4"];
spec_dihed["TYR"] = ["TY1", "TY2", "TY3", "TY4"];


#######################################################

# CHARGE HASH
charge={}
charge["LY2"] = 0.1118
charge["AR2"] = 0.1118
charge["HI1"] = 0.00
charge["HI2"] = 0.00
charge["ASP"] = -0.1118
charge["GLU"] = -0.1118
charge["GBT"] = 0.1118
charge["ABT"] = 0.1118
# MASS HASH
bd_mass={}
bd_mass["GBB"] = 56.0385
bd_mass["GBM"] = 56.0385
bd_mass["GBT"] = 56.0385
bd_mass["ABB"] = 71.0732 
bd_mass["ABT"] = 71.0732 
bd_mass["VAL"] = 43.0883
bd_mass["LEU"] = 57.1151
bd_mass["ILE"] = 57.1151
bd_mass["PRO"] = 42.0804
bd_mass["MET"] = 75.1543

bd_mass["TR1"] = 39.0567
bd_mass["TR2"] = 38.0488
bd_mass["TR3"] = 53.0634

bd_mass["SER"] = 31.0287
bd_mass["THR"] = 45.0555
bd_mass["CYS"] = 46.0928
bd_mass["ASN"] = 58.0543
bd_mass["GLN"] = 72.0811

bd_mass["PH1"] = 33.5575
bd_mass["PH2"] = 19.5275 
bd_mass["PH3"] = 19.5275
bd_mass["PH4"] = 19.5275 

bd_mass["TY1"] = 33.5575
bd_mass["TY2"] = 19.5275
bd_mass["TY3"] = 35.5275
bd_mass["TY4"] = 19.5275

bd_mass["ASP"] = 58.0258
bd_mass["GLU"] = 72.0526
bd_mass["LY1"] = 42.0804
bd_mass["LY2"] = 31.0572
bd_mass["AR1"] = 42.0804
bd_mass["AR2"] = 59.0706

bd_mass["HI1"] = 26.0378
bd_mass["HI2"] = 27.0256
bd_mass["HI3"] = 28.0335

def read_pdb(infile,pdb_list,cryst):
	natom = 0
	try:
		f = open(infile,"r")
	except:
		print "ERROR: FILE",infile,"IS NOT FOUND"
		sys.exit(0)
	line = f.readline()
	while line:
	  recname = line[0:6].strip()
	  if recname == "CRYST1":
	    cryst.append(line.strip())
	  elif recname == "ATOM":
	    natom += 1
	    index    = line[6:11]
	    atmname = line[12:16]
	    indicat  = line[16:17]
	    resname =  line[17:21]
	    chainid =  line[21:22]
	    resid   =  line[22:26]
	    code    =  line[26:27]
	    posX = float(line[30:38])
	    posY = float(line[38:46])
	    posZ = float(line[46:54])
	    occup = line[54:60]
	    Tfact = line[60:66]
	    segid = line[72:76]
	    elesym = line[76:78].split("\n")[0]
	    #charge = line[78:80]
	    pdb_list.append([recname,index,atmname,indicat,resname,chainid,resid,code,posX,posY,posZ,occup,Tfact,segid,elesym])
	  line = f.readline()
	f.close()
	return natom

########################################################
# MAIN CODE STARTS HERE..... THE PDB FILE WAS READ ABOVE
# START BY LOOPING OVER EACH RESIDUE (DOMAIN)
########################################################
if __name__ == "__main__":
	args = sys.argv
	if len(args) != 3:
		print "USAGE: gen_elastic_network.py <cg.pdb filename> <cg.top filename (output)>"
		sys.exit(0)
	infile  = args[1]
	outfile = args[2]
	cryst             = []
	pdb_data          = []
    	natom             = read_pdb(infile,pdb_data,cryst)
	nat = 0
	nbb = 0
	bbndx   = []
	iat2ibb = []
	bbtype  = []
	bBackbone = []
	bPH1TY1   = []
	name    = []
	resname = []
	resid = []
	x = []
	y = []
	z = []
	
    	for i in range(natom):
      		tmp_name    = pdb_data[i][PDB_ATMNAME].strip()[:2]
      		tmp_resname = pdb_data[i][PDB_RESNAME].strip()
      
      		if tmp_name == "GB" or tmp_name == "AB":
        		if tmp_resname == "ALA":
              			name.append("ABB")
          		elif tmp_resname == "GLY":
              			name.append("GBB")
          		else :
              			name.append("GBM")
          

         		bbtype.append(pdb_data[i][PDB_RESNAME].strip())
         		bbndx.append(nat)

         		x.append(pdb_data[i][PDB_POSX])
         		y.append(pdb_data[i][PDB_POSY])
         		z.append(pdb_data[i][PDB_POSZ])
         		resname.append(tmp_resname)
         		resid.append(float(pdb_data[i][PDB_RESID]))
         		iat2ibb.append(nbb)
         		bBackbone.append(1)
         		bPH1TY1.append(0)
        		nat+=1; nbb+=1

         		if tmp_resname != "GLY" and tmp_resname != "ALA":
	  			for j in res_map[tmp_resname]:
					name.append(j)
         				x.append(None)
         				y.append(None)
         				z.append(None)
					resname.append(tmp_resname)
         				resid.append(float(pdb_data[i][PDB_RESID]))
         				iat2ibb.append(None)
					bBackbone.append(0)
					if j == "PH2" or j == "TY2" or j == "PH4" or j == "TY4":
						bPH1TY1.append(1)
					else :
						bPH1TY1.append(0)
					nat+=1

######################################################
#########        PRINT OUT TOP FILE      #############
######################################################
    	try:
		ftop = open(outfile,"w")
	except:
		print "ERROR: FILE",outfile,"CANNOT OPEN"
		sys.exit(0)
    	ndx = 0

    	for i in range(nat):
		I = i+1
		if name[i] in charge:
			thischrg = charge[name[i]]
		else :
			thischrg = 0.0
		
		wtype = this_type[name[i]]
		mass  = bd_mass[name[i]]

		if bbndx[0] == i:
			# backbone N-terminal
			thischrg = charge["GBT"]
			#print "charge of terminal is",thischrg

			if resname[i] == "ALA":
				print >> ftop,"atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
					%(I,resname[i],"ABT","ABT",mass,thischrg)
			else :
				print >> ftop,"atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
					%(I,resname[i],"GBT","GBT",mass,thischrg)
			
		elif bbndx[nbb-1] == i:
			# backbone C-terminal
			thischrg = -1*charge["GBT"]
			
			if resname[i] == "ALA":
				print >> ftop,"atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
					%(I,resname[i],"ABT","ABT",mass,thischrg)
			else :
				print >> ftop,"atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
					%(I,resname[i],"GBT","GBT",mass,thischrg)
				
		else :
			print >> ftop,"atom %5d %5s %5s %5s %8.4f   %8.4f  P" \
				%(I,resname[i],name[i],wtype,mass,thischrg)
			
	print >> ftop

##########################
# SET UP BONDS NOW
##########################

#######################################
# Elastic network model for Backbone
########################################
	print 
	print "******** Elastic network info ********"
	print "Cutoff distance ... ",MAXdr,"A"
	print "Force constant ... ",kENM,"kcal/A^2"
	print "**************************************"
	print

	for i1 in range(nat): 
	        I1 = i1 + 1
		for i2 in range(i1 + 1,nat):
		        I2 = i2 + 1
			if bBackbone[i1] == 1 and bBackbone[i2] == 1:
				dx = x[i1] - x[i2]
				dy = y[i1] - y[i2]
				dz = z[i1] - z[i2]
				dr = math.sqrt(dx*dx + dy*dy + dz*dz)
				if dr < MAXdr and resid[i2]-resid[i1] >= 3:
					print >> ftop, "bondparam %5d %5d   %f %f # %s-%s" \
						%(I1,I2,kENM,dr,name[i1],name[i2])
	
	print >> ftop

#######################################
# PRINT THE BACKBONE BONDS FIRST 
#######################################
	total_bonds = 0
	bond_index1 = []
	bond_index2 = []

	for i in range(nbb - 1):
		total_bonds+=1
		print >> ftop, "bond %5d %5d # %s-%s" \
			%(bbndx[i]+1,bbndx[i+1]+1,name[bbndx[i]],name[bbndx[i+1]])
		bond_index1.append(bbndx[i])
		bond_index2.append(bbndx[i+1])

	print >> ftop

#######################################
# DO SIDECHAINS. SPEC_BOND ABOVE CAN BE USED 
# TO DESCRIBE SPECIAL BONDS THAT WOULD BE
# DIFFICULT TO GUESS 
#######################################
	for i in range(nbb):
		if bbtype[i] != "GLY" and bbtype[i] != "ALA":
			# check to see if it is a special bond
			if bbtype[i] in spec_bond:
				for j in spec_bond[bbtype[i]]:
					bndx1 = 0
					bndx2 = 0
					for k in j:
						if bndx1 == 0:
							tmp_ndx = bbndx[i]
							while name[tmp_ndx] != k:
								tmp_ndx+=1
							bndx1 = tmp_ndx
						else :
							tmp_ndx = bbndx[i]
							while name[tmp_ndx] != k:
								tmp_ndx+=1
							bndx2 = tmp_ndx
					print >> ftop, "bond %5d %5d # %s-%s" \
						%(bndx1+1,bndx2+1,name[bndx1],name[bndx2])
					total_bonds+=1
					bond_index1.append(bndx1)
					bond_index2.append(bndx2)

			else :
				print >> ftop, "bond %5d %5d # %s-%s" \
					%(bbndx[i]+1,bbndx[i]+1+1,name[bbndx[i]],name[bbndx[i]+1])
				total_bonds+=1
				bond_index1.append(bbndx[i])
				bond_index2.append(bbndx[i]+1)

	print >> ftop

#########################################
# NOW TAKE CARE OF ANGLES IF WANTED
# GENERATE EVERY COMBINATION FROM THE
# BOND LIST UNLESS SOME ARE NOT SPECIFIED 
# ON THE COMMAND LINE.
#########################################

	for i1 in range(total_bonds):
		for i2 in range(i1 + 1,total_bonds):
			if bond_index1[i1] == bond_index1[i2]:
				andx1 = bond_index2[i1]
				andx2 = bond_index1[i1]
				andx3 = bond_index2[i2]
		      		# GB-GB-GB is non zero.
		      		# GB-GB-SC is non zero.
		      		# GB-SC-SC is non zero.
		      		# SC-TY1-SC is non zero.
		      		# SC-PH1-SC is non zero.
		      		# SC-SC-SC is zero.
				if bBackbone[andx1] + bBackbone[andx2] + bBackbone[andx3] != 0 or bPH1TY1[andx2] == 1:
					print >> ftop, "angle %5d %5d %5d # %s %s %s" \
						%(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3])
			elif bond_index1[i1] == bond_index2[i2]:
				andx1 = bond_index2[i1]
				andx2 = bond_index1[i1]
				andx3 = bond_index1[i2]
				if bBackbone[andx1] + bBackbone[andx2] + bBackbone[andx3] != 0 or bPH1TY1[andx2] == 1:
					print >> ftop, "angle %5d %5d %5d # %s %s %s" \
						%(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3])
			elif bond_index2[i1] == bond_index1[i2]:
				andx1 = bond_index1[i1]
				andx2 = bond_index2[i1]
				andx3 = bond_index2[i2]
				if bBackbone[andx1] + bBackbone[andx2] + bBackbone[andx3] != 0 or bPH1TY1[andx2] == 1:
					print >> ftop, "angle %5d %5d %5d # %s %s %s" \
						%(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3])
			elif bond_index2[i1] == bond_index2[i2]:
				andx1 = bond_index1[i1]
				andx2 = bond_index2[i1]
				andx3 = bond_index1[i2]
				if bBackbone[andx1] + bBackbone[andx2] + bBackbone[andx3] != 0 or bPH1TY1[andx2] == 1:
					print >> ftop, "angle %5d %5d %5d # %s %s %s" \
						%(andx1+1,andx2+1,andx3+1,name[andx1],name[andx2],name[andx3])

	print >> ftop

######################################
#  WRITE OUT ANY DIHEDRALS FOR RINGS
######################################
	for i in range(nbb):
		#MAKE SURE WE HAVE A SIDECHAIN
		if bbtype[i] != "GLY" and bbtype[i] != "ALA":
		#CHECK TO SEE IF IT IS A SPECIAL DIHEDRAL
			if bbtype[i] in spec_dihed:
				ndx1 = bbndx[i] + 1
				ndx2 = ndx1 + 1
				ndx3 = ndx2 + 1
				ndx4 = ndx3 + 1
				print >> ftop, "dihedralparam  %5d %5d %5d %5d   50.0 1 180 0.0 # %s-%s-%s-%s" \
					%(ndx1+1,ndx2+1,ndx3+1,ndx4+1,name[bbndx[i]+1],name[bbndx[i]+2],name[bbndx[i]+3],name[bbndx[i]+4])

	ftop.close()
	print "Normal Termination."
	print
