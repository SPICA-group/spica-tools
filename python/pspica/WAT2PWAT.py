#!/usr/bin/env python
import sys
from math import sin,cos
import random
import numpy as np

def gen_line(infile):
	f = open(infile,"r")
	n = sum(1 for line in f)
	f.close()
        return n

def read_pdb(infile,pdb_list):
	f = open(infile,"r")
	line = f.readline()
	while line:
	  recname = line[0:6].strip()
	  if recname == "CRYST1":
	    Lx = float(line[6:15])
	    Ly = float(line[15:24])
	    Lz = float(line[24:33])
	    Aa = float(line[33:40])
	    Ab = float(line[40:47])
	    Ag = float(line[47:54])
	    Lstr = line[55:66]
	    Zval = line[66:70].split("\n")[0]
	    pdb_list.append([recname,Lx,Ly,Lz,Aa,Ab,Ag,Lstr,Zval])
	  elif recname == "ATOM":
	    index    = line[6:11]
	    atmname = line[12:16]
	    indicat  = line[16:17]
	    resname =  line[17:21]
	    chainid =  line[21:22]
	    resid   =  int(line[22:26])
	    code    =  line[26:27]
	    posX = float(line[30:38])
	    posY = float(line[38:46])
	    posZ = float(line[46:54])
	    occup = line[54:60].split("\n")[0]
	    Tfact = line[60:66].split("\n")[0]
	    segid = line[72:76].split("\n")[0]
	    elesym = line[76:78].split("\n")[0]
	    #charge = line[78:80]
	    pdb_list.append([recname,index,atmname,indicat,resname,chainid,resid,code,posX,posY,posZ,occup,Tfact,segid,elesym])
	  elif recname == "END":
	    pdb_list.append([recname])
	    return 0
	  else:
	    print "Skip",recname
	  line = f.readline()
	f.close()

def write_pdb(f,pdb_list,l):
	recname = pdb_list[l][0]
	if recname == "CRYST1":
	    Lx = pdb_list[l][1]
	    Ly = pdb_list[l][2]
	    Lz = pdb_list[l][3]
	    Aa = pdb_list[l][4]
	    Ab = pdb_list[l][5]
	    Ag = pdb_list[l][6]
	    Lstr = pdb_list[l][7]
	    Zval = pdb_list[l][8]
	    print >> f,'%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4s' %(recname,Lx,Ly,Lz,Aa,Ab,Ag,Lstr,Zval)
	elif recname == "ATOM":
	    index    = str(pdb_list[l][1])
	    atmname  = pdb_list[l][2]
	    indicat  = pdb_list[l][3]
	    resname  = pdb_list[l][4]
	    chainid  = pdb_list[l][5]
	    resid    = pdb_list[l][6]
	    code     = pdb_list[l][7]
	    posX     = pdb_list[l][8]
	    posY     = pdb_list[l][9]
	    posZ     = pdb_list[l][10]
	    occup    = pdb_list[l][11]
	    Tfact    = pdb_list[l][12]
	    segid    = pdb_list[l][13]
	    elesym   = pdb_list[l][14]
	    #charge  = line[78:80]
	    print >> f,'%-6s%5s %4s%1s%4s%1s%4d%1s   %8.3f%8.3f%8.3f%6s%6s      %-4s%2s' \
                    %(recname,index[0:5],atmname,indicat,resname,chainid,resid,code,posX,posY,posZ,occup,Tfact,segid,elesym)
	elif recname == "END":
	      print >> f,'END'
	else:
	      print "Error.",recname

def write_pdb_tmp(f,pdb_list):
	recname = pdb_list[0]
	if recname == "CRYST1":
	    Lx = pdb_list[1]
	    Ly = pdb_list[2]
	    Lz = pdb_list[3]
	    Aa = pdb_list[4]
	    Ab = pdb_list[5]
	    Ag = pdb_list[6]
	    Lstr = pdb_list[7]
	    Zval = pdb_list[8]
	    print >> f,'%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4s' %(recname,Lx,Ly,Lz,Aa,Ab,Ag,Lstr,Zval)
	elif recname == "ATOM":
	    index    = str(pdb_list[1])
	    atmname  = pdb_list[2]
	    indicat  = pdb_list[3]
	    resname  = pdb_list[4]
	    chainid  = pdb_list[5]
	    resid    = pdb_list[6]
	    code     = pdb_list[7]
	    posX     = pdb_list[8]
	    posY     = pdb_list[9]
	    posZ     = pdb_list[10]
	    occup    = pdb_list[11]
	    Tfact    = pdb_list[12]
	    segid    = pdb_list[13]
	    elesym   = pdb_list[14]
	    #charge  = line[78:80]
	    print >> f,'%-6s%5s %4s%1s%4s%1s%4d%1s   %8.3f%8.3f%8.3f%6s%6s      %-4s%-2s' \
                    %(recname,index[0:5],atmname,indicat,resname,chainid,resid,code,posX,posY,posZ,occup,Tfact,segid,elesym)
	elif recname == "END":
	      print >> f,'END'
	else:
	      print "Error.",recname

def gen_angle(angle_3):
        for i in range(3):
	     angle_3[i] = random.uniform(0,359)

def gen_Rmat(evec,x):
        a = np.array([evec[0],evec[1],evec[2]])
        b = a/np.linalg.norm(a)
        R = np.array([[b[0]*b[0]*(1.0 - cos(x)) + cos(x),\
                      b[0]*b[1]*(1.0 - cos(x)) - b[2]*sin(x),\
                      b[0]*b[2]*(1.0 - cos(x)) + b[1]*sin(x)],\
                     [b[0]*b[1]*(1.0 - cos(x)) + b[2]*sin(x),\
                      b[1]*b[1]*(1.0 - cos(x)) + cos(x),\
                      b[1]*b[2]*(1.0 - cos(x)) - b[0]*sin(x)],\
                     [b[0]*b[2]*(1.0 - cos(x)) - b[1]*sin(x),\
                      b[1]*b[2]*(1.0 - cos(x)) + b[0]*sin(x),\
                      b[2]*b[2]*(1.0 - cos(x)) + cos(x)]])
        return R

def rotate(pos,angle):
        evec1 = [1.,0.,0.]
        evec2 = [0.,1.,0.]
        evec3 = [0.,0.,1.]
        Rx = np.identity(3)
        Ry = np.identity(3)
        Rz = np.identity(3)
        Rx = gen_Rmat(evec1,angle[0])
        Ry = gen_Rmat(evec2,angle[1])
        Rz = gen_Rmat(evec3,angle[2])
        r = np.array([pos[0],pos[1],pos[2]])
        r = np.dot(Rz,r)
        r = np.dot(Ry,r)
        r = np.dot(Rx,r)
        return r

def trl(pos,dx,dy,dz):
	a = np.array([dx,dy,dz])
	pos = pos + a
	return pos

### main ###
args = sys.argv
if len(args) != 3:
	print "USAGE: WAT2PWAT.py input_filename output_filename"
	sys.exit(0)
infile = args[1]
outfile= args[2]
fo = open(outfile,"w")

cnt=1
l = 1.1
half_l = 0.5*1.1
angle_3 = [0.,0.,0.]
pdb_list=[]

nl = gen_line(infile)
read_pdb(infile,pdb_list)
for i in xrange(len(pdb_list)):
      if pdb_list[i][0].find("CRYST1") != -1:
      	pdb_list_tmp = pdb_list[i]
	write_pdb_tmp(fo,pdb_list_tmp)
      elif pdb_list[i][0].find("ATOM") != -1:
	if pdb_list[i][4].find("WAT")!= -1 or pdb_list[i][4].find("TIP3")!= -1:
		pos_n = [0.,0.,0.]
		pos_p = [l,0.,0.]
		dx = pdb_list[i][8]
		dy = pdb_list[i][9]
		dz = pdb_list[i][10]
		gen_angle(angle_3)
		pos_p = rotate(pos_p,angle_3)
		pos_n = trl(pos_n,dx,dy,dz)
		pos_p = trl(pos_p,dx,dy,dz)
		pdb_list_tmp = [pdb_list[i][0],cnt,"W1"," ","PWAT",pdb_list[i][5],pdb_list[i][6]," ",pos_n[0],pos_n[1],pos_n[2],\
				pdb_list[i][11],pdb_list[i][12],pdb_list[i][13],pdb_list[i][14]]
 	        write_pdb_tmp(fo,pdb_list_tmp)
		cnt += 1
		pdb_list_tmp = [pdb_list[i][0],cnt,"W2"," ","PWAT",pdb_list[i][5],pdb_list[i][6]," ",pos_p[0],pos_p[1],pos_p[2],\
				pdb_list[i][11],pdb_list[i][12],pdb_list[i][13],pdb_list[i][14]]
 	        write_pdb_tmp(fo,pdb_list_tmp)
		cnt += 1
	elif pdb_list[i][4].find("SOD") != -1:
		pos_n = [0.,0.,0.]
		pos_p = [l,0.,0.]
		dx = pdb_list[i][8]
		dy = pdb_list[i][9]
		dz = pdb_list[i][10]
		gen_angle(angle_3)
		pos_p = rotate(pos_p,angle_3)
		pos_n = trl(pos_n,dx,dy,dz)
		pos_p = trl(pos_p,dx,dy,dz)
		pdb_list_tmp = [pdb_list[i][0],cnt,"SOD1"," ","PSOD",pdb_list[i][5],pdb_list[i][6]," ",pos_n[0],pos_n[1],pos_n[2],\
				pdb_list[i][11],pdb_list[i][12],pdb_list[i][13],pdb_list[i][14]]
 	        write_pdb_tmp(fo,pdb_list_tmp)
		cnt += 1
		pdb_list_tmp = [pdb_list[i][0],cnt,"SOD2"," ","PSOD",pdb_list[i][5],pdb_list[i][6]," ",pos_p[0],pos_p[1],pos_p[2],\
				pdb_list[i][11],pdb_list[i][12],pdb_list[i][13],pdb_list[i][14]]
 	        write_pdb_tmp(fo,pdb_list_tmp)
		cnt += 1
	elif pdb_list[i][4].find("CLA") != -1:
		pos_n = [0.,0.,0.]
		pos_p = [l,0.,0.]
		dx = pdb_list[i][8]
		dy = pdb_list[i][9]
		dz = pdb_list[i][10]
		gen_angle(angle_3)
		pos_p = rotate(pos_p,angle_3)
		pos_n = trl(pos_n,dx,dy,dz)
		pos_p = trl(pos_p,dx,dy,dz)
		pdb_list_tmp = [pdb_list[i][0],cnt,"CLA1"," ","PCLA",pdb_list[i][5],pdb_list[i][6]," ",pos_n[0],pos_n[1],pos_n[2],\
				pdb_list[i][11],pdb_list[i][12],pdb_list[i][13],pdb_list[i][14]]
 	        write_pdb_tmp(fo,pdb_list_tmp)
		cnt += 1
		pdb_list_tmp = [pdb_list[i][0],cnt,"CLA2"," ","PCLA",pdb_list[i][5],pdb_list[i][6]," ",pos_p[0],pos_p[1],pos_p[2],\
				pdb_list[i][11],pdb_list[i][12],pdb_list[i][13],pdb_list[i][14]]
 	        write_pdb_tmp(fo,pdb_list_tmp)
		cnt += 1
	elif pdb_list[i][4].find("CAL") != -1:
		pos_n = [0.,0.,0.]
		pos_p = [l,0.,0.]
		dx = pdb_list[i][8]
		dy = pdb_list[i][9]
		dz = pdb_list[i][10]
		gen_angle(angle_3)
		pos_p = rotate(pos_p,angle_3)
		pos_n = trl(pos_n,dx,dy,dz)
		pos_p = trl(pos_p,dx,dy,dz)
		pdb_list_tmp = [pdb_list[i][0],cnt,"CAL1"," ","PCAL",pdb_list[i][5],pdb_list[i][6]," ",pos_n[0],pos_n[1],pos_n[2],\
				pdb_list[i][11],pdb_list[i][12],pdb_list[i][13],pdb_list[i][14]]
 	        write_pdb_tmp(fo,pdb_list_tmp)
		cnt += 1
		pdb_list_tmp = [pdb_list[i][0],cnt,"CAL2"," ","PCAL",pdb_list[i][5],pdb_list[i][6]," ",pos_p[0],pos_p[1],pos_p[2],\
				pdb_list[i][11],pdb_list[i][12],pdb_list[i][13],pdb_list[i][14]]
 	        write_pdb_tmp(fo,pdb_list_tmp)
		cnt += 1
	else:
      		pdb_list_tmp = pdb_list[i]
		pdb_list_tmp[1] = cnt
 	        write_pdb_tmp(fo,pdb_list_tmp)
		cnt += 1
      elif pdb_list[i][0] == "END":
      	pdb_list_tmp = pdb_list[i]
 	write_pdb_tmp(fo,pdb_list_tmp)
      else:
	print "Writing Error"	
print "Fin."
fo.close()
