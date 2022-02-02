/*
	Apr-10-2019 setup_lammps.cpp 
		Modified by Sangjae Seo
			Combined with previous versions.
			Added feature of writing PSF 

	Aug-25-2017 setup_lammps.cpp 
		Modified by Sangjae Seo
			Fix typo for dihedral. 
			Distinguish SDK and harmonic angle type
	

	Mar-20-2014 setup_lammps.cpp
		Modified by Shuhei Kawamoto.
		Added -1 value for bond and angle parameter for protein.
			(When parameter value is -1 for bond length or angle degree, 
			 the values are extracted from pdb file.)
		Fixed Bugs of memory allocation and pointers.
		
	2010?	setup_lammps-wc.c
		Created by Russell DeVane.
*/
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<stdarg.h>
#include<math.h>

const int MaxNBondPerAtom = 20;
const double PI = 3.141592653589793;

void printf1(const int id, const char* format, ...){
	{ // if id is found in id list, skip printf
		static int *ids = 0;
		static int n = 0;
		for (int i=0; i<n; i++){
			if (ids[i]==id) return;
		}
		// if id is not found in id list, add the new id to id list
		int *idsnew = new int [n+1];
		for (int i=0; i<n; i++){
			idsnew[i] = ids[i];
		}
		idsnew[n] = id;
		n++;
		if (ids) delete [] ids;
		ids = idsnew;
	}
	{ // print
		va_list arg;
		va_start(arg, format);
		vprintf(format, arg);
		va_end(arg);
	}
}
double abs(const double a[3]){
	return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}
double inner(const double a[3], const double b[3]){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
void outer(double axb[3], const double a[3], const double b[3]){
	for (int d=0; d<3; d++){
		int d1=(d+1)%3;
		int d2=(d+2)%3;
		axb[d] = a[d1]*b[d2]-b[d1]*a[d2];
	}
}
double GetAngle(const double *r1, const double *r2, const double *r3){
	const double a[3] = {r1[0]-r2[0], r1[1]-r2[1], r1[2]-r2[2]};
	const double b[3] = {r3[0]-r2[0], r3[1]-r2[1], r3[2]-r2[2]};
	const double ab = inner(a,b);
	const double lallbl = abs(a)*abs(b);
	const double cos = ab/lallbl;
	if (cos<-1.0)	  return PI;
	else if (cos>1.0) return 0;
	else			  return acos(cos);
}
double GetDihedral(const double *r1, const double *r2, const double *r3, const double *r4){
	const double r21[3] = {r1[0]-r2[0], r1[1]-r2[1], r1[2]-r2[2]};
	const double r23[3] = {r3[0]-r2[0], r3[1]-r2[1], r3[2]-r2[2]};
	const double r34[3] = {r4[0]-r3[0], r4[1]-r3[1], r4[2]-r3[2]};
	const double r32[3] = {r2[0]-r3[0], r2[1]-r3[1], r2[2]-r3[2]};
	const double cos21 = inner(r21,r23)/(abs(r21)*abs(r23));
	const double cos34 = inner(r34,r32)/(abs(r34)*abs(r32));
	double p1[3], p4[3];
	for (int d=0; d<3; d++){
		p1[d] = r21[d]-cos21*r23[d]/abs(r23)*abs(r21);
		p4[d] = r34[d]-cos34*r32[d]/abs(r32)*abs(r34);
	}
	const double cos = inner(p1,p4)/(abs(p1)*abs(p4));
	double p1xp4[3], sign;
	outer(p1xp4,p1,p4);
	if (inner(p1xp4,r23)>0.0) sign =  1;
	else                      sign = -1; 
	if (cos<-1.0)	  return sign*PI;
	else if (cos>1.0) return sign*0;
	else			  return sign*acos(cos);
}
typedef struct sysdat{
	int nats,nbnds,nangs,nimprops,ndiheds,ntops,total_ats,total_bnds,total_angs,total_improps,total_diheds,foundatoms,boxinfo;
	int uniq_nats,uniq_nbnds,uniq_nangs, uniq_nimprops,uniq_ndiheds;
	int *param_bnds,*param_angs;
	double *coordx,*coordy,*coordz;
	double boxx,boxy,boxz;
} SYSDAT;
typedef struct topdat{
	int *bndndx1,*bndndx2,*bndtype,*angndx1,*angndx2,*angndx3,*angtype;
	int *improp_func,*impropndx1,*impropndx2,*impropndx3,*impropndx4, *improptype;
	int *dihed_func,*dihedndx1,*dihedndx2,*dihedndx3,*dihedndx4, *dihedtype, *dihedn;
	int *dihedpset, *improppset, *bndpset, *angpset;
	int *index, nat,nbnd,nang, nimprop, nmol,*parm_atomtype,ndihed, *dihedeq;
	double *dihedfk, *dihedof;
	double *mass,*charge, *bndfk, *bndeq, *angfk, *angeq, *impropfk, *impropeq;
	char (*atomname)[5],(*type)[5],(*segid)[5],(*resname)[5];
}TOPDAT;
typedef struct database{
	double *fbnd,*bnde,*fang,*ange,*eps,*sig, *angsdk;
	int nvdwtype,nbndtype,nangtype;
	char (*vdwtype1)[5],(*vdwtype2)[5],(*vdwstyle)[7];
	char (*bndtype1)[5],(*bndtype2)[5];
	char (*angtype1)[5],(*angtype2)[5],(*angtype3)[5];
}DATABASE;

void read_pdb(char *filename,SYSDAT *sysdat){
	FILE *fpin;
	const int Nbuf=100;
	char buf[Nbuf];
	char numstr[10];
	sysdat->foundatoms=0;
	sysdat->boxx=0.0;

	if((fpin = fopen(filename,"r")) == NULL){
		fprintf(stderr,"ERROR: can't open infile %s\n",filename);
		exit(1);
	}

	while(fgets(buf,Nbuf,fpin)) {
		if (!strncmp(buf,"CRYST1",6)){
			printf1(__LINE__,"FOUND BOXSIZE DATA\n");
			sscanf(&buf[7],"%lf	%lf	%lf", &sysdat->boxx,&sysdat->boxy,&sysdat->boxz);
		}
		else if (!strncmp(buf,"ATOM ",5) || !strncmp(buf,"HETATM",6)){
			if (sysdat->foundatoms >= sysdat->total_ats){
				printf1(__LINE__,"ERROR: foundatoms in pdb file >= total_ats in top files.\n");
				exit(1);
			}
			memset(numstr, 0, sizeof(numstr));
			strncpy(numstr, buf+30, 8);		sysdat->coordx[sysdat->foundatoms] = atof(numstr);
			strncpy(numstr, buf+38, 8);		sysdat->coordy[sysdat->foundatoms] = atof(numstr);
			strncpy(numstr, buf+46, 8);		sysdat->coordz[sysdat->foundatoms] = atof(numstr);
			sysdat->foundatoms++;
		}
	}
	if(sysdat->foundatoms==0){
		printf1(__LINE__,"ERROR: DID NOT FIND ANY ATOMS IN THE PDB FILE\n");
		exit(1);
	}
	if(sysdat->boxx==0.0){
		printf1(__LINE__,"WARNING: DID NOT FIND CELL SIZE!!\n");
		printf1(__LINE__,"BOX SIZE WILL HAVE BE SET BY HAND\n");
	}

}

void read_coords(char *filename,DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat, int ISCHARGED){
	FILE *fpout;
	int i,j,k,atindex,molindex,bondindex,angleindex,offset,impropindex,dihedindex;

	/* will call modules depending on file type to be read */
	/* right now only pdb capability */

	if(ISCHARGED)
		printf1(__LINE__,"FOUND CHARGES!!\n");
	else
		printf1(__LINE__,"DID NOT FIND CHARGES!!\n");


	if((fpout = fopen("DATA.FILE","w")) == NULL){
		fprintf(stderr,"ERROR: can't open DATA.FILE\n");
		exit(1);
	}
	/* call for reading pdb */
//	read_pdb(filename,sysdat);
//	printf1(__LINE__,"READ PDB FILE\n");

	if(sysdat->total_ats!=sysdat->foundatoms){
		printf1(__LINE__,"ERROR: NUMBER OF ATOMS READ FROM TOPOLOGY AND COMMAND LINE\n");
		printf1(__LINE__,"DOES NOT MATCH THE NUMBER OF ATOMS FOUND IN THE COORD FILE\n");
		exit(0);
	}

	fprintf(fpout,"LAMMPS description\n");
	fprintf(fpout,"\n");
	fprintf(fpout,"%-7d atoms\n",sysdat->total_ats);

	if(sysdat->nbnds>0)
		fprintf(fpout,"%-7d bonds\n",sysdat->total_bnds);

	if(sysdat->nangs>0)
		fprintf(fpout,"%-7d angles\n",sysdat->total_angs);

	if(sysdat->total_diheds>0)
		fprintf(fpout,"%-7d dihedrals\n",sysdat->total_diheds);

	if(sysdat->total_improps>0)
		fprintf(fpout,"%-7d impropers\n",sysdat->total_improps);

	fprintf(fpout,"\n");

	fprintf(fpout,"%-7d atom types\n",sysdat->uniq_nats);

	if(sysdat->nbnds>0)
		fprintf(fpout,"%-7d bond types\n",sysdat->uniq_nbnds);

	if(sysdat->nangs>0)
		fprintf(fpout,"%-7d angle types\n",sysdat->uniq_nangs);

	if(sysdat->total_diheds>0)
		fprintf(fpout,"%-7d dihedral types\n", sysdat->uniq_ndiheds);

	if(sysdat->total_improps>0)
		fprintf(fpout,"%-7d improper types\n", sysdat->uniq_nimprops);

	fprintf(fpout,"\n");

	fprintf(fpout,"%lf %lf xlo xhi\n",-sysdat->boxx/2,sysdat->boxx/2);
	fprintf(fpout,"%lf %lf ylo yhi\n",-sysdat->boxy/2,sysdat->boxy/2);
	fprintf(fpout,"%lf %lf zlo zhi\n",-sysdat->boxz/2,sysdat->boxz/2);

	fprintf(fpout,"\n");
	fprintf(fpout,"Atoms\n");
	fprintf(fpout,"\n");
	atindex=0;
	molindex=0;

	for(i=0;i<sysdat->ntops;i++) {
		for(j=0;j<topdat[i].nmol;j++){
			molindex++;
			for(k=0;k<topdat[i].nat;k++){
				atindex++;

//				if(ISCHARGED){
					/*           if(sysdat->total_bnds>0){ */
					fprintf(fpout,"%d %d %d %7.4f %7.4f %7.4f %7.4f # %s\n",atindex,molindex,topdat[i].parm_atomtype[k]+1,topdat[i].charge[k] ,sysdat->coordx[atindex-1],
						sysdat->coordy[atindex-1],sysdat->coordz[atindex-1],topdat[i].type[k]);
					/*           } */
					/*           else{ */
					/*             fprintf(fpout,"%d %d %7.4f %7.4f %7.4f %7.4f # %s\n",atindex,topdat[i].parm_atomtype[k]+1,topdat[i].charge[k],sysdat->coordx[atindex-1], */
					/*                sysdat->coordy[atindex-1],sysdat->coordz[atindex-1],topdat[i].type[k]); */
					/*           } */
//				} else {
//					if(sysdat->total_bnds>0){
//						fprintf(fpout,"%d %d %d 0.0 %7.4f %7.4f %7.4f # %s\n",atindex,molindex,topdat[i].parm_atomtype[k]+1 ,sysdat->coordx[atindex-1],
//							sysdat->coordy[atindex-1],sysdat->coordz[atindex-1],topdat[i].type[k]);
//					}
//					else{
//						fprintf(fpout,"%d %d %7.4f %7.4f %7.4f # %s\n",atindex,topdat[i].parm_atomtype[k]+1,sysdat->coordx[atindex-1],
//							sysdat->coordy[atindex-1],sysdat->coordz[atindex-1],topdat[i].type[k]);
//					}
//				}
			}
		}
	}

	if(sysdat->total_bnds>0){
		fprintf(fpout,"\n");
		fprintf(fpout,"Bonds\n");
		fprintf(fpout,"\n");
		bondindex=0;
		offset=0;
		for(i=0;i<sysdat->ntops;i++) {
			for(j=0;j<topdat[i].nmol;j++){
				molindex++;
				for(k=0;k<topdat[i].nbnd;k++){
					bondindex++;
					if (bondindex == 3547)
						k *= 1;//debug
					fprintf(fpout,"%d %d %d %d # %s %s\n",bondindex,topdat[i].bndtype[k]+1,
						topdat[i].bndndx1[k]+(j*topdat[i].nat)+offset,topdat[i].bndndx2[k]+(j*topdat[i].nat)+offset,
						topdat[i].type[topdat[i].bndndx1[k]-1],topdat[i].type[topdat[i].bndndx2[k]-1]);
				}
			}
			offset+=topdat[i].nmol*topdat[i].nat;
		}
	}

	/* if ther are angles....print them out */

	if(sysdat->total_angs>0){
		fprintf(fpout,"\n");
		fprintf(fpout,"Angles\n");
		fprintf(fpout,"\n");
		angleindex=0;
		offset=0;
		for(i=0;i<sysdat->ntops;i++) {
			for(j=0;j<topdat[i].nmol;j++){
				molindex++;
				for(k=0;k<topdat[i].nang;k++){
					angleindex++;
					fprintf(fpout,"%d %d %d %d %d # %s %s %s\n",angleindex,topdat[i].angtype[k]+1,
						topdat[i].angndx1[k]+(j*topdat[i].nat)+offset,topdat[i].angndx2[k]+(j*topdat[i].nat)+offset,
						topdat[i].angndx3[k]+(j*topdat[i].nat)+offset,topdat[i].type[topdat[i].angndx1[k]-1],
						topdat[i].type[topdat[i].angndx2[k]-1],topdat[i].type[topdat[i].angndx3[k]-1]);
				}
			}
			offset+=topdat[i].nmol*topdat[i].nat;
		}
	}

	/* if ther are dihedrals....print them out */

	if(sysdat->total_diheds>0){
		fprintf(fpout,"\n");
		fprintf(fpout,"Dihedrals\n");
		fprintf(fpout,"\n");
		dihedindex=0;
		offset=0;
		for(i=0;i<sysdat->ntops;i++) {
			for(j=0;j<topdat[i].nmol;j++){
				molindex++;
				for(k=0;k<topdat[i].ndihed;k++){
					dihedindex++;
					fprintf(fpout,"%d %d %d %d %d %d # %s %s %s %s\n",dihedindex,topdat[i].dihedtype[k]+1,
						topdat[i].dihedndx1[k]+(j*topdat[i].nat)+offset,topdat[i].dihedndx2[k]+(j*topdat[i].nat)+offset,
						topdat[i].dihedndx3[k]+(j*topdat[i].nat)+offset,topdat[i].dihedndx4[k]+(j*topdat[i].nat)+offset , topdat[i].type[topdat[i].dihedndx1[k]-1],
						topdat[i].type[topdat[i].dihedndx2[k]-1],topdat[i].type[topdat[i].dihedndx3[k]-1],topdat[i].type[topdat[i].dihedndx4[k]-1]);
				}
			}
			offset+=topdat[i].nmol*topdat[i].nat;
		}
	}

	/* if ther are impropers....print them out */

	if(sysdat->total_improps>0){
		fprintf(fpout,"\n");
		fprintf(fpout,"Impropers\n");
		fprintf(fpout,"\n");
		impropindex=0;
		offset=0;
		for(i=0;i<sysdat->ntops;i++) {
			for(j=0;j<topdat[i].nmol;j++){
				molindex++;
				for(k=0;k<topdat[i].nimprop;k++){
					impropindex++;
					fprintf(fpout,"%d %d %d %d %d %d # %s %s %s %s\n",impropindex,topdat[i].improptype[k]+1,
						topdat[i].impropndx1[k]+(j*topdat[i].nat)+offset,topdat[i].impropndx2[k]+(j*topdat[i].nat)+offset,
						topdat[i].impropndx3[k]+(j*topdat[i].nat)+offset,topdat[i].impropndx4[k]+(j*topdat[i].nat)+offset , topdat[i].type[topdat[i].impropndx1[k]-1],
						topdat[i].type[topdat[i].impropndx2[k]-1],topdat[i].type[topdat[i].impropndx3[k]-1],topdat[i].type[topdat[i].impropndx4[k]-1]);
				}
			}
			offset+=topdat[i].nmol*topdat[i].nat;
		}
	}


}
int min(const int a, const int b){
	if (a<b) return a;
	else     return b;
}
void write_psf(char *filename,DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat, int ISCHARGED){
	// this function is modified by read_coords()

	FILE *fpout;
	if((fpout = fopen("out.psf","w")) == NULL){fprintf(stderr,"ERROR: can't open out.psf\n");exit(1);}
	fprintf(fpout,"PSF \n");
	fprintf(fpout,"\n");
	fprintf(fpout,"       2 !NTITLE\n");
	fprintf(fpout,"* created by setup_lammps\n");
	fprintf(fpout,"* dummy\n");
	fprintf(fpout,"\n");
	fprintf(fpout,"%8d !NATOM\n",sysdat->total_ats);

	int atindex  = 0;
	int molindex = 0;
	for(int i=0;i<sysdat->ntops;i++) {
		for(int j=0;j<topdat[i].nmol;j++){
			molindex++;
			for(int k=0;k<topdat[i].nat;k++){
				atindex++;
//				      96 DMPC    8 DMPC C11  C      0.000000       42.0804
				fprintf(fpout,"%8d %-4s%5d %-4s %-4s %-4s  %9.6f  %12.4f\n",
					atindex, topdat[i].resname[k], min(9999,molindex), topdat[i].resname[k], topdat[i].atomname[k],topdat[i].type[k],
					topdat[i].charge[k], topdat[i].mass[k]);
			}
		}
	}
	fprintf(fpout,"\n");


	fprintf(fpout,"%8d !NBOND: bonds\n",sysdat->total_bnds);
	int bondindex=0;
	int offset   =0;
	for(int i=0;i<sysdat->ntops;i++) {
		for(int j=0;j<topdat[i].nmol;j++){
			for(int k=0;k<topdat[i].nbnd;k++){
				bondindex++;
				fprintf(fpout,"%8d%8d",topdat[i].bndndx1[k]+(j*topdat[i].nat)+offset,topdat[i].bndndx2[k]+(j*topdat[i].nat)+offset);
				if (bondindex%4==0) fprintf(fpout,"\n");
			}
		}
		offset+=topdat[i].nmol*topdat[i].nat;
	}
	fprintf(fpout,"\n");


	fprintf(fpout,"%-7d !NTHETA: angles\n",sysdat->total_angs);
	fprintf(fpout,"%-7d !NPHI: dihedrals\n",sysdat->total_diheds);
	fprintf(fpout,"%-7d !NIMPHI: impropers\n",sysdat->total_improps);
	//fprintf(fpout,"       0 !NTHETA: angles\n\n",sysdat->total_angs);
	//fprintf(fpout,"       0 !NPHI: dihedrals\n\n",sysdat->total_diheds);
	//fprintf(fpout,"       0 !NIMPHI: impropers\n\n",sysdat->total_improps);
	fprintf(fpout,"       0 !NDON: donors\n\n");
	fprintf(fpout,"       0 !NACC: acceptors\n\n");

	fclose(fpout);
}

// assuming three charactor string, "GB*".
int strcmp_wc(const char *s1, const char *s2){
	if (s1[2]=='*' || s2[2]=='*'){
		return memcmp(s1, s2, 2);
	}else{
		return strcmp(s1, s2);
	}
}

/* Decide which params will be needed */
void get_unique(DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat,int ISCHARGED){
	FILE *fpout;
	int i,j,k;
	char (*uniq_atype)[5];
	double *uniq_mass,*uniq_charge;
	int uniq_nats=0,uniq_bnds=0,uniq_angs=0,uniq_improps=0,uniq_diheds=0,ikeep,ifound,vdwtmp;
	int *bnd_params,*ang_params,*ang_vdw,datndx;

	if((fpout = fopen("PARM.FILE","w")) == NULL){
		fprintf(stderr,"ERROR: can't open PARM.FILE\n");
		exit(1);
	}

	uniq_atype=(char (*)[5])calloc(sysdat->nats,sizeof(*uniq_atype)); 
	uniq_mass=(double *)calloc(sysdat->nats,sizeof(double));
	uniq_charge=(double *)calloc(sysdat->nats,sizeof(double));
	bnd_params=(int *)calloc(MaxNBondPerAtom*sysdat->nats,sizeof(int));
	ang_params=(int *)calloc(MaxNBondPerAtom*sysdat->nats,sizeof(int));
	ang_vdw=(int *)calloc(MaxNBondPerAtom*sysdat->nats,sizeof(int));

	fprintf(fpout,"# Generated by setup_lammps\n");
	fprintf(fpout,"\n");

	if(ISCHARGED){
		fprintf(fpout,"pair_style      lj/sdk/coul/long        15.0\n");
		//fprintf(fpout,"kspace_style    pppm/cg 1.0e-5 1.0e-6\n");
        } else
		fprintf(fpout,"pair_style      cg/cmm        15.0\n");

	if(sysdat->nbnds > 0)
		fprintf(fpout,"bond_style      harmonic\n");

	if(sysdat->nangs > 0)
		fprintf(fpout,"angle_style     hybrid sdk harmonic\n");

	if(sysdat->total_diheds>0)
		fprintf(fpout,"dihedral_style  charmm\n");

	if(sysdat->total_improps>0)
		fprintf(fpout,"improper_style  harmonic\n");

	fprintf(fpout,"special_bonds   lj 0.0 0.0 1.0\n");

	fprintf(fpout,"\n");

	/* first gather the unique atom types */

	uniq_nats=0;
	for(i=0;i<sysdat->ntops;i++){
		for(j=0;j<topdat[i].nat;j++){
			ikeep=1;
			for(k=0;k<uniq_nats;k++){
				if(!strcmp(topdat[i].type[j],uniq_atype[k])){
					ikeep=0;
					topdat[i].parm_atomtype[j]=k;
					k=uniq_nats;
				}
			}
			if(ikeep==1) {
				strcpy(uniq_atype[uniq_nats],topdat[i].type[j]);
				uniq_mass[uniq_nats]=topdat[i].mass[j];
				uniq_charge[uniq_nats]=topdat[i].charge[j];
				topdat[i].parm_atomtype[j]=uniq_nats;
				uniq_nats++;
			}
		}
	}

	printf1(__LINE__,"FOUND %d UNIQUE ATOMS\n",uniq_nats);
	sysdat->uniq_nats=uniq_nats;
	for(i=0;i<uniq_nats;i++)
		fprintf(fpout,"mass   %-5d  %6.4f # %s\n",i+1, uniq_mass[i],uniq_atype[i]);

	fprintf(fpout,"\n");

	/* get pair interactions */


	for(i=0;i<uniq_nats;i++){
		for(j=i;j<uniq_nats;j++){
			ifound=0;
			for(k=0;k<database->nvdwtype;k++){
				if(!strcmp(database->vdwtype1[k],uniq_atype[i]) && !strcmp(database->vdwtype2[k],uniq_atype[j])) {
					/*           printf1(__LINE__,"Found database entry\n"); */
					ifound=1;
					vdwtmp=k;
					k=database->nvdwtype;
				}
				else if (!strcmp(database->vdwtype2[k],uniq_atype[i]) && !strcmp(database->vdwtype1[k],uniq_atype[j])) {
					/*           printf1(__LINE__,"Found database entry\n"); */
					ifound=1;
					vdwtmp=k;
					k=database->nvdwtype;
				}
			}
			if(ifound==0){
				printf1(__LINE__,"*********************\n");
				printf1(__LINE__,"WARNING:No params for VDW interaction between %s and %s\n",uniq_atype[i],uniq_atype[j]);
				printf1(__LINE__,"WILL CHECK FOR WILDCARDS!!!\n");
				/*         exit(1); */
			} else if (ifound==1){
				/*         printf1(__LINE__,"Found database entry for VDW\n"); */
				fprintf(fpout,"pair_coeff  %-5d %-5d %-6s %5.4f %5.4f # %-4s %-4s\n",i+1,j+1,database->vdwstyle[vdwtmp],
					database->eps[vdwtmp],database->sig[vdwtmp], database->vdwtype1[vdwtmp], database->vdwtype2[vdwtmp]);
			}
		}
	}

	/* BONDSSSSSSSSSS */

	/* get bond interactions  */
	fprintf(fpout,"\n");

	if(sysdat->nbnds > 0) {
		int index0 = 0;
		for(i=0;i<sysdat->ntops;i++){
			for(j=0;j<topdat[i].nbnd;j++){
				datndx=-1;

				/* AT THIS POINT WE WILL CHECK TO SEE IF THE PARAMS WERE GIVEN */
				/* IN THE TOP FILE....IF SO WE WILL SKIP A LOT OF THIS AND ADD */
				/* THIS AS A UNIQUE BOND....IF NOT WE GO THROUGH THE PROCEDURE */
				if ( topdat[i].bndpset[j] == false ) {
					/* now compare to the database */
					for(k=0;k<database->nbndtype;k++){
						if(!strcmp_wc(database->bndtype1[k], topdat[i].type[topdat[i].bndndx1[j]-1])
						&& !strcmp_wc(database->bndtype2[k], topdat[i].type[topdat[i].bndndx2[j]-1])){
								datndx=k;
								k=database->nbndtype;
						}
						if(!strcmp_wc(database->bndtype2[k], topdat[i].type[topdat[i].bndndx1[j]-1])
						&& !strcmp_wc(database->bndtype1[k], topdat[i].type[topdat[i].bndndx2[j]-1])){
								datndx=k;
								k=database->nbndtype;
						}
					}
/*
					// If we did not find database params specifically for this interaction, look for wildcards 
					// Wildcards are of the type GB*-LYS  
					if(datndx==-1){
						for(k=0;k<database->nbndtype;k++){

							// Check to see if database entry is a wildcard type...does it have a "*" 
							char *pt1=database->bndtype1[k];
							char *pt2=database->bndtype2[k];
							const int n1 = strlen(pt1);
							const int n2 = strlen(pt2);
							char *ptb1=topdat[i].type[topdat[i].bndndx1[j]-1];
							char *ptb2=topdat[i].type[topdat[i].bndndx2[j]-1];


							if (!strcmp(&pt1[n1-1], "*" ) || !strcmp(&pt2[n2-1], "*" )) {
								if(!memcmp(pt1, ptb1, n1-1) && !memcmp(pt2, ptb2, n2-1)) {
									printf1(__LINE__,"WARNING: USING WILDCARD FOR BOND!\n");
									printf1(__LINE__,"%s %s FOR %s %s\n", &pt1[0], &pt2[0], &ptb1[0], &ptb2[0]);
									printf1(__LINE__,"\n");
									datndx=k;
									k=database->nbndtype;
								} else if(!memcmp(pt1, ptb2, n1-1) && !memcmp(pt2, ptb1, n2-1)) {
									printf1(__LINE__,"WARNING: USING WILDCARD FOR BOND!\n");
									printf1(__LINE__,"%s %s FOR %s %s\n", &pt1[0], &pt2[0], &ptb1[0], &ptb2[0]);
									printf1(__LINE__,"\n");
									datndx=k;
									k=database->nbndtype;
								}
							} 
						}
					}	  
*/
					if(datndx==-1){
						printf1(__LINE__,"ERROR: DID NOT FIND BOND PARAMETERS IN DATABASE %d %d %s %s\n",
							topdat[i].bndndx1[j],
							topdat[i].bndndx2[j],
							topdat[i].type[topdat[i].bndndx1[j]-1],
							topdat[i].type[topdat[i].bndndx2[j]-1]);
						exit(1);
					}

					//// Check backbone
					//bool bBackboneBond;
					//if (((!memcmp(pt1,"GB",2) || !memcmp(pt1,"AB",2)))
					//||	((!memcmp(pt2,"GB",2) || !memcmp(pt2,"AB",2))) ) {
					//	bBackboneBond = true;
					//}else{
					//	bBackboneBond = false;
					//}

					///* Now make sure we do not already know we have this interaction */
					ikeep = 1;
					if (database->bnde[datndx] > 0){
						for(k=0;k<uniq_bnds;k++){
							if(datndx==bnd_params[k]){
								// found a replica
								ikeep=0;
								topdat[i].bndtype[j]=k;
								// kill the for loop
								k=uniq_bnds;
							}
						}
					}


					/* IKEEP = 1 IF WE FOUND A NEW ONE */
					if(ikeep==1) {
						bnd_params[uniq_bnds]=datndx;
						sysdat->param_bnds[uniq_bnds]=datndx;
						topdat[i].bndtype[j]=uniq_bnds;
						uniq_bnds++;
						printf1(__LINE__,"#####  %d\n", datndx);
						if (database->bnde[datndx] > 0){
							fprintf(fpout,"bond_coeff  %-6d %8.4f %8.4f # %s %s\n",uniq_bnds,database->fbnd[bnd_params[uniq_bnds-1]],
								database->bnde[bnd_params[uniq_bnds-1]],database->bndtype1[bnd_params[uniq_bnds-1]],database->bndtype2[bnd_params[uniq_bnds-1]]);
						}else{
							const int i1 = topdat[i].bndndx1[j]-1 + index0;
							const int i2 = topdat[i].bndndx2[j]-1 + index0;
							const double dx = sysdat->coordx[i1] - sysdat->coordx[i2];
							const double dy = sysdat->coordy[i1] - sysdat->coordy[i2];
							const double dz = sysdat->coordz[i1] - sysdat->coordz[i2];
							const double bond_in_pdb = sqrt(dx*dx+dy*dy+dz*dz);
							fprintf(fpout,"bond_coeff  %-6d %8.4f %8.4f # %s %s\n",uniq_bnds,database->fbnd[bnd_params[uniq_bnds-1]],
								bond_in_pdb,database->bndtype1[bnd_params[uniq_bnds-1]],database->bndtype2[bnd_params[uniq_bnds-1]]);
						}
					}
				} else {
					/* THE PARAMS WERE GIVEN IN THE TOP FILE SO LETS ADD IT TO THE PARAM FILE */
					topdat[i].bndtype[j]=uniq_bnds;
					bnd_params[uniq_bnds] = -1;
					uniq_bnds++;

					fprintf(fpout,"bond_coeff  %-6d %8.4f %8.4f # %s %s FROM TOP\n",uniq_bnds, topdat[i].bndfk[j],topdat[i].bndeq[j],
						topdat[i].type[topdat[i].bndndx1[j]-1], topdat[i].type[topdat[i].bndndx2[j]-1] );

				}
			}
			index0 += topdat[i].nat*topdat[i].nmol;
		}

		/* FINISHED LOOPING OVER TOP FILES */
		/* NOW DUMP THE BOND PARAMS */

		/*     for(i=0;i<uniq_bnds;i++){ */
		/*       fprintf(fpout,"bond_coeff  %-6d %8.4f %8.4f # %s %s\n",i+1,database->fbnd[bnd_params[i]], */
		/*          database->bnde[bnd_params[i]],database->bndtype1[bnd_params[i]],database->bndtype2[bnd_params[i]]); */
		/*     } */
	}
	sysdat->uniq_nbnds=uniq_bnds;


	/* ANGLESSSSSSSSSSS */

	/* get angle interactions if needed */
	fprintf(fpout,"\n");

	if(sysdat->nangs > 0) {

		uniq_angs=0;
		int index0 = 0;

		for(i=0;i<sysdat->ntops;i++){
			for(j=0;j<topdat[i].nang;j++){
				datndx=-1;

				/* AT THIS POINT WE WILL CHECK TO SEE IF THE PARAMS WERE GIVEN */
				/* IN THE TOP FILE....IF SO WE WILL SKIP A LOT OF THIS AND ADD */
				/* THIS AS A UNIQUE BOND....IF NOT WE GO THROUGH THE PROCEDURE */

				if ( topdat[i].angpset[j] == -1 ) {

					/* now compare to the database */

					for(k=0;k<database->nangtype;k++){
						if(!strcmp_wc(database->angtype2[k], topdat[i].type[topdat[i].angndx2[j]-1])){
							if(!strcmp_wc(database->angtype1[k], topdat[i].type[topdat[i].angndx1[j]-1])
								&& !strcmp_wc(database->angtype3[k], topdat[i].type[topdat[i].angndx3[j]-1])){
									datndx=k;
									k=database->nangtype;
							}
							if(!strcmp_wc(database->angtype3[k], topdat[i].type[topdat[i].angndx1[j]-1])
								&& !strcmp_wc(database->angtype1[k], topdat[i].type[topdat[i].angndx3[j]-1])){
									datndx=k;
									k=database->nangtype;
							}
						}
					}
/*
					if(datndx==-1){
						printf1(__LINE__,"WARNING: DID NOT FIND EXPLICIT ANGLE PARAMETERS IN DATABASE %s %s %s (%d %d %d )\n",topdat[i].type[topdat[i].angndx1[j]-1],
							topdat[i].type[topdat[i].angndx2[j]-1],topdat[i].type[topdat[i].angndx3[j]-1],topdat[i].angndx1[j],topdat[i].angndx2[j],topdat[i].angndx3[j] );
						printf1(__LINE__,"WILL CHECK FOR WILDCARDS!\n");

						// ######## WILD CARDS ############# 

						for(k=0;k<database->nangtype;k++){

							char *pt1=database->angtype1[k];
							char *pt2=database->angtype2[k];	      
							char *pt3=database->angtype3[k];
							int n1=strlen(pt1);
							int n2=strlen(pt2);
							int n3=strlen(pt3);
							char *ptb1=topdat[i].type[topdat[i].angndx1[j]-1];
							char *ptb2=topdat[i].type[topdat[i].angndx2[j]-1];
							char *ptb3=topdat[i].type[topdat[i].angndx3[j]-1];

							if (!strcmp(&pt1[n1-1], "*" ) || !strcmp(&pt2[n2-1], "*" ) || !strcmp(&pt3[n3-1], "*" )) {
								// 		printf1(__LINE__,"FOUND ANGLE WILDCARD: %s %s %s\n", &pt1[0], &pt2[0], &pt3[0] ); 

								if(!memcmp(pt2, ptb2, n2-1)) {
									// 		  printf1(__LINE__,"FOUND ANGLE WILDCARD: %s %s %s\n", &pt1[0], &pt2[0], &pt3[0] ); 
									if(!memcmp(pt1, ptb1, n1-1) && !memcmp(pt3, ptb3, n3-1) ) {
										datndx=k;
										k=database->nangtype;
									}  else if(!memcmp(pt1, ptb3, n1-1)  && !memcmp(pt3, ptb1, n3-1) ) {
										datndx=k;
										k=database->nangtype;
									}
								}
							}
						}
					}
*/
					/* RHD Get the VDW for the CG angles */
					ifound=0;
					for(k=0;k<database->nvdwtype;k++){
						if(!strcmp(database->vdwtype1[k],topdat[i].type[topdat[i].angndx1[j]-1]) && !strcmp(database->vdwtype2[k],
							topdat[i].type[topdat[i].angndx3[j]-1])) {
								ifound=1;
								vdwtmp=k;
								k=database->nvdwtype;
						}
						else if (!strcmp(database->vdwtype1[k],topdat[i].type[topdat[i].angndx3[j]-1]) && !strcmp(database->vdwtype2[k],
							topdat[i].type[topdat[i].angndx1[j]-1])) {
								ifound=1;
								vdwtmp=k;
								k=database->nvdwtype;
						}
					}

					if(ifound==0){
						printf1(__LINE__,"*********************\n");
						printf1(__LINE__,"ERROR:No params for VDW interaction between %s and %s for angle (database)\n",topdat[i].type[topdat[i].angndx1[j]-1],
							topdat[i].type[topdat[i].angndx3[j]-1]);
						printf1(__LINE__,"UPDATE DATABASE!!!\n");
						/*             exit(1); */
					}

					/* end VDW for CG angles */

					/* No params for this interaction in the database */

					if(datndx==-1){
						printf1(__LINE__,"ERROR: DID NOT FIND ANGLE PARAMETERS IN DATABASE %s %s %s (%d %d %d )\n",topdat[i].type[topdat[i].angndx1[j]-1],
							topdat[i].type[topdat[i].angndx2[j]-1],topdat[i].type[topdat[i].angndx3[j]-1],topdat[i].angndx1[j],topdat[i].angndx2[j],topdat[i].angndx3[j] );
						exit(1);
					}

					/* ############################ */



					/* Now make sure we do not already have this one */
					ikeep = 1;
					if (database->ange[datndx] > 0){
						for(k=0;k<uniq_angs;k++){
							if(datndx==ang_params[k]){
								ikeep=0;
								topdat[i].angtype[j]=k;
								k=uniq_angs;
							}
						}
					}

					if(ikeep==1) {
						ang_params[uniq_angs]=datndx;
						ang_vdw[uniq_angs]=vdwtmp;
						sysdat->param_angs[uniq_angs]=datndx;
						topdat[i].angtype[j]=uniq_angs;
						uniq_angs++;
						if (database->ange[datndx] > 0){
							if(database->angsdk[ang_params[ uniq_angs-1  ]]){
								fprintf(fpout,"angle_coeff %-10d      sdk  %8.4f %8.4f %s %8.4f %8.4f # %s %s %s\n",uniq_angs,database->fang[ang_params[uniq_angs-1]],
									database->ange[ang_params[ uniq_angs-1   ]],database->vdwstyle[ang_vdw[ uniq_angs-1  ]],
									database->eps[ang_vdw[ uniq_angs-1  ]],database->sig[ang_vdw[ uniq_angs-1  ]],
									database->angtype1[ang_params[ uniq_angs-1  ]],database->angtype2[ang_params[ uniq_angs-1  ]],
									database->angtype3[ang_params[ uniq_angs-1  ]]);
							}else{
								fprintf(fpout,"angle_coeff %-10d harmonic  %8.4f %8.4f # %s %s %s\n",
									uniq_angs,database->fang[ang_params[uniq_angs-1]],database->ange[ang_params[ uniq_angs-1   ]],
									database->angtype1[ang_params[ uniq_angs-1  ]],database->angtype2[ang_params[ uniq_angs-1  ]],
									database->angtype3[ang_params[ uniq_angs-1  ]]);
							}
						}else{
							const int i1 = topdat[i].angndx1[j]-1 + index0;
							const int i2 = topdat[i].angndx2[j]-1 + index0;
							const int i3 = topdat[i].angndx3[j]-1 + index0;
							const double r1[3] = {sysdat->coordx[i1],sysdat->coordy[i1],sysdat->coordz[i1]};
							const double r2[3] = {sysdat->coordx[i2],sysdat->coordy[i2],sysdat->coordz[i2]};
							const double r3[3] = {sysdat->coordx[i3],sysdat->coordy[i3],sysdat->coordz[i3]};
							const double angle_in_pdb = 180.0/PI*GetAngle(r1,r2,r3);
							if(database->angsdk[ang_params[ uniq_angs-1  ]]){
								fprintf(fpout,"angle_coeff %-10d      sdk  %8.4f %8.4f %s %8.4f %8.4f # %s %s %s\n",
									uniq_angs,
									database->fang[ang_params[uniq_angs-1]],
									angle_in_pdb,
									database->vdwstyle[ang_vdw[ uniq_angs-1]],
									database->eps[ang_vdw[uniq_angs-1]],
									database->sig[ang_vdw[uniq_angs-1]],
									database->angtype1[ang_params[uniq_angs-1]],
									database->angtype2[ang_params[uniq_angs-1]],
									database->angtype3[ang_params[uniq_angs-1]]);
							}else{
								fprintf(fpout,"angle_coeff %-10d harmonic  %8.4f %8.4f # %s %s %s\n",
									uniq_angs,
									database->fang[ang_params[uniq_angs-1]],
									angle_in_pdb,
									database->angtype1[ang_params[uniq_angs-1]],
									database->angtype2[ang_params[uniq_angs-1]],
									database->angtype3[ang_params[uniq_angs-1]]);
							}
						}
					}
				} else {
					/* THIS PARAM WAS SET IN THE TOP FILE */

					/* still need vdw stuff */
					/* RHD Get the VDW for the CG angles */
					ifound=0;
					for(k=0;k<database->nvdwtype;k++){
						if(!strcmp(database->vdwtype1[k],topdat[i].type[topdat[i].angndx1[j]-1]) && !strcmp(database->vdwtype2[k],
							topdat[i].type[topdat[i].angndx3[j]-1])) {
								ifound=1;
								vdwtmp=k;
								k=database->nvdwtype;
						}
						else if (!strcmp(database->vdwtype1[k],topdat[i].type[topdat[i].angndx3[j]-1]) && !strcmp(database->vdwtype2[k],
							topdat[i].type[topdat[i].angndx1[j]-1])) {
								ifound=1;
								vdwtmp=k;
								k=database->nvdwtype;
						}
					}
					if(ifound==0){
						printf1(__LINE__,"*********************\n");
						printf1(__LINE__,"ERROR:No params for VDW interaction between %s and %s for angle (topfile)\n",topdat[i].type[topdat[i].angndx1[j]-1],
							topdat[i].type[topdat[i].angndx3[j]-1]);
						printf1(__LINE__,"UPDATE DATABASE!!!\n");
						exit(1);
					}

					/* end VDW for CG angles */

					topdat[i].angtype[j]=uniq_angs;
					uniq_angs++;
					if(database->angsdk[ang_params[ uniq_angs-1  ]]){
						fprintf(fpout,"angle_coeff %-10d      sdk  %8.4f %8.4f %s %8.4f %8.4f # %s %s %s FROM TOP\n",
						uniq_angs,topdat[i].angfk[j],topdat[i].angeq[j],
						database->vdwstyle[vdwtmp],database->eps[vdwtmp],database->sig[vdwtmp],
						topdat[i].type[topdat[i].angndx1[j]-1],topdat[i].type[topdat[i].angndx2[j]-1] ,topdat[i].type[topdat[i].angndx3[j]-1]);
					}else{
						fprintf(fpout,"angle_coeff %-10d harmonic  %8.4f %8.4f # %s %s %s FROM TOP\n",
						uniq_angs,topdat[i].angfk[j],topdat[i].angeq[j],
						topdat[i].type[topdat[i].angndx1[j]-1],topdat[i].type[topdat[i].angndx2[j]-1] ,topdat[i].type[topdat[i].angndx3[j]-1]);
					}
				}
			}
			index0 += topdat[i].nat*topdat[i].nmol;
		}
	}
	sysdat->uniq_nangs=uniq_angs;

	/* NOW LET HANDLE THE DIHEDRAL PARAMS THE ONLY WAY WE DO...THEY HAVE TO BE SPECIFIED IN THE TOP FILE */
	fprintf(fpout,"\n");

	printf1(__LINE__,"DIHEDS TEST %d\n",sysdat->total_diheds);
	if(sysdat->total_diheds > 0) {
		uniq_diheds=0;
		int index0 = 0;

		for(i=0;i<sysdat->ntops;i++){
			for(j=0;j<topdat[i].ndihed;j++){

				if ( topdat[i].dihedpset[j] == 1 ) {
					topdat[i].dihedtype[j]=uniq_diheds;         
					uniq_diheds++;
					///* RHDDDD */
					//fprintf(fpout,"dihedral_coeff %-10d %8.4f %-3d %-3d %2.1f # %s %s %s %s FROM TOP\n",uniq_diheds,topdat[i].dihedfk[j],topdat[i].dihedn[j],
					//	topdat[i].dihedeq[j],topdat[i].dihedof[j],topdat[i].type[topdat[i].dihedndx1[j]-1],topdat[i].type[topdat[i].dihedndx2[j]-1],
					//	topdat[i].type[topdat[i].dihedndx3[j]-1],topdat[i].type[topdat[i].dihedndx4[j]-1]);
					const int i1 = topdat[i].dihedndx1[j]-1 + index0;
					const int i2 = topdat[i].dihedndx2[j]-1 + index0;
					const int i3 = topdat[i].dihedndx3[j]-1 + index0;
					const int i4 = topdat[i].dihedndx4[j]-1 + index0;
					const double r1[3] = {sysdat->coordx[i1],sysdat->coordy[i1],sysdat->coordz[i1]};
					const double r2[3] = {sysdat->coordx[i2],sysdat->coordy[i2],sysdat->coordz[i2]};
					const double r3[3] = {sysdat->coordx[i3],sysdat->coordy[i3],sysdat->coordz[i3]};
					const double r4[3] = {sysdat->coordx[i4],sysdat->coordy[i4],sysdat->coordz[i4]};
					const double dihedral_in_pdb = 180.0 + 180.0/PI*GetDihedral(r1,r2,r3,r4);
					int idihedral_in_pdb;
					if (dihedral_in_pdb>0.0) idihedral_in_pdb = int(dihedral_in_pdb+0.5);
					else                     idihedral_in_pdb = int(dihedral_in_pdb-0.5);
					fprintf(fpout,
					"dihedral_coeff %-10d %8.4f %-3d %-3d %2.1f # %s %s %s %s\n",
					uniq_diheds,
					topdat[i].dihedfk[j],
					topdat[i].dihedn[j],
//					idihedral_in_pdb,
					topdat[i].dihedeq[j],
					topdat[i].dihedof[j],
					topdat[i].type[topdat[i].dihedndx1[j]-1],
					topdat[i].type[topdat[i].dihedndx2[j]-1],
					topdat[i].type[topdat[i].dihedndx3[j]-1],
					topdat[i].type[topdat[i].dihedndx4[j]-1]);

				}
			}
			index0 += topdat[i].nat*topdat[i].nmol;
		}
	}
	sysdat->uniq_ndiheds=uniq_diheds;

	/* NOW LET HANDLE THE IMPROPS PARAMS THE ONLY WAY WE DO...THEY HAVE TO BE SPECIFIED IN THE TOP FILE */
	fprintf(fpout,"\n");

	printf1(__LINE__,"IMPROPS TEST %d\n",sysdat->total_improps);
	if(sysdat->total_improps > 0) {
		uniq_improps=0;

		for(i=0;i<sysdat->ntops;i++){
			for(j=0;j<topdat[i].nimprop;j++){

				if ( topdat[i].improppset[j] == 1 ) {
					topdat[i].improptype[j]=uniq_improps;         
					uniq_improps++;
					/* RHDDDD */
					fprintf(fpout,"improper_coeff %-10d %8.4f %8.4f # %s %s %s %s FROM TOP\n",uniq_improps,topdat[i].impropfk[j],topdat[i].impropeq[j],
						topdat[i].type[topdat[i].impropndx1[j]-1],topdat[i].type[topdat[i].impropndx2[j]-1] ,topdat[i].type[topdat[i].impropndx3[j]-1],
						topdat[i].type[topdat[i].impropndx4[j]-1]);

				}
			}
		}
	}
	sysdat->uniq_nimprops=uniq_improps;

}

/* Read the database file and store unique params */
/* Warn if you find duplicates */
void read_database(char *filename, DATABASE *database) {
	FILE *fpin;
	int nvdw,nbnd,nang,i,ikeep;
	char col1[100],vdwtype1[5],vdwtype2[5],vdwstyle[7];
	char bndtype1[5],bndtype2[5],angtype1[5],angtype2[5],angtype3[5];
	double eps,sig,fbnd,bnde,fang,ange, angsdk;
	const int Nbuf=500;
	char buf[Nbuf];
	char *line;

	if((fpin = fopen(filename,"r")) == NULL){
		fprintf(stderr,"ERROR: can't open infile %s\n",filename);
		exit(1);
	}
	nvdw=0;
	nbnd=0;
	nang=0;

	while(fgets(buf,Nbuf,fpin)){
		strcpy(col1,""); sscanf(buf,"%s", col1);
		if(!strcmp(col1, "pair")){

			line = &buf[strlen(col1)];
			sscanf(line,"%s %s %s %lf %lf",vdwtype1,vdwtype2,vdwstyle,&eps,&sig);

			ikeep=1;
			for(i=0;i<nvdw;i++){
				if(!strcmp(vdwtype1,database->vdwtype1[i]) && !strcmp(vdwtype2,database->vdwtype2[i])){
					printf1(__LINE__,"WARNING: FOUND DUP VDW PARAM %s %s\n",vdwtype1,vdwtype2);
					ikeep=0;
				}
				else if (!strcmp(vdwtype1,database->vdwtype2[i]) && !strcmp(vdwtype2,database->vdwtype1[i])){
					printf1(__LINE__,"WARNING: FOUND DUP VDW PARAM %s %s\n",vdwtype1,vdwtype2);
					ikeep=0;
				}
			}
			if(ikeep==1){
				strcpy(database->vdwtype1[nvdw],vdwtype1);
				strcpy(database->vdwtype2[nvdw],vdwtype2);
				strcpy(database->vdwstyle[nvdw],vdwstyle);
				database->eps[nvdw]=eps;
				database->sig[nvdw]=sig;
				nvdw++;
			}
		}
		if(!strcmp(col1, "bond")){
			line = &buf[strlen(col1)];
			sscanf(line,"%s %s %lf %lf",bndtype1,bndtype2,&fbnd,&bnde);

			ikeep=1;
			for(i=0;i<nbnd;i++){

				if(!strcmp(bndtype1,database->bndtype1[i]) && !strcmp(bndtype2,database->bndtype2[i])){
					printf1(__LINE__,"WARNING: FOUND DUP BOND PARAM %s %s\n",bndtype1,bndtype2);
					ikeep=0;
				}
				else if (!strcmp(bndtype1,database->bndtype2[i]) && !strcmp(bndtype2,database->bndtype1[i])){
					printf1(__LINE__,"WARNING: FOUND DUP BOND PARAM %s %s\n",bndtype1,bndtype2);
					ikeep=0;
				}
			}
			if(ikeep==1){
				strcpy(database->bndtype1[nbnd],bndtype1);
				strcpy(database->bndtype2[nbnd],bndtype2);
				database->fbnd[nbnd]=fbnd;
				database->bnde[nbnd]=bnde;
				nbnd++;
			}
		}

		if(!strcmp(col1, "angle")){
			
			line = &buf[strlen(col1)];
			if(strstr(line, "harmonic") != NULL) {
				angsdk = 0;
			}
			else {
				angsdk = 1;
			}
			sscanf(line,"%s %s %s %lf %lf",angtype1,angtype2,angtype3,&fang,&ange);

			ikeep=1;
			for(i=0;i<nang;i++){
				if(!strcmp(angtype2,database->angtype2[i])) {
					if(!strcmp(angtype1,database->angtype1[i]) && !strcmp(angtype3,database->angtype3[i])){
						printf1(__LINE__,"WARNING: FOUND DUP ANGLE PARAM %s %s %s\n",angtype1,angtype2,angtype3);
						ikeep=0;
					}
					else if (!strcmp(angtype3,database->angtype1[i]) && !strcmp(angtype1,database->angtype3[i])){
						printf1(__LINE__,"WARNING: FOUND DUP ANGLE PARAM %s %s %s\n",angtype1,angtype2,angtype3);
						ikeep=0;
					}
				}
			}
			if(ikeep==1){
				strcpy(database->angtype1[nang],angtype1);
				strcpy(database->angtype2[nang],angtype2);
				strcpy(database->angtype3[nang],angtype3);
				database->fang[nang]=fang;
				database->ange[nang]=ange;
				database->angsdk[nang]=angsdk;
				nang++;
			}
		}
	}
	fclose(fpin);

	database->nvdwtype=nvdw;
	database->nbndtype=nbnd;
	database->nangtype=nang;

}

/* Count the number of params in the database so we can allocate for storage */
void count_params(char *filename, DATABASE *database){
	FILE *fpin;
	const int Nbuf=500;
	char buf[Nbuf];
	char col1[Nbuf];

	if((fpin = fopen(filename,"r")) == NULL){
		fprintf(stderr,"ERROR: can't open infile %s\n",filename);
		exit(1);
	}
	database->nvdwtype=0;
	database->nbndtype=0;
	database->nangtype=0;

	while(fgets(buf,Nbuf,fpin)){
		strcpy(col1,""); sscanf(buf,"%s", col1);
		if(!strcmp(col1, "pair")){
			database->nvdwtype++;
		}
		if(!strcmp(col1, "bond")){
			database->nbndtype++;
		}
		if(!strcmp(col1, "angle")){
			database->nangtype++;
		}
	}
	fclose(fpin);
}

/* Read the topology file and store the data */
void read_top(char *filename, TOPDAT *topdat, int ntop, int *ISCHARGED){
	FILE *fpin;
	int ndx,bndx,andx,indx,dndx;
	const int Nbuf=500;
	char buf[Nbuf];
	char col1[Nbuf];
	char *line;

	if((fpin = fopen(filename,"r")) == NULL){
		fprintf(stderr,"ERROR: can't open infile %s\n",filename);
		exit(1);
	}
	ndx=0;
	bndx=0;
	andx=0;
	indx=0;
	dndx=0;

	printf1(__LINE__,"################\n");
	printf1(__LINE__,"##### READING %s\n", filename);

	while(fgets(buf,Nbuf,fpin)){
		strcpy(col1,""); sscanf(buf,"%s", col1);
		line = &buf[strlen(col1)];

		if(!strcmp(col1, "atom")){
			//char atomname[7];
			if(sscanf(line,"%d %s %s %s %lf %lf %s",
				&topdat[ntop].index[ndx],topdat[ntop].resname[ndx], topdat[ntop].atomname[ndx], 
				topdat[ntop].type[ndx], &topdat[ntop].mass[ndx],
				&topdat[ntop].charge[ndx], topdat[ntop].segid[ndx])==EOF){
				printf1(__LINE__,"Error at FILE %s, line %d\n", __FILE__,__LINE__);
			}
			if((topdat[ntop].charge[ndx]*topdat[ntop].charge[ndx]) > 0.00001){
				printf1(__LINE__,"CHARGE IN TOP FILE: %s %lf\n", filename, topdat[ntop].charge[ndx]);
				*ISCHARGED=1;
			}

			ndx++;
		}

		if(!strcmp(col1, "bond")){
			if (sscanf(line,"%d %d", &topdat[ntop].bndndx1[bndx], &topdat[ntop].bndndx2[bndx])==EOF){
				printf1(__LINE__,"Error at FILE %s, line %d\n", __FILE__,__LINE__);
			}
			topdat[ntop].bndpset[bndx] = false;
			bndx++;
		}

		if(!strcmp(col1, "bondparam")){
			printf1(__LINE__,"WARNING!!!!!!!: USING BOND PARAMETERS FROM THE TOP FILE\n");
			if(sscanf(line,"%d %d %lf %lf", &topdat[ntop].bndndx1[bndx], &topdat[ntop].bndndx2[bndx], &topdat[ntop].bndfk[bndx],
				&topdat[ntop].bndeq[bndx]) != 4) {

					printf1(__LINE__,"ERROR: NOT ENOUGH ARGS FOR BONDPARAM: MUST BE: NDX1 NDX2 FK EQ\n");
					exit(1);
			}

			topdat[ntop].bndpset[bndx] = true;
			bndx++;
		}

		if(!strcmp(col1, "angle")){
			if (sscanf(line,"%d %d %d", &topdat[ntop].angndx1[andx], &topdat[ntop].angndx2[andx], &topdat[ntop].angndx3[andx])==EOF){
				printf1(__LINE__,"Error at FILE %s, line %d\n", __FILE__,__LINE__);
			}
			topdat[ntop].angpset[andx] = -1;
			andx++;
		}

		if(!strcmp(col1, "angleparam")){
			printf1(__LINE__,"WARNING!!!!!!!: USING ANGLE PARAMETERS FROM THE TOP FILE\n");
			if(sscanf(line,"%d %d %d %lf %lf", &topdat[ntop].angndx1[andx], &topdat[ntop].angndx2[andx], &topdat[ntop].angndx3[andx],
				&topdat[ntop].angfk[andx],&topdat[ntop].angeq[andx]) != 5) {
					printf1(__LINE__,"ERROR: NOT ENOUGH ARGS FOR ANGLEPARAM: MUST BE: NDX1 NDX2 NDX3 FK EQ\n");
					exit(1);
			}

			topdat[ntop].angpset[andx] = 1;
			andx++;
		}

		if(!strcmp(col1, "improper")){
			printf1(__LINE__,"WARNING!!!!!!!: THIS IS NOT IMPLEMENTED.  MUST USE improperparam AND ASSIGN IMPROPER PARAMETERS IN THE TOP FILE\n");
			if (sscanf(line,"%d %d %d %d", &topdat[ntop].impropndx1[indx], &topdat[ntop].impropndx2[indx], &topdat[ntop].impropndx3[indx],&topdat[ntop].impropndx4[indx])==EOF){
				printf1(__LINE__,"Error at FILE %s, line %d\n", __FILE__,__LINE__);
			}
			topdat[ntop].improppset[indx] = -1;
			indx++;
		}

		if(!strcmp(col1, "improperparam")){
			printf1(__LINE__,"WARNING!!!!!!!: USING IMPROPER PARAMETERS FROM THE TOP FILE\n");
			if(sscanf(line,"%d %d %d %d %lf %lf", &topdat[ntop].impropndx1[indx], &topdat[ntop].impropndx2[indx],
				&topdat[ntop].impropndx3[indx],&topdat[ntop].impropndx4[indx],&topdat[ntop].impropfk[indx],&topdat[ntop].impropeq[indx]) != 6) {
					printf1(__LINE__,"ERROR: NOT ENOUGH ARGS FOR IMPROPERPARAM: MUST BE: NDX1 NDX2 NDX3 NDX4 FK EQ\n");
					exit(1);
			}
			topdat[ntop].improppset[indx] = 1;
			indx++;
		}

		if(!strcmp(col1, "dihedralparam")){
			printf1(__LINE__,"WARNING!!!!!!!: USING DIHEDRAL PARAMETERS FROM THE TOP FILE\n");
			if(sscanf(line,"%d %d %d %d %lf %d %d %lf", &topdat[ntop].dihedndx1[dndx], &topdat[ntop].dihedndx2[dndx],
				&topdat[ntop].dihedndx3[dndx],&topdat[ntop].dihedndx4[dndx],&topdat[ntop].dihedfk[dndx],&topdat[ntop].dihedn[dndx], &topdat[ntop].dihedeq[dndx], &topdat[ntop].dihedof[dndx]) != 8) {
					printf1(__LINE__,"ERROR: NOT ENOUGH ARGS FOR CHARMM DIHEDRALPARAM: MUST BE: NDX1 NDX2 NDX3 NDX4 FK n EQ ONEFOUR\n");
					exit(1);
			}
			topdat[ntop].dihedpset[dndx] = 1;
			dndx++;
		}


	}
	fclose(fpin);
}

/* count the number of things in the topology files so we can allocate */
void count_atoms(char *filename, TOPDAT *topdat, int ntop){
	FILE *fpin;
	const int Nbuf=500;
	char buf[Nbuf];
	char col1[Nbuf];

	if((fpin = fopen(filename,"r")) == NULL){
		fprintf(stderr,"ERROR: can't open infile %s\n",filename);
		exit(1);
	}
	topdat[ntop].nat=0;
	topdat[ntop].nbnd=0;
	topdat[ntop].nang=0;
	topdat[ntop].nimprop=0;
	topdat[ntop].ndihed=0;

	while(fgets(buf,Nbuf,fpin)){
		strcpy(col1,""); sscanf(buf,"%s", col1);
		if(!strcmp(col1, "atom")){
			topdat[ntop].nat++;
		}
		if(!strcmp(col1, "bond") || !strcmp(col1, "bondparam")){
			topdat[ntop].nbnd++;
		}
		if(!strcmp(col1, "angle") || !strcmp(col1, "angleparam")){
			topdat[ntop].nang++;
		}
		if(!strcmp(col1, "improperparam") || !strcmp(col1, "improper")){
			topdat[ntop].nimprop++;
		}
		if(!strcmp(col1, "dihedralparam") || !strcmp(col1, "dihedral")){
			topdat[ntop].ndihed++;
		}
	}
	
	if (topdat[ntop].nat==0){
		printf1(__LINE__,"ERROR: natom in %s is zero.\n", filename);
		exit(0);
	}
	
	fclose(fpin);
}


/* Main routine. Call and allocate                                       */
/* The idea is to read in the topologies and then check the database for */
/* all of the required interaction params.                               */


int main(int argc, char **argv ){
	int ntops,rdtp,i,ISCHARGED;
	TOPDAT topdat[1000];
	DATABASE database;
	SYSDAT sysdat;
	memset(&sysdat,0,sizeof(SYSDAT));

	if(argc < 5){
		printf1(__LINE__,"USAGE: setup_lammps <topfile 1> <nmol 1> [ <topfile 2> <nmol 2> ..... <topfile n> <nmol n>] <param file> <coordfile>\n");
		printf1(__LINE__,"Prints out input files for a lammps run.  Takes a pdb file as the coordfile\n");
		exit(1);
	}

	ISCHARGED=0;

	ntops=(argc-3)/2;
	printf("WILL READ %d TOPOLOGY FILE(S).\n\n",ntops);
	rdtp=0;

	/* Loop through the topologies and count the number of atoms, bonds and bends */

	while(rdtp < ntops){
		topdat[rdtp].nmol=atoi(argv[(2*rdtp)+2]);
		count_atoms(argv[(2*rdtp)+1],topdat,rdtp);
		rdtp++;
	}
	sysdat.ntops=ntops;
	sysdat.nats=0;
	sysdat.nbnds=0;
	sysdat.nangs=0;
	sysdat.total_ats=0;
	sysdat.total_bnds=0;
	sysdat.total_angs=0;
	sysdat.total_improps=0;
	sysdat.total_diheds=0;

	for(i=0;i<ntops;i++){
		printf("BOOKKEEPING: %d/%d\n", i,ntops);
		printf("TOPFILE %s\n",argv[(2*i)+1]);
		printf("FOUND: %d atoms\n",topdat[i].nat);
		printf("FOUND: %d bonds\n",topdat[i].nbnd);
		printf("FOUND: %d angles\n",topdat[i].nang);
		printf("FOUND: %d improperss\n",topdat[i].nimprop);
		printf("FOUND: %d dihedrals\n",topdat[i].ndihed);
		printf("\n");

		topdat[i].index   = (int *)calloc(topdat[i].nat,sizeof(int));
		topdat[i].mass    = (double *)calloc(topdat[i].nat,sizeof(double));
		topdat[i].charge  = (double *)calloc(topdat[i].nat,sizeof(double));

		topdat[i].bndndx1 = (int *)calloc(topdat[i].nbnd,sizeof(int));
		topdat[i].bndndx2 = (int *)calloc(topdat[i].nbnd,sizeof(int));
		topdat[i].bndtype = (int *)calloc(topdat[i].nbnd,sizeof(int));
		topdat[i].bndpset = (int *)calloc(topdat[i].nbnd,sizeof(int));
		topdat[i].bndfk = (double *)calloc(topdat[i].nbnd,sizeof(double));
		topdat[i].bndeq = (double *)calloc(topdat[i].nbnd,sizeof(double));

		topdat[i].angndx1 = (int *)calloc(topdat[i].nang,sizeof(int));
		topdat[i].angndx2 = (int *)calloc(topdat[i].nang,sizeof(int));
		topdat[i].angndx3 = (int *)calloc(topdat[i].nang,sizeof(int));
		topdat[i].angtype = (int *)calloc(topdat[i].nang,sizeof(int));
		topdat[i].angpset = (int *)calloc(topdat[i].nang,sizeof(int));
		topdat[i].angfk = (double *)calloc(topdat[i].nang,sizeof(double));
		topdat[i].angeq = (double *)calloc(topdat[i].nang,sizeof(double));

		topdat[i].impropndx1 = (int *)calloc(topdat[i].nimprop,sizeof(int));
		topdat[i].impropndx2 = (int *)calloc(topdat[i].nimprop,sizeof(int));
		topdat[i].impropndx3 = (int *)calloc(topdat[i].nimprop,sizeof(int));
		topdat[i].impropndx4 = (int *)calloc(topdat[i].nimprop,sizeof(int));
		topdat[i].improptype = (int *)calloc(topdat[i].nimprop,sizeof(int));
		topdat[i].improppset = (int *)calloc(topdat[i].nimprop,sizeof(int));
		topdat[i].impropfk = (double *)calloc(topdat[i].nimprop,sizeof(double));
		topdat[i].impropeq = (double *)calloc(topdat[i].nimprop,sizeof(double));

		topdat[i].dihedndx1 = (int *)calloc(topdat[i].ndihed,sizeof(int));
		topdat[i].dihedndx2 = (int *)calloc(topdat[i].ndihed,sizeof(int));
		topdat[i].dihedndx3 = (int *)calloc(topdat[i].ndihed,sizeof(int));
		topdat[i].dihedndx4 = (int *)calloc(topdat[i].ndihed,sizeof(int));
		topdat[i].dihedtype = (int *)calloc(topdat[i].ndihed,sizeof(int));
		topdat[i].dihedpset = (int *)calloc(topdat[i].ndihed,sizeof(int));
		topdat[i].dihedfk = (double *)calloc(topdat[i].ndihed,sizeof(double));
		topdat[i].dihedeq = (int *)calloc(topdat[i].ndihed,sizeof(int));
		topdat[i].dihedof = (double *)calloc(topdat[i].ndihed,sizeof(double));
		topdat[i].dihedn = (int *)calloc(topdat[i].ndihed,sizeof(int));

		topdat[i].parm_atomtype = (int *)calloc(topdat[i].nat,sizeof(int));
		topdat[i].atomname= (char (*)[5])calloc(topdat[i].nat,sizeof(*topdat[i].atomname));
		topdat[i].type    = (char (*)[5])calloc(topdat[i].nat,sizeof(*topdat[i].type));
		topdat[i].resname = (char (*)[5])calloc(topdat[i].nat,sizeof(*topdat[i].resname));
		topdat[i].segid   = (char (*)[5])calloc(topdat[i].nat,sizeof(*topdat[i].segid));

		sysdat.nats     += topdat[i].nat;
		sysdat.nbnds    += topdat[i].nbnd;
		sysdat.nangs    += topdat[i].nang;
		sysdat.nimprops += topdat[i].nimprop;
		sysdat.ndiheds  += topdat[i].ndihed;

		sysdat.total_ats     += topdat[i].nat*topdat[i].nmol;
		sysdat.total_bnds    += topdat[i].nbnd*topdat[i].nmol;
		sysdat.total_angs    += topdat[i].nang*topdat[i].nmol;
		sysdat.total_improps += topdat[i].nimprop*topdat[i].nmol;
		sysdat.total_diheds  += topdat[i].ndihed*topdat[i].nmol;
	}
	/*     printf1(__LINE__,"BOOKKEEPING:\n"); */
	/*     printf1(__LINE__,"TOTAL NONUNIQUE\n",i); */
	/*     printf1(__LINE__,"FOUND: %d atoms\n",sysdat.nats); */
	/*     printf1(__LINE__,"FOUND: %d bonds\n",sysdat.nbnds); */

	printf("TOTALS:\n");
	printf("FOUND: %d impropers\n",sysdat.total_improps);
	printf("FOUND: %d dihderalss\n",sysdat.total_diheds);


	rdtp=0;  
	while(rdtp < ntops){
		topdat[rdtp].nmol=atoi(argv[(2*rdtp+2)]);
		read_top(argv[(2*rdtp)+1],topdat,rdtp,&ISCHARGED);
		rdtp++;
	} 

	count_params(argv[argc-2],&database);

	database.fbnd=(double *)calloc(database.nbndtype,sizeof(double)); // changed by shuhei
	database.bnde=(double *)calloc(database.nbndtype,sizeof(double));
	database.fang=(double *)calloc(database.nangtype,sizeof(double));
	database.ange=(double *)calloc(database.nangtype,sizeof(double));
	database.angsdk=(double *)calloc(database.nangtype,sizeof(double));
	database.eps=(double *)calloc(database.nvdwtype,sizeof(double));
	database.sig=(double *)calloc(database.nvdwtype,sizeof(double));
	database.vdwtype1=(char (*)[5])calloc(database.nvdwtype,sizeof(*database.vdwtype1));
	database.vdwtype2=(char (*)[5])calloc(database.nvdwtype,sizeof(*database.vdwtype2));
	database.vdwstyle=(char (*)[7])calloc(database.nvdwtype,sizeof(*database.vdwstyle));
	database.bndtype1=(char (*)[5])calloc(database.nbndtype,sizeof(*database.bndtype1));
	database.bndtype2=(char (*)[5])calloc(database.nbndtype,sizeof(*database.bndtype2));
	database.angtype1=(char (*)[5])calloc(database.nangtype,sizeof(*database.angtype1));
	database.angtype2=(char (*)[5])calloc(database.nangtype,sizeof(*database.angtype2));
	database.angtype3=(char (*)[5])calloc(database.nangtype,sizeof(*database.angtype3));
//	sysdat.param_bnds=(int *)calloc(database.nbndtype,sizeof(int)); // changed by shuhei
//	sysdat.param_angs=(int *)calloc(database.nangtype,sizeof(int));
	sysdat.param_bnds=(int *)calloc(MaxNBondPerAtom*sysdat.nats,sizeof(int));
	sysdat.param_angs=(int *)calloc(MaxNBondPerAtom*sysdat.nats,sizeof(int));
	sysdat.coordx=(double *)calloc(sysdat.total_ats,sizeof(double));
	sysdat.coordy=(double *)calloc(sysdat.total_ats,sizeof(double));
	sysdat.coordz=(double *)calloc(sysdat.total_ats,sizeof(double));

	read_database(argv[argc-2],&database);

	printf("###########################\n");
	printf("####  DATABASE SUMMARY ####\n");
	printf("FOUND %d UNIQUE VDW PAIR PARAMS\n",database.nvdwtype);
	printf("FOUND %d UNIQUE BOND PARAMS\n",database.nbndtype);
	printf("FOUND %d UNIQUE ANGLE PARAMS\n",database.nangtype);


	read_pdb(argv[argc-1],&sysdat);
	get_unique(&database,topdat,&sysdat,ISCHARGED);	// write PARM.FILE

	read_coords(argv[argc-1],&database,topdat,&sysdat,ISCHARGED); // write DATA.FILE
	write_psf(argv[argc-1],&database,topdat,&sysdat,ISCHARGED);
	return 0;
}

