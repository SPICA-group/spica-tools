#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#define MAX 128
#define RECORD 80
#define BUFFER 83

#define PI 3.14159265358979

typedef struct sysdat
{
  int nats,nbnds,nangs,nimprops,ndiheds,ntops,total_ats,total_bnds,total_angs,total_improps,total_diheds,foundatoms,boxinfo;
  int uniq_nats,uniq_nbnds,uniq_nangs, uniq_nimprops,uniq_ndiheds;
  int *param_bnds,*param_angs;
  double *coordx,*coordy,*coordz;
  double boxx,boxy,boxz;
} SYSDAT;

typedef struct topdat
{
  int *bndndx1,*bndndx2,*bndtype,*angndx1,*angndx2,*angndx3,*angtype;
  int *improp_func,*impropndx1,*impropndx2,*impropndx3,*impropndx4, *improptype;
  int *dihed_func,*dihedndx1,*dihedndx2,*dihedndx3,*dihedndx4, *dihedtype, *dihedn;
  int *dihedpset, *improppset, *bndpset, *angpset;
  int *index, nat,nbnd,nang, nimprop, nmol,*parm_atomtype,ndihed, *dihedeq;
  double *dihedfk, *dihedof;
  double *mass,*charge, *bndfk, *bndeq, *angfk, *angeq, *impropfk, *impropeq;
  char (*type)[5],(*name)[5],(*segid)[5],(*resname)[5];
}TOPDAT;

typedef struct database
{
  double *fbnd,*bnde,*fang,*ange,*eps,*sig,*angsdk;
  int nvdwtype,nbndtype,nangtype;
  char (*vdwtype1)[5],(*vdwtype2)[5],(*vdwstyle)[7];
  char (*bndtype1)[5],(*bndtype2)[5];
  char (*angtype1)[5],(*angtype2)[5],(*angtype3)[5];
}DATABASE;

/* Function from setup_lammps*/
int strcmp_wc(const char *s1, const char *s2);
double absv(const double a[3]);
double inner(const double a[3], const double b[3]);
void outer(double axb[3], const double a[3], const double b[3]);
double GetAngle(const double *r1, const double *r2, const double *r3);
void read_pdb(char *filename,SYSDAT *sysdat);
void read_coords(DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat, int ISCHARGED);
void get_unique(DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat,int ISCHARGED);
void read_database(char *filename, DATABASE *database);
void count_params(char *filename, DATABASE *database);
void read_top(char *filename, TOPDAT *topdat, int ntop, int *ISCHARGED);
void count_atoms(char *filename, TOPDAT *topdat, int ntop);
int min(const int a, const int b);
void write_psf(DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat, int ISCHARGED);
void make_ndx(DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat, int ISCHARGED);
int find_solv(const char *s);
int find_prot(const char *s);
void make_top(DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat, int ISCHARGED);

/* Main routine. Call and allocate                                       */
/* The idea is to read in the topologies and then check the database for */
/* all of the required interaction params.                               */

int main(int argc, char **argv )
{
  int ntops,rdtp,i,ISCHARGED,opt;
  TOPDAT topdat[1000];
  DATABASE database;
  SYSDAT sysdat;
  memset(&sysdat,0,sizeof(SYSDAT));
  
  opterr = 0;
  //while((opt = getopt(argc,argv,"p"))!=-1){
  opt = getopt(argc,argv,"p");
  switch (opt) {
    case 'p':
      if( argc < 5){
        printf("DUMPS input files for a GROMACS run.\n");
        printf("USAGE: setup_gromacs -p <topfile 1> <nmol 1> [ <topfile 2> <nmol 2> ..... <topfile n> <nmol n>] <database> <pdbfile>\n");
        printf("Takes at least four arguments (one component system): 1) Topology, 2) number of molecules, 3) parameter database, 4) PDB.\n");
        exit(1);
      }
      printf("SETUP_GROMACS for SPICA PROTEIN model.\n");
      break;

    default:
      if( argc < 4 || argc%2 != 0 ){
        printf("DUMPS input files for a GROMACS run.\n");
        printf("USAGE: setup_gromacs [-p] <topfile 1> <nmol 1> [ <topfile 2> <nmol 2> ..... <topfile n> <nmol n>] <database>\n");
        printf("Takes at least three arguments (one component system): 1) Topology, 2) number of molecules, 3) parameter database.\n");
        exit(1);
      }
      printf("SETUP_GROMACS for SPICA.\n");
  }
  //printf("FOUND %d ARGS\n",argc);
  ISCHARGED=0;

/*   printf ("argc=%d\n",argc); 
   printf("0= %s\n",argv[0]); 
   printf("1= %s\n",argv[1]); 
   printf("2= %s\n",argv[2]); 
   printf("3= %s\n",argv[3]); 
   printf("4= %s\n",argv[4]); */
  switch (opt) {
    case 'p':
      ntops=(argc-3)/2;
      break;
    default:
      ntops=(argc-2)/2;
  }
  printf("WILL READ %d TOPOLOGY FILE(S).\n",ntops);
  rdtp=0;

/* Loop through the topologies and count the number of atoms, bonds and bends */
  
  switch (opt) {
    case 'p':
      while(rdtp < ntops){
        topdat[rdtp].nmol=atoi(argv[(2*rdtp)+3]);
        count_atoms(argv[(2*rdtp)+2],topdat,rdtp);
        rdtp++;
      }
      break;
    default:
      while(rdtp < ntops){
        topdat[rdtp].nmol=atoi(argv[(2*rdtp)+2]);
        count_atoms(argv[(2*rdtp)+1],topdat,rdtp);
        rdtp++;
      }
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
    printf("BOOKKEEPING:\n");
  //  printf("TOPFILE %s\n",argv[(2*i)+1]);
    printf("FOUND: %d atoms\n",topdat[i].nat);
    printf("FOUND: %d bonds\n",topdat[i].nbnd);
    printf("FOUND: %d angles\n",topdat[i].nang);
    printf("FOUND: %d improperss\n",topdat[i].nimprop);
    
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
    topdat[i].type    = (char (*)[5])calloc(topdat[i].nat,sizeof(*topdat[i].type));
    topdat[i].name    = (char (*)[5])calloc(topdat[i].nat,sizeof(*topdat[i].name));
    topdat[i].resname = (char (*)[5])calloc(topdat[i].nat,sizeof(*topdat[i].resname));
    topdat[i].segid   = (char (*)[5])calloc(topdat[i].nat,sizeof(*topdat[i].segid));

    sysdat.nats	     +=topdat[i].nat;
    sysdat.nbnds     +=topdat[i].nbnd;
    sysdat.nangs     +=topdat[i].nang;
    sysdat.nimprops  +=topdat[i].nimprop;
    sysdat.ndiheds   +=topdat[i].ndihed;

    sysdat.total_ats	 +=topdat[i].nat*topdat[i].nmol;
    sysdat.total_bnds	 +=topdat[i].nbnd*topdat[i].nmol;
    sysdat.total_angs	 +=topdat[i].nang*topdat[i].nmol;
    sysdat.total_improps +=topdat[i].nimprop*topdat[i].nmol;
    sysdat.total_diheds	 +=topdat[i].ndihed*topdat[i].nmol;
  }
/*   printf("BOOKKEEPING:\n"); */
/*   printf("TOTAL NONUNIQUE\n",i); */
/*   printf("FOUND: %d atoms\n",sysdat.nats); */
/*   printf("FOUND: %d bonds\n",sysdat.nbnds); */
/*   printf("FOUND: %d impropers\n",sysdat.total_improps); */
/*   printf("FOUND: %d dihderalss\n",sysdat.total_diheds); */

  rdtp=0;  
  switch(opt){
    case 'p':
      while(rdtp < ntops){
        topdat[rdtp].nmol=atoi(argv[(2*rdtp+3)]);
        read_top(argv[(2*rdtp)+2],topdat,rdtp,&ISCHARGED);
        rdtp++;
      }
      break;
    default:
      while(rdtp < ntops){
        topdat[rdtp].nmol=atoi(argv[(2*rdtp+2)]);
        read_top(argv[(2*rdtp)+1],topdat,rdtp,&ISCHARGED);
        rdtp++;
    }
  } 
  
  switch(opt){
    case 'p':
      count_params(argv[argc-2],&database);
      break;
    default:
      count_params(argv[argc-1],&database);
  }
  
  database.fbnd=(double *)calloc(database.nbndtype,sizeof(double));
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
  sysdat.param_bnds=(int *)calloc(database.nbndtype,sizeof(int));
  sysdat.param_angs=(int *)calloc(database.nangtype,sizeof(int));
  sysdat.coordx=(double *)calloc(sysdat.total_ats,sizeof(double));
  sysdat.coordy=(double *)calloc(sysdat.total_ats,sizeof(double));
  sysdat.coordz=(double *)calloc(sysdat.total_ats,sizeof(double));
  
  switch(opt){
    case 'p':
      read_database(argv[argc-2],&database);
      break;
    default:
      read_database(argv[argc-1],&database);
  }

  printf("FOUND %d UNIQUE VDW PAIR PARAMS\n",database.nvdwtype);
  printf("FOUND %d UNIQUE BOND PARAMS\n",database.nbndtype);
  printf("FOUND %d UNIQUE ANGLE PARAMS\n",database.nangtype);
  switch(opt){
    case 'p':
      if (strstr(argv[argc-1],".pdb" )){
        printf("Takes angles from %s.\n",argv[argc-1]);
        read_pdb(argv[argc-1],&sysdat);
      }
      else{
        printf("NO PDB FILES to take angles!\n");
	exit(1);
      }
      break;
  }

  
  get_unique(&database,topdat,&sysdat,ISCHARGED);
  read_coords(&database,topdat,&sysdat,ISCHARGED);
  make_top(&database,topdat,&sysdat,ISCHARGED);
  make_ndx(&database,topdat,&sysdat,ISCHARGED);
  write_psf(&database,topdat,&sysdat,ISCHARGED);

  return 0;
  
}
//
// END main
//
int strcmp_wc(const char *s1, const char *s2){
	if (s1[2]=='*' || s2[2]=='*'){
		return memcmp(s1, s2, 2);
	}else{
		return strcmp(s1, s2);
	}
}
double absv(const double a[3]){
	return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}
double inner(const double a[3], const double b[3]){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
void outer(double axb[3], const double a[3], const double b[3]){
	int d;
	for (d=0; d<3; d++){
		int d1=(d+1)%3;
		int d2=(d+2)%3;
		axb[d] = a[d1]*b[d2]-b[d1]*a[d2];
	}
}
double GetAngle(const double *r1, const double *r2, const double *r3){
	const double a[3] = {r1[0]-r2[0], r1[1]-r2[1], r1[2]-r2[2]};
        const double b[3] = {r3[0]-r2[0], r3[1]-r2[1], r3[2]-r2[2]};
        const double ab = inner(a,b);
        const double lallbl = absv(a)*absv(b);
        const double cos = ab/lallbl;
        if (cos<-1.0)     return PI;
        else if (cos>1.0) return 0;
        else                      return acos(cos);
}
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
			printf("FOUND BOXSIZE DATA\n");
			sscanf(&buf[7],"%lf	%lf	%lf", &sysdat->boxx,&sysdat->boxy,&sysdat->boxz);
		}
		else if (!strncmp(buf,"ATOM ",5) || !strncmp(buf,"HETATM",6)){
			if (sysdat->foundatoms >= sysdat->total_ats){
				printf("ERROR: foundatoms in pdb file >= total_ats in top files.\n");
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
		printf("ERROR: DID NOT FIND ANY ATOMS IN THE PDB FILE\n");
		exit(1);
	}
	if(sysdat->boxx==0.0){
		printf("WARNING: DID NOT FIND CELL SIZE!!\n");
		printf("BOX SIZE WILL HAVE BE SET BY HAND\n");
	}
	fclose(fpin);

}
/* */
void read_coords(DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat, int ISCHARGED)
{
  FILE *fpmol;
  int i,j,k,atindex,molindex,bondindex,angleindex,offset,impropindex,dihedindex,indcnt;
  double angle_in_pdb,eps,sig,tmpcalca,tmpcalcb;
  int ifound, datndx, vdwtmp, l, ikeep;
  int sol_flag,oth_flag;
  sol_flag = oth_flag = 0;
  //int *ang_params, *bnd_params;
  
/* will call modules depending on file type to be read */
/* right now only pdb capability */
  
  
  if((fpmol = fopen("molecule.itp","w")) == NULL)
    {
      fprintf(stderr,"ERROR: can't open molecule.itp\n");
      exit(1);
    }

  fprintf(fpmol,"; generated by setup_gromacs\n");
  fprintf(fpmol,"; RH DEVANE\n");

  
  atindex = 0;
  for(i=0;i<sysdat->ntops;i++) {
    fprintf(fpmol,"\n");
    fprintf(fpmol,"[ moleculetype ]\n");
    fprintf(fpmol,"; name nrexcl\n");
    fprintf(fpmol,"%s      2\n",topdat[i].resname[0]);
    fprintf(fpmol,"\n");
    fprintf(fpmol,"[ atoms ]\n");
    fprintf(fpmol,"; nr    type    resnr    residu   atom   cgnr   charge  mass\n");
    atindex=0;
    molindex=0;
    
    for(k=0;k<topdat[i].nat;k++){
      atindex++;
      fprintf(fpmol,"%6d %6s    1    %6s %6s %6d %8.4f\n",atindex,topdat[i].type[k],topdat[i].resname[k],topdat[i].name[k],k+1, topdat[i].charge[k]);
    }
    
    if(topdat[i].nbnd>0){
      fprintf(fpmol,"\n");
      fprintf(fpmol,"[ bonds ]\n");
      fprintf(fpmol,"\n");
      for(k=0;k<topdat[i].nbnd;k++){
        if ( topdat[i].bndpset[k] == 1 ) {
	  tmpcalca=topdat[i].bndeq[k]/10.0;
	  tmpcalcb=topdat[i].bndfk[k]*4.184*2.0*100;
	  fprintf(fpmol,"%5d %5d    6  %8.4f  %8.4f ; EN specified by bondparam \n",topdat[i].bndndx1[k],topdat[i].bndndx2[k], tmpcalca, tmpcalcb);
        }
	else if ( topdat[i].bndpset[k] == -1 ) {
	  fprintf(fpmol,"%5d %5d    1\n", topdat[i].bndndx1[k],topdat[i].bndndx2[k]);
        }
	else {
	  printf("NO BOND directive!\n");
	  exit(1);
	}
        bondindex++;
      }
    }

    if(topdat[i].nang>0){
      fprintf(fpmol,"\n");
      fprintf(fpmol,"[ angles ]\n");
      fprintf(fpmol,"\n");
      angleindex=0;
      offset=0;
      for(k=0;k<topdat[i].nang;k++){
        angleindex++;
	if(topdat[i].angpset[k] == -1){
/* now compare to the database */
          for(l=0;l<database->nangtype;l++){
            if(!strcmp_wc(database->angtype2[l], topdat[i].type[topdat[i].angndx2[k]-1]))
            {
              if(!strcmp_wc(database->angtype1[l], topdat[i].type[topdat[i].angndx1[k]-1])
                 && !strcmp_wc(database->angtype3[l], topdat[i].type[topdat[i].angndx3[k]-1])){
                datndx=l;
                l=database->nangtype;
              }
              if(!strcmp_wc(database->angtype3[l], topdat[i].type[topdat[i].angndx1[k]-1])
                 && !strcmp_wc(database->angtype1[l], topdat[i].type[topdat[i].angndx3[k]-1])){
                datndx=l;
                l=database->nangtype;
              }
            }
          }
         // datndx is determined , changing -1 to something
          ifound=0;
          for(l=0;l<database->nvdwtype;l++){
            if(!strcmp_wc(database->vdwtype1[l],topdat[i].type[topdat[i].angndx1[k]-1]) && !strcmp_wc(database->vdwtype2[l],
                  topdat[i].type[topdat[i].angndx3[k]-1])) {
              ifound=1;
              vdwtmp=l;
              l=database->nvdwtype;
            }
            else if (!strcmp_wc(database->vdwtype1[l],topdat[i].type[topdat[i].angndx3[k]-1]) && !strcmp_wc(database->vdwtype2[l],
                topdat[i].type[topdat[i].angndx1[k]-1])) {
              ifound=1;
              vdwtmp=l;
              l=database->nvdwtype;
            }
          }
	  eps=database->eps[vdwtmp]*4.184;
	  sig=database->sig[vdwtmp]/10.0;
	  tmpcalca=database->ange[datndx];
	  tmpcalcb=database->fang[datndx]*4.184*2.0;
	  if (tmpcalca < 0){
	    const int i1 = topdat[i].angndx1[k]-1;
	    const int i2 = topdat[i].angndx2[k]-1;
	    const int i3 = topdat[i].angndx3[k]-1;
	    const double r1[3] = {sysdat->coordx[i1],sysdat->coordy[i1],sysdat->coordz[i1]};
	    const double r2[3] = {sysdat->coordx[i2],sysdat->coordy[i2],sysdat->coordz[i2]};
	    const double r3[3] = {sysdat->coordx[i3],sysdat->coordy[i3],sysdat->coordz[i3]};
	    angle_in_pdb = 180.0/PI*GetAngle(r1,r2,r3);
            fprintf(fpmol,"%5d %5d %5d    1  %8.4f %8.4f ; %6s %6s %6s taken from pdb\n", topdat[i].angndx1[k],topdat[i].angndx2[k],topdat[i].angndx3[k], angle_in_pdb, tmpcalcb, database->angtype1[datndx],database->angtype2[datndx],database->angtype3[datndx]);
	  }
	  else{
	   	if (database->angsdk[datndx] == 0){
	    		fprintf(fpmol,"%5d %5d %5d    1\n",topdat[i].angndx1[k],topdat[i].angndx2[k],topdat[i].angndx3[k]);
		} else {
	    		fprintf(fpmol,"%5d %5d %5d    5\n",topdat[i].angndx1[k],topdat[i].angndx2[k],topdat[i].angndx3[k]);
		}
          }
        }
      }	
      
    }

    if(topdat[i].ndihed>0){
      fprintf(fpmol,"\n");
      fprintf(fpmol,"[ dihedrals ]\n");
      fprintf(fpmol,"\n");
      dihedindex=0;
      offset=0;
      for(k=0;k<topdat[i].ndihed;k++){
	fprintf(fpmol,"%5d %5d %5d %5d  1  %-3d %8.4f %-3d ; FROM TOP \n",topdat[i].dihedndx1[k],topdat[i].dihedndx2[k],topdat[i].dihedndx3[k],topdat[i].dihedndx4[k],
			topdat[i].dihedeq[k],topdat[i].dihedfk[k]*4.184,topdat[i].dihedn[k]);
      }
    }
    /* if ther are impropers....print them out */
  
    if(topdat[i].nimprop>0){
      fprintf(fpmol,"\n");
      fprintf(fpmol,"[ impropers ]\n");
      fprintf(fpmol,"\n");
      impropindex=0;
      offset=0;
      for(k=0;k<topdat[i].nimprop;k++){
	fprintf(fpmol,"%5d %5d %5d %5d  2  %8.4f %8.4f ; FROM TOP\n",topdat[i].impropndx1[k],topdat[i].impropndx2[k],topdat[i].impropndx3[k],topdat[i].impropndx4[k],
			topdat[i].impropfk[k]*4.184*2.0,topdat[i].impropeq[k]);
      }
    }
  }
  fclose(fpmol);
  
}


/* Decide which params will be needed */


void get_unique(DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat,int ISCHARGED)
{
  FILE *fpout;
  int i,j,k;
  char (*uniq_atype)[5],(*bnd_name1)[5],(*bnd_name2)[5],(*b1tmp)[5],(*b2tmp)[5];
  double *uniq_mass,*uniq_charge;
  int uniq_nats,uniq_bnds,uniq_angs,uniq_improps,uniq_diheds,ikeep,ifound,vdwtmp;
  int *bnd_params,*ang_params,*ang_vdw,datndx;
  double tmpcalca, tmpcalcb;
  double sig, sigsq, sigcub, eps, disp, repul, pf124, pf96, pf125, pf126, c5, c4, c12, c6, c9;
  
  if((fpout = fopen("SPICA.itp","w")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open SPICA.itp\n");
    exit(1);
  }
  
  uniq_atype=(char (*)[5])calloc(sysdat->nats,sizeof(*uniq_atype)); 
  uniq_mass=(double *)calloc(sysdat->nats,sizeof(double));
  uniq_charge=(double *)calloc(sysdat->nats,sizeof(double));
  bnd_params=(int *)calloc(database->nbndtype,sizeof(int));
  ang_params=(int *)calloc(database->nangtype,sizeof(int));
  ang_vdw=(int *)calloc(database->nangtype,sizeof(int));
  bnd_name1=(char (*)[5])calloc(database->nbndtype,sizeof(*bnd_name1));
  bnd_name2=(char (*)[5])calloc(database->nbndtype,sizeof(*bnd_name2));
  b1tmp=(char (*)[5])calloc(database->nbndtype,sizeof(*b1tmp));
  b2tmp=(char (*)[5])calloc(database->nbndtype,sizeof(*b2tmp));

  fprintf(fpout,"; Generated by setup_gromacs\n");
  fprintf(fpout,"; RH DEVANE\n");
  fprintf(fpout,"\n");


  /* first gather the unique atom types */
  
  uniq_nats=0;
  for(i=0;i<sysdat->ntops;i++){
    for(j=0;j<topdat[i].nat;j++){
      ikeep=1;
      for(k=0;k<uniq_nats;k++){
        if(!strcmp(topdat[i].type[j],uniq_atype[k]))
        {
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
  
  sysdat->uniq_nats=uniq_nats;
  fprintf(fpout,"[ defaults ]\n");
  fprintf(fpout,"; nbfun    comb-rule    gen-pairs    fudgeLJ fudgeQQ\n");
  fprintf(fpout,"1           1             no\n");
  fprintf(fpout,"\n");
  fprintf(fpout,"[ atomtypes ]\n");
  fprintf(fpout,"; name   mass   charge   ptype   sigma   epsilon\n");

  /* I SET THE VDW PARAMS TO ZERO HERE BECAUSE I PUT THEM IN PAIRS...ALSO THIS MAY PREVENT */
  /* GROMACS FROM TRYING TO GENERATE PAIRS....THERE IS ALOS A SETTING TO MAKE SURE GROMACS */
  /* DOES NOT GENERATE PAIRS, BUT IF THAT IS NOT SET PROPERLY THIS WOULD CAUSE A CRASH AT WORSE. */

  for(i=0;i<uniq_nats;i++){
    fprintf(fpout,"%8s %8.4f %8.4f    A    0.0    0.0\n",uniq_atype[i], uniq_mass[i], uniq_charge[i] );
  }

  fprintf(fpout,"\n");
  
/* get pair interactions */
  
  fprintf(fpout,"[ nonbond_params ]\n");
  fprintf(fpout,"; i     j    func   C    A\n");
  
  for(i=0;i<uniq_nats;i++){
    for(j=i;j<uniq_nats;j++){
      ifound=0;
      for(k=0;k<database->nvdwtype;k++){
        if(!strcmp(database->vdwtype1[k],uniq_atype[i]) && !strcmp(database->vdwtype2[k],uniq_atype[j])) {
/*           printf("Found database entry\n"); */
          ifound=1;
          vdwtmp=k;
          k=database->nvdwtype;
        }
        else if (!strcmp(database->vdwtype2[k],uniq_atype[i]) && !strcmp(database->vdwtype1[k],uniq_atype[j])) {
/*           printf("Found database entry\n"); */
          ifound=1;
          vdwtmp=k;
          k=database->nvdwtype;
        }
      }
      if(ifound==0){
        printf("*********************\n");
        printf("ERROR:No params for VDW interaction between %s and %s\n",uniq_atype[i],uniq_atype[j]);
        printf("UPDATE DATABASE!!!\n");
        exit(1);
      } else if (ifound==1){
	/* convert our params to GROMACS C6 C12 equivalents */
	eps=database->eps[vdwtmp];
	sig=database->sig[vdwtmp]/10.0;
	sigsq=sig*sig;
	sigcub=sig*sig*sig;

        if(!strcmp( database->vdwstyle[vdwtmp], "lj12_4")) {
	  pf124=4.184*3.0*1.73205080757/2.0;
	  disp=pf124*eps*sigsq*sigsq;
	  repul=pf124*eps*sigcub*sigcub*sigcub*sigcub;
	  //fprintf(fpout,"; %s %5.6f %5.6f\n",database->vdwstyle[vdwtmp], eps, sig);
	} else if(!strcmp( database->vdwstyle[vdwtmp], "lj9_6")) {
	  pf96=4.184*27.0/4.0;
	  disp=pf96*eps*sigcub*sigcub;
	  repul=pf96*eps*sigcub*sigcub*sigcub;
	 // fprintf(fpout,"; %s %5.6f %5.6f\n",database->vdwstyle[vdwtmp], eps, sig);
	} else if(!strcmp( database->vdwstyle[vdwtmp], "lj12_5")) {
	  pf125=4.184*3.203779841;
	  disp=pf125*eps*sigcub*sigsq;
	  repul=pf125*eps*sigcub*sigcub*sigcub*sigcub;
	} else if(!strcmp( database->vdwstyle[vdwtmp], "lj12_6")) {
	  pf126=4.184*4.0;
	  disp=pf126*eps*sigcub*sigcub;
	  repul=pf126*eps*sigcub*sigcub*sigcub*sigcub;
	} else {
	  printf("ERROR: WRITE CORRECT LJform %s %s (e.g.)lj12_4\n",&database->vdwtype1[vdwtmp],&database->vdwtype2[vdwtmp]);
	  exit(0);
	}
	/* 	&database->vdwstyle[vdwtmp] */
	/*         printf("Found database entry for VDW\n"); */
        fprintf(fpout,"%6s %6s    1  %14.5e %14.5e ;  %s %5.6f %5.6f \n",&database->vdwtype1[vdwtmp], &database->vdwtype2[vdwtmp], disp, repul,database->vdwstyle[vdwtmp],eps,sig);
      }
    }
  }

  /* BONDSSSSSSSSSS */
  
/* get bond interactions  */
  fprintf(fpout,"\n");

  if(sysdat->nbnds > 0) {
    uniq_bnds=0;
    fprintf(fpout,"[ bondtypes ]\n");
    fprintf(fpout,"; i     j   funct   length  force.c\n");
    

    for(i=0;i<sysdat->ntops;i++){
      for(j=0;j<topdat[i].nbnd;j++){
        ikeep=1;
        datndx=-1;
        
        /* AT THIS POINT WE WILL CHECK TO SEE IF THE PARAMS WERE GIVEN */
        /* IN THE TOP FILE....IF SO WE WILL SKIP A LOT OF THIS AND ADD */
        /* THIS AS A UNIQUE BOND....IF NOT WE GO THROUGH THE PROCEDURE */
        if ( topdat[i].bndpset[j] == -1 ) {
          
/* now compare to the database */
          for(k=0;k<database->nbndtype;k++){
            if(!strcmp_wc(database->bndtype1[k], topdat[i].type[topdat[i].bndndx1[j]-1])
               && !strcmp_wc(database->bndtype2[k], topdat[i].type[topdat[i].bndndx2[j]-1])){
              datndx=k;
              k=database->nbndtype;
	      strcpy(b1tmp[j],topdat[i].type[topdat[i].bndndx1[j]-1]);
	      strcpy(b2tmp[j],topdat[i].type[topdat[i].bndndx2[j]-1]);
            }
            if(!strcmp_wc(database->bndtype2[k], topdat[i].type[topdat[i].bndndx1[j]-1])
               && !strcmp_wc(database->bndtype1[k], topdat[i].type[topdat[i].bndndx2[j]-1])){
              datndx=k;
              k=database->nbndtype;
	      strcpy(b2tmp[j],topdat[i].type[topdat[i].bndndx1[j]-1]);
	      strcpy(b1tmp[j],topdat[i].type[topdat[i].bndndx2[j]-1]);
            }
          }
/* did not find any params!!! TIME TO DIE  */
          if(datndx==-1){
            printf("ERROR: DID NOT FIND BOND PARAMETERS IN DATABASE %d %d %s %s\n",topdat[i].bndndx1[j],
               topdat[i].bndndx2[j],topdat[i].type[topdat[i].bndndx1[j]-1],topdat[i].type[topdat[i].bndndx2[j]-1]);
            printf("%d %s\n",topdat[i].bndndx1[j],topdat[i].type[0]);
            
            exit(1);
          }
          
/* Now make sure we do not already know we have this interaction */
          
          for(k=0;k<uniq_bnds;k++){
            if(!strcmp(bnd_name1[k], b1tmp[j])
               && !strcmp(bnd_name2[k], b2tmp[j])){
	      ikeep=0;
              topdat[i].bndtype[j]=k;
              k=uniq_bnds;
	    }
            if(!strcmp(bnd_name2[k], b1tmp[j])
               && !strcmp(bnd_name1[k], b2tmp[j])){
	      ikeep=0;
              topdat[i].bndtype[j]=k;
              k=uniq_bnds;
	    }
//
           // if(datndx==bnd_params[k])
           // {
           //   ikeep=0;
           //   topdat[i].bndtype[j]=k;
           //   k=uniq_bnds;
           // }
          }
          /* IKEEP = 1 IF WE FOUND A NEW ONE */
          if(ikeep==1) {
            bnd_params[uniq_bnds]=datndx;
	    strcpy(bnd_name1[uniq_bnds],b1tmp[j]);
	    strcpy(bnd_name2[uniq_bnds],b2tmp[j]);
            sysdat->param_bnds[uniq_bnds]=datndx;
            topdat[i].bndtype[j]=uniq_bnds;
            uniq_bnds++;
	    tmpcalca=database->bnde[bnd_params[uniq_bnds-1]]/10.0;
	    tmpcalcb=database->fbnd[bnd_params[uniq_bnds-1]]*4.184*2.0*100;
	    printf("%s %s\n",database->bndtype1[bnd_params[uniq_bnds-1]],database->bndtype2[bnd_params[uniq_bnds-1]]);
	    //if ( (strcmp(database->bndtype1[bnd_params[uniq_bnds-1]],"GB*")==0 && database->bndtype2[bnd_params[uniq_bnds-1]][2] != '*') || 
	    //     (strcmp(database->bndtype2[bnd_params[uniq_bnds-1]],"GB*")==0 && database->bndtype1[bnd_params[uniq_bnds-1]][2] != '*') ){
	    //  fprintf(fpout,"   GBT %6s    1  %8.4f  %8.4f\n",database->bndtype2[bnd_params[uniq_bnds-1]], tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   GBM %6s    1  %8.4f  %8.4f\n",database->bndtype2[bnd_params[uniq_bnds-1]], tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   GBB %6s    1  %8.4f  %8.4f\n",database->bndtype2[bnd_params[uniq_bnds-1]], tmpcalca, tmpcalcb);
	    //} else if (strcmp(database->bndtype1[bnd_params[uniq_bnds-1]],"GB*")==0 && strcmp(database->bndtype2[bnd_params[uniq_bnds-1]],"GB*")==0){
	    //  fprintf(fpout,"   GBT    GBT    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   GBM    GBT    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   GBB    GBT    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   GBM    GBM    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   GBB    GBM    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   GBB    GBB    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //} else if ( (strcmp(database->bndtype1[bnd_params[uniq_bnds-1]],"GB*")==0 && strcmp(database->bndtype2[bnd_params[uniq_bnds-1]],"AB*")==0) || 
	    //            (strcmp(database->bndtype1[bnd_params[uniq_bnds-1]],"AB*")==0 && strcmp(database->bndtype2[bnd_params[uniq_bnds-1]],"GB*")==0) ){
	    //  fprintf(fpout,"   GBT    ABT    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   GBM    ABT    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   GBB    ABT    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   GBT    ABB    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   GBM    ABB    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   GBB    ABB    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //} else if (strcmp(database->bndtype1[bnd_params[uniq_bnds-1]],"AB*")==0 && strcmp(database->bndtype2[bnd_params[uniq_bnds-1]],"AB*")==0){
	    //  fprintf(fpout,"   ABT    ABT    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   ABB    ABT    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //  fprintf(fpout,"   ABB    ABB    1  %8.4f  %8.4f\n", tmpcalca, tmpcalcb);
	    //} else {
	    fprintf(fpout,"%6s %6s    1  %8.4f  %8.4f\n",b1tmp[j],b2tmp[j], tmpcalca, tmpcalcb);
	    //fprintf(fpout,"%6s %6s    1  %8.4f  %8.4f\n",database->bndtype1[bnd_params[uniq_bnds-1]],database->bndtype2[bnd_params[uniq_bnds-1]], tmpcalca, tmpcalcb);
	    //}
          }
        } 
	/*else {
          // THE PARAMS WERE GIVEN IN THE TOP FILE SO LETS ADD IT TO THE PARAM FILE 
          topdat[i].bndtype[j]=uniq_bnds;
          uniq_bnds++;
	  tmpcalca=topdat[i].bndeq[j]/10.0;
	  tmpcalcb=topdat[i].bndfk[j]*4.184*2.0*100;
          fprintf(fpout,"%6s %6s    1  %8.4f  %8.4f\n",topdat[i].type[topdat[i].bndndx1[j]-1], topdat[i].type[topdat[i].bndndx2[j]-1], tmpcalca, tmpcalcb);

        }*/
      }
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
    fprintf(fpout,"[ angletypes ]\n");
    fprintf(fpout,"; i     j     k    funct   angle  force.c  sigma   eps\n");
    
    uniq_angs=0;
    
    for(i=0;i<sysdat->ntops;i++){
      for(j=0;j<topdat[i].nang;j++){
        ikeep=1;
        datndx=-1;

        /* AT THIS POINT WE WILL CHECK TO SEE IF THE PARAMS WERE GIVEN */
        /* IN THE TOP FILE....IF SO WE WILL SKIP A LOT OF THIS AND ADD */
        /* THIS AS A UNIQUE BOND....IF NOT WE GO THROUGH THE PROCEDURE */

        if ( topdat[i].angpset[j] == -1 ) {
          
/* now compare to the database */
        
          for(k=0;k<database->nangtype;k++){
            if(!strcmp_wc(database->angtype2[k], topdat[i].type[topdat[i].angndx2[j]-1]))
            {
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
         // datndx is determined , changing -1 to something
	 //
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
	 // when ifound=1, vdwtmp is determined
/* 	  printf("%d %d %d %d\n",ifound,topdat[i].angndx1[j],topdat[i].angndx2[j],topdat[i].angndx3[j]); */
          if(ifound==0){
            printf("*********************\n");
            printf("ERROR:No params for VDW interaction between %s and %s for angle (database)\n",topdat[i].type[topdat[i].angndx1[j]-1],
               topdat[i].type[topdat[i].angndx3[j]-1]);
            printf("UPDATE DATABASE!!!\n");
            exit(1);
          }
          
/* end VDW for CG angles */
          
/* No params for this interaction in the database */
          
          if(datndx==-1){
            printf("ERROR: DID NOT FIND ANGLE PARAMETERS IN DATABASE %s %s %s (%d %d %d )\n",topdat[i].type[topdat[i].angndx1[j]-1],
		   topdat[i].type[topdat[i].angndx2[j]-1],topdat[i].type[topdat[i].angndx3[j]-1],topdat[i].angndx1[j],topdat[i].angndx2[j],topdat[i].angndx3[j] );
            exit(1);
          }
          
/* Now make sure we do not already have this one */
          //currently, ikeep is 1, which means we dont have this one. 
          for(k=0;k<uniq_angs;k++){
            if(datndx==ang_params[k])
	      {
		ikeep=0;
		topdat[i].angtype[j]=k;
		k=uniq_angs;
	      }
          }
	  //ikeep is 1, which means we dont have this one.
          if(ikeep==1) {
	    eps=database->eps[vdwtmp]*4.184;
	    sig=database->sig[vdwtmp]/10.0;
            ang_params[uniq_angs]=datndx;
            ang_vdw[uniq_angs]=vdwtmp;
            sysdat->param_angs[uniq_angs]=datndx;
            topdat[i].angtype[j]=uniq_angs;
            uniq_angs++;
	    tmpcalca=database->ange[ang_params[uniq_angs-1]];
	    tmpcalcb=database->fang[ang_params[uniq_angs-1]]*4.184*2.0;
	    if (tmpcalca > 0)
	   	if (database->angsdk[datndx] == 0){
	    		fprintf(fpout,"%6s %6s %6s    1  %8.4f %8.4f\n",database->angtype1[ang_params[ uniq_angs-1]],database->angtype2[ang_params[ uniq_angs-1]],database->angtype3[ang_params[ uniq_angs-1]], tmpcalca, tmpcalcb);
		} else {
              		fprintf(fpout,"%6s %6s %6s    5  %8.4f %8.4f %8.4f %8.4f\n", database->angtype1[ang_params[ uniq_angs-1]],database->angtype2[ang_params[ uniq_angs-1]],
		    database->angtype3[ang_params[ uniq_angs-1]], tmpcalca, tmpcalcb, sig,eps);
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
            printf("*********************\n");
            printf("ERROR:No params for VDW interaction between %s and %s for angle (topfile)\n",topdat[i].type[topdat[i].angndx1[j]-1],
               topdat[i].type[topdat[i].angndx3[j]-1]);
            printf("UPDATE DATABASE!!!\n");
            exit(1);
          }
          
/* end VDW for CG angles */

          topdat[i].angtype[j]=uniq_angs;
          uniq_angs++;
	  eps=database->eps[vdwtmp]*4.184;
	  sig=database->sig[vdwtmp]/10.0;
	  tmpcalca=topdat[i].angeq[j];
	  tmpcalcb=topdat[i].angfk[j]*4.184*2.0;
	  if (tmpcalca > 0)
            fprintf(fpout,"%6s %6s %6s    5  %8.4f %8.4f %8.4f %8.4f ; FROM TOP \n",topdat[i].type[topdat[i].angndx1[j]-1],topdat[i].type[topdat[i].angndx2[j]-1],
		  topdat[i].type[topdat[i].angndx3[j]-1], tmpcalca, tmpcalcb,sig,eps);
        }
      }
    }
  }
  sysdat->uniq_nangs=uniq_angs;

/* NOW LET HANDLE THE DIHEDRAL PARAMS THE ONLY WAY WE DO...THEY HAVE TO BE SPECIFIED IN THE TOP FILE */
  fprintf(fpout,"\n");
  fclose(fpout);

#if 0
  if(sysdat->total_diheds > 0) {
    uniq_diheds=0;
    fprintf(fpout,"[ dihedraltypes ]\n");
    fprintf(fpout,"; i     j     k     l    funct   phi  force.c  multiplicity\n");
    
    for(i=0;i<sysdat->ntops;i++){
      for(j=0;j<topdat[i].ndihed;j++){

        if ( topdat[i].dihedpset[j] == 1 ) {
          topdat[i].dihedtype[j]=uniq_diheds;         
          uniq_diheds++;
	  ikeep = 1;
	  /* RHDDDD */
         /* for(k=0;k<database->nangtype;k++){
            if(!strcmp(database->angtype2[k], topdat[i].type[topdat[i].angndx2[j]-1]))
            {
              if(!strcmp(database->angtype1[k], topdat[i].type[topdat[i].angndx1[j]-1])
                 && !strcmp(database->angtype3[k], topdat[i].type[topdat[i].angndx3[j]-1])){
                datndx=k;
                k=database->nangtype;
              }
              if(!strcmp(database->angtype3[k], topdat[i].type[topdat[i].angndx1[j]-1])
                 && !strcmp(database->angtype1[k], topdat[i].type[topdat[i].angndx3[j]-1])){
                datndx=k;
                k=database->nangtype;
              }
            }
          }
          //currently, ikeep is 1, which means we dont have this one. 
          for(k=0;k<uniq_diheds;k++){
            if(datndx==ang_params[k])
	      {
		ikeep=0;
		topdat[i].dihedtype[j]=k;
		k=uniq_diheds;
	      }
          }*/
          if (ikeep == 1){
          fprintf(fpout,"%6s %6s %6s %6s  1  %-3d %8.4f %-3d \n",topdat[i].type[topdat[i].dihedndx1[j]-1],topdat[i].type[topdat[i].dihedndx2[j]-1],
		  topdat[i].type[topdat[i].dihedndx3[j]-1],topdat[i].type[topdat[i].dihedndx4[j]-1], topdat[i].dihedeq[j],topdat[i].dihedfk[j]*4.184,topdat[i].dihedn[j]);
          }
        }
      }
    }
  }
  sysdat->uniq_ndiheds=uniq_diheds;

/* NOW LET HANDLE THE IMPROPS PARAMS THE ONLY WAY WE DO...THEY HAVE TO BE SPECIFIED IN THE TOP FILE */
  fprintf(fpout,"\n");
  
  printf("IMPROPS TEST %d\n",sysdat->total_improps);
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
#endif
  
}


/* Read the database file and store unique params */
/* Warn if you find duplicates */


void read_database(char *filename, DATABASE *database) 
{
  FILE *fpin;
  int c,nvdw,nbnd,nang,i,ikeep;
  char col1[256],vdwtype1[5],vdwtype2[5],vdwstyle[7];
  char bndtype1[5],bndtype2[5],angtype1[5],angtype2[5],angtype3[5];
  double eps,sig,fbnd,bnde,fang,ange,angsdk;
  char line[350];

  if((fpin = fopen(filename,"r")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open DATABASE to read: %s\n",&filename);
    exit(1);
  }
  nvdw=0;
  nbnd=0;
  nang=0;

  while(c=fscanf(fpin,"%s",&col1)!= EOF)
  {
    if(col1[0]=='#')
	    continue;
    if(!strcmp((char *)&col1, "pair"))
    {
      fscanf(fpin,"%[^\n]",line);

      if(strchr(line,'#'))
        sscanf(line,"%s %s %s %lf %lf %*[^\n]",&vdwtype1,&vdwtype2,&vdwstyle,&eps,&sig);
      else
        sscanf(line,"%s %s %s %lf %lf",&vdwtype1,&vdwtype2,&vdwstyle,&eps,&sig);

      ikeep=1;
      for(i=0;i<nvdw;i++){
        if(!strcmp(vdwtype1,database->vdwtype1[i]) && !strcmp(vdwtype2,database->vdwtype2[i]))
        {
          printf("WARNING: FOUND DUP VDW PARAM %s %s\n",&vdwtype1,&vdwtype2);
          ikeep=0;
        }
        else if (!strcmp(vdwtype1,database->vdwtype2[i]) && !strcmp(vdwtype2,database->vdwtype1[i]))
        {
          printf("WARNING: FOUND DUP VDW PARAM %s %s\n",&vdwtype1,&vdwtype2);
          ikeep=0;
        }
      }
      if(ikeep==1)
      {
        strcpy(database->vdwtype1[nvdw],vdwtype1);
        strcpy(database->vdwtype2[nvdw],vdwtype2);
        strcpy(database->vdwstyle[nvdw],vdwstyle);
        database->eps[nvdw]=eps;
        database->sig[nvdw]=sig;
        nvdw++;
      }
    }
    if(!strcmp((char *)&col1, "bond"))
    {

      fscanf(fpin,"%[^\n]",line);
      
      if(strchr(line,'#'))
        sscanf(line,"%s %s %lf %lf %*[^\n]",&bndtype1,&bndtype2,&fbnd,&bnde);
      else
        sscanf(line,"%s %s %lf %lf",&bndtype1,&bndtype2,&fbnd,&bnde);

      ikeep=1;
      for(i=0;i<nbnd;i++){
        if(!strcmp(bndtype1,database->bndtype1[i]) && !strcmp(bndtype2,database->bndtype2[i]))
        {
          printf("WARNGIN: FOUND DUP BOND PARAM %s %s\n",&bndtype1,&bndtype2);
          ikeep=0;
        }
        else if (!strcmp(bndtype1,database->bndtype2[i]) && !strcmp(bndtype2,database->bndtype1[i]))
        {
          printf("WARNING: FOUND DUP BOND PARAM %s %s\n",&bndtype1,&bndtype2);
          ikeep=0;
        }
      }
      if(ikeep==1)
      {
        strcpy(database->bndtype1[nbnd],bndtype1);
        strcpy(database->bndtype2[nbnd],bndtype2);
        database->fbnd[nbnd]=fbnd;
        database->bnde[nbnd]=bnde;
        nbnd++;
      }
    }

    if(!strcmp((char *)&col1, "angle"))
    {

      fscanf(fpin,"%[^\n]",line);
      if(strstr(line, "harmonic") != NULL) {
	      angsdk = 0;
      }
      else {
	      angsdk = 1;
      }
      if(strchr(line,'#'))
        sscanf(line,"%s %s %s %lf %lf %*[^\n]",&angtype1,&angtype2,&angtype3,&fang,&ange);
      else
        sscanf(line,"%s %s %s %lf %lf",&angtype1,&angtype2,&angtype3,&fang,&ange);

      ikeep=1;
      for(i=0;i<nang;i++){
        if(!strcmp(angtype2,database->angtype2[i])) {
          if(!strcmp(angtype1,database->angtype1[i]) && !strcmp(angtype3,database->angtype3[i]))
          {
            printf("WARNGIN: FOUND DUP ANGLE PARAM %s %s %s\n",&angtype1,&angtype2,&angtype3);
            ikeep=0;
          }
          else if (!strcmp(angtype3,database->angtype1[i]) && !strcmp(angtype1,database->angtype3[i]))
          {
            printf("WARNING: FOUND DUP ANGLE PARAM %s %s %s\n",&angtype1,&angtype2,&angtype3);
            ikeep=0;
          }
        }
      }
      if(ikeep==1)
      {
        strcpy(database->angtype1[nang],angtype1);
        strcpy(database->angtype2[nang],angtype2);
        strcpy(database->angtype3[nang],angtype3);
        database->fang[nang]=fang;
        database->ange[nang]=ange;
	database->angsdk[nang]=angsdk;
        nang++;
	//printf("%d\n",nang);
	//printf("%d\n",database->nangtype);
      }
    }
  }
  fclose(fpin);
  
  database->nvdwtype=nvdw;
  database->nbndtype=nbnd;
  database->nangtype=nang;
  
}


/* Count the number of params in the database so we can allocate for storage */


void count_params(char *filename, DATABASE *database)
{
  FILE *fpin;
  int c;
  char col1[256];
  
  if((fpin = fopen(filename,"r")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open DATABASE to parse %s\n",&filename);
    exit(1);
  }
  database->nvdwtype=0;
  database->nbndtype=0;
  database->nangtype=0;
  
  while(c=fscanf(fpin,"%s",&col1)!= EOF)
  {
    if(col1[0]=='#')
	    continue;
    if(!strcmp((char *)&col1, "pair"))
    {

      fscanf(fpin,"%*[^\n]");
      database->nvdwtype++;
    }
    if(!strcmp((char *)&col1, "bond"))
    {

      fscanf(fpin,"%*[^\n]");
      database->nbndtype++;
    }
    if(!strcmp((char *)&col1, "angle"))
    {
      
      fscanf(fpin,"%*[^\n]");
      database->nangtype++;
    }
  }
  fclose(fpin);
}


/* Read the topology file and store the data */


void read_top(char *filename, TOPDAT *topdat, int ntop, int *ISCHARGED)
{
  FILE *fpin;
  char *line;
  int c, index,ndx,bndx,andx,indx,dndx;
  const int Nbuf=256;
  double tmpfl, tmp1, tmp2;
  char buf[Nbuf];
  char col1[Nbuf];
  
  if((fpin = fopen(filename,"r")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open infile %s\n",&filename);
    exit(1);
  }
  ndx=0;
  bndx=0;
  andx=0;
  indx=0;
  dndx=0;

  printf("READ %s\n", filename);
  
  //while(c=fscanf(fpin,"%s",&col1)!= EOF)
  while(fgets(buf,Nbuf,fpin))
  {
    strcpy(col1,""); sscanf(buf,"%s",col1);
    line = &buf[strlen(col1)];
    if(!strcmp((char *)&col1, "atom"))
    {
      sscanf(line,"%d %s %s %s %lf %lf %s",&topdat[ntop].index[ndx], &topdat[ntop].resname[ndx], &topdat[ntop].name[ndx], &topdat[ntop].type[ndx],&topdat[ntop].mass[ndx],
         &topdat[ntop].charge[ndx],&topdat[ntop].segid[ndx]);
      if((topdat[ntop].charge[ndx]*topdat[ntop].charge[ndx]) > 0.00001){
        printf("CHARGE IN TOP FILE: %s %lf\n", filename, topdat[ntop].charge[ndx]);
        *ISCHARGED=1;
      }
      
        ndx++;
    }
    
    if(!strcmp((char *)&col1, "bond"))
    {
      sscanf(line,"%d %d", &topdat[ntop].bndndx1[bndx], &topdat[ntop].bndndx2[bndx]);
      topdat[ntop].bndpset[bndx] = -1;
      bndx++;
    }
    
    if(!strcmp((char *)&col1, "bondparam"))
    {
      printf("WARNING!!!!!!!: USING BOND PARAMETERS FROM THE TOP FILE\n");
      if(sscanf(line,"%d %d %lf %lf", &topdat[ntop].bndndx1[bndx], &topdat[ntop].bndndx2[bndx], &topdat[ntop].bndfk[bndx],
            &topdat[ntop].bndeq[bndx]) != 4) {
        
        printf("ERROR: NOT ENOUGH ARGS FOR BONDPARAM: MUST BE: NDX1 NDX2 FK EQ\n");
        exit(1);
      }
      
      topdat[ntop].bndpset[bndx] = 1;
      bndx++;
    }
    
    if(!strcmp((char *)&col1, "angle"))
    {
      sscanf(line,"%d %d %d", &topdat[ntop].angndx1[andx], &topdat[ntop].angndx2[andx], &topdat[ntop].angndx3[andx]);
      topdat[ntop].angpset[andx] = -1;
      andx++;
    }
    
    if(!strcmp((char *)&col1, "angleparam"))
    {
      printf("WARNING!!!!!!!: USING ANGLE PARAMETERS FROM THE TOP FILE\n");
      if(sscanf(line,"%d %d %d %lf %lf", &topdat[ntop].angndx1[andx], &topdat[ntop].angndx2[andx], &topdat[ntop].angndx3[andx],
            &topdat[ntop].angfk[andx],&topdat[ntop].angeq[andx]) != 5) {
        printf("ERROR: NOT ENOUGH ARGS FOR ANGLEPARAM: MUST BE: NDX1 NDX2 NDX3 FK EQ\n");
        exit(1);
      }
      
      topdat[ntop].angpset[andx] = 1;
      andx++;
    }
    
    if(!strcmp((char *)&col1, "improper"))
    {
      printf("WARNING!!!!!!!: THIS IS NOT IMPLEMENTED.  MUST USE improperparam AND ASSIGN IMPROPER PARAMETERS IN THE TOP FILE\n");
      sscanf(line,"%d %d %d %d", &topdat[ntop].impropndx1[indx], &topdat[ntop].impropndx2[indx], &topdat[ntop].impropndx3[indx],&topdat[ntop].impropndx4[indx]);
      topdat[ntop].improppset[indx] = -1;
      indx++;
    }
    
    if(!strcmp((char *)&col1, "improperparam"))
    {
      printf("WARNING!!!!!!!: USING IMPROPER PARAMETERS FROM THE TOP FILE\n");
      if( sscanf(line,"%d %d %d %d %lf %lf", &topdat[ntop].impropndx1[indx], &topdat[ntop].impropndx2[indx],
            &topdat[ntop].impropndx3[indx],&topdat[ntop].impropndx4[indx],&topdat[ntop].impropfk[indx],&topdat[ntop].impropeq[indx]) != 6) {
        
        printf("ERROR: NOT ENOUGH ARGS FOR IMPROPERPARAM: MUST BE: NDX1 NDX2 NDX3 NDX4 FK EQ\n");
        exit(1);
      }
      topdat[ntop].improppset[indx] = 1;
      indx++;
    }

    if(!strcmp((char *)&col1, "dihedralparam"))
    {
      printf("WARNING!!!!!!!: USING DIHEDRAL PARAMETERS FROM THE TOP FILE\n");
      if( sscanf(line,"%d %d %d %d %lf %d %d %lf", &topdat[ntop].dihedndx1[dndx], &topdat[ntop].dihedndx2[dndx],
		 &topdat[ntop].dihedndx3[dndx],&topdat[ntop].dihedndx4[dndx],&topdat[ntop].dihedfk[dndx],&topdat[ntop].dihedn[dndx], &topdat[ntop].dihedeq[dndx], &topdat[ntop].dihedof[dndx]) != 8) {
        
        printf("ERROR: NOT ENOUGH ARGS FOR CHARMM DIHEDRALPARAM: MUST BE: NDX1 NDX2 NDX3 NDX4 FK n EQ ONEFOUR\n");
        exit(1);
      }
      topdat[ntop].dihedpset[dndx] = 1;
      dndx++;
    }

  }
  fclose(fpin);
}

/* count the number of things in the topology files so we can allocate */

void count_atoms(char *filename, TOPDAT *topdat, int ntop)
{
  FILE *fpin;
  int c;
  char col1[20];
  
  if((fpin = fopen(filename,"r")) == NULL)
  {
    fprintf(stderr,"ERROR: can't open infile %s\n",&filename);
    exit(1);
  }
  topdat[ntop].nat=0;
  topdat[ntop].nbnd=0;
  topdat[ntop].nang=0;
  topdat[ntop].nimprop=0;
  topdat[ntop].ndihed=0;
  
  while(c=fscanf(fpin,"%s",&col1)!= EOF)
  {
    
    if(!strcmp((char *)&col1, "atom"))
    {
      fscanf(fpin,"%*[^\n]");
      topdat[ntop].nat++;
    }

    if(!strcmp((char *)&col1, "bond") || !strcmp((char *)&col1, "bondparam"))
    {
      fscanf(fpin,"%*[^\n]");
      topdat[ntop].nbnd++;
    }

    if(!strcmp((char *)&col1, "angle") || !strcmp((char *)&col1, "angleparam"))
    {
      fscanf(fpin,"%*[^\n]");
      topdat[ntop].nang++;
    }

    if(!strcmp((char *)&col1, "improperparam") || !strcmp((char *)&col1, "improper"))
    {
      fscanf(fpin,"%*[^\n]");
      topdat[ntop].nimprop++;
    }

    if(!strcmp((char *)&col1, "dihedralparam") || !strcmp((char *)&col1, "dihdedral"))
    {
      fscanf(fpin,"%*[^\n]");
      topdat[ntop].ndihed++;
    }
    
  }
  fclose(fpin);
}
int min(const int a, const int b){
	if (a<b) return a;
	else     return b;
}
void write_psf(DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat, int ISCHARGED){
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
					atindex, topdat[i].resname[k], min(9999,molindex), topdat[i].resname[k], topdat[i].name[k],topdat[i].type[k],
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
void make_ndx(DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat, int ISCHARGED){
  FILE *fpind;
  int i,j,k,atindex,indcnt;
  int sol_flag,pro_flag,oth_flag,chr_flag;
  sol_flag = pro_flag = oth_flag = chr_flag = 0;

  if((fpind = fopen("CGindex.ndx","w")) == NULL)
    {
      fprintf(stderr,"ERROR: can't open CGindex.ndx\n");
      exit(1);
    }

  atindex = 0;
  indcnt = 0;
  for(i=0;i<sysdat->ntops;i++) {
   for(j=0;j<topdat[i].nmol;j++){
   for(k=0;k<topdat[i].nat;k++) {
    atindex++;
    if(find_solv(topdat[i].resname[k]) == 0 ){
      if (sol_flag == 0){
       fprintf(fpind,"[ SOLVENT ]\n");
       sol_flag = 1;
      }
      fprintf(fpind,"%5d ",atindex);
      indcnt++;
      if(indcnt%15 == 0)
	      fprintf(fpind,"\n");
    }
   }
   }
  }
  if (sol_flag == 1){
   fprintf(fpind,"\n\n");
  }

  atindex = 0;
  indcnt = 0;
  for(i=0;i<sysdat->ntops;i++) {
   for(j=0;j<topdat[i].nmol;j++){
   for(k=0;k<topdat[i].nat;k++) {
    atindex++;
    if(find_prot(topdat[i].resname[k]) == 0 ){
      if (pro_flag == 0){
       fprintf(fpind,"[ PROTEIN ]\n");
       pro_flag = 1;
      }
      fprintf(fpind,"%5d ",atindex);
      indcnt++;
      if(indcnt%15 == 0)
	      fprintf(fpind,"\n");
    }
   }
   }
  }
  if (pro_flag == 1){
   fprintf(fpind,"\n\n");
  }

  atindex = 0;
  indcnt = 0;
  for(i=0;i<sysdat->ntops;i++) {
   for(j=0;j<topdat[i].nmol;j++){
   for(k=0;k<topdat[i].nat;k++) {
    atindex++;
    if(find_solv(topdat[i].resname[k]) != 0 && find_prot(topdat[i].resname[k]) != 0){
      if (oth_flag == 0){
       fprintf(fpind,"[ OTHERS ]\n");
       oth_flag = 1;
      }
      fprintf(fpind,"%5d  ",atindex);
      indcnt++;
      if(indcnt%15 == 0)
	      fprintf(fpind,"\n");
    }
   }
   }
  }
  if (oth_flag == 1){
   fprintf(fpind,"\n\n");
  }

  atindex = 0;
  indcnt = 0;
  for(i=0;i<sysdat->ntops;i++) {
   for(j=0;j<topdat[i].nmol;j++){
   for(k=0;k<topdat[i].nat;k++) {
    atindex++;
    if(find_solv(topdat[i].resname[k]) != 0 && topdat[i].charge[k] != 0.){
      if (chr_flag == 0){
       fprintf(fpind,"[ CHARGED ]\n");
       chr_flag = 1;
      }
      fprintf(fpind,"%5d  ",atindex);
      indcnt++;
      if(indcnt%15 == 0)
	      fprintf(fpind,"\n");
    }
   }
   }
  }
  if (chr_flag == 1){
   fprintf(fpind,"\n\n");
  }
  fclose(fpind);
}

int find_solv(const char *s){
    const int Ntype = 3;
    const int Nbuf = 4;
    char mol_list[Ntype][Nbuf] = { "WAT", "SOD", "CLA" };
    int i;
    for (i=0;i<Ntype;i++){
      if(strncmp(s,mol_list[i],3) == 0)
	      return 0;
    }
    return 1;
}
int find_prot(const char *s){
    const int Ntype = 22;
    const int Nbuf = 4;
    char mol_list[Ntype][Nbuf] = { "GLY", "ALA", "VAL", "LEU", "ILE",
   				"SER", "THR", "MET", "CYS", "TYR",
   				"PHE", "TRP", "ASN", "GLN", "ASP",
   				"GLU", "ARG", "LYS", "PRO", "HIS",
   				"HSE", "HSD" };
    int i;
    for (i=0;i<Ntype;i++){
      if(strncmp(s,mol_list[i],3) == 0)
	      return 0;
    }
    return 1;
}

void make_top(DATABASE *database,TOPDAT *topdat,SYSDAT *sysdat, int ISCHARGED)
{
  FILE *fpout;
  int i;
  
  if((fpout = fopen("topol.top","w")) == NULL)
    {
      fprintf(stderr,"ERROR: can't open topol.top\n");
      exit(1);
    }
  fprintf(fpout,"; generated by setup_gromacs\n");
  fprintf(fpout,"; RH DEVANE\n");
  fprintf(fpout,"#include \"SPICA.itp\"\n");
  fprintf(fpout,"#include \"molecule.itp\"\n");
  fprintf(fpout,"\n");
  fprintf(fpout,"[ system ]\n");
  fprintf(fpout,"; Name\n");
  fprintf(fpout,"CG\n");
  fprintf(fpout,"\n");
  fprintf(fpout,"[ molecules ]\n");
  fprintf(fpout,"; Compound   #mols\n");

  for(i=0;i<sysdat->ntops;i++) 
    fprintf(fpout,"%s     %d\n",topdat[i].resname[0],topdat[i].nmol);

  fclose(fpout);
  
}
