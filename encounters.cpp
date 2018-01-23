/*
#########################################################
#    _ _       __                __         		#
#   (_) |     / /___ _____  ____/ /__  _____		#
#  / /| | /| / / __ `/ __ \/ __  / _ \/ ___/		#
# / / | |/ |/ / /_/ / / / / /_/ /  __/ /    		#
#/_/  |__/|__/\__,_/_/ /_/\__,_/\___/_/     		#
# Dynamics of Interestellar Wanderers			#
# Jorge I. Zuluaga et al. [)] 2017			#
# http://github.com/seap-udea/iWander.git		#
#########################################################
# Compute encounters
#########################################################
*/
#include <iwander.cpp>
using namespace std;
#define VERBOSE 0

int main(int argc,char* argv[])
{
  /*
    Example: ./encounters.exe

    Function: 

    This program perform two different tasks:

    1) Compute the LMA minimum distance and time to all stars in the
       AstroRV catalogue..

    2) Select the progenitor candidates.

    Input: 
    * wanderer.csv

    Output: 

    * encounters.csv: all the columns of the input catalog (AstroRV)
      plus:

      Cols:
          0: n
	  1-6: Position and velocity of the star for LMA purposes
	  7: Initial distance of the star, d
	  8: Minimum LMA distance, dmin
	  9: Minimum LMA time, tmin
	  10-13: Relative velocity computed with LMA vrelx,vrely,vrelz,vrel,
	  14-...: All fields in AstroRV catalog

    * candidates.csv

      Cols:
          0: n
	  1-6: Position and velocity of the star for LMA purposes
	  7: Initial distance of the star, d
	  8: Minimum LMA distance, dmin
	  9: Minimum LMA time, tmin
	  10-13: Relative velocity computed with LMA vrelx,vrely,vrelz,vrel,
	  14-...: All fields in AstroRV catalog
  */

  ////////////////////////////////////////////////////
  //CONFIGURATION
  ////////////////////////////////////////////////////
  #include <iwander.conf>
  #include <encounters.conf>

  ////////////////////////////////////////////////////
  //INITIALIZE iWANDER
  ////////////////////////////////////////////////////
  initWander();
  
  ////////////////////////////////////////////////////
  //GENERAL VARIABLES
  ////////////////////////////////////////////////////
  double tmp,telaps;
  char ctmp[100],line[10000],aline[10000],head[10000];
  double posbody[6],tbody,direction;
  char **fields=charMatrixAllocate(MAXCOLS,MAXTEXT);
  int i,j,k,n,nfields;
  double p1[3],*d1;

  int Nstars,Naccept=0,Ncand=0;

  //VARIABLES
  double tstar,dt;
  double raep,decep;
  double mura,dmura,mudec,dmudec;
  double vr,dvr;
  double ra,dra,dec,ddec;
  double postar[6];
  double gmag,gMag;
  double lep,bep,l,b;
  double par,dpar,d,dd;
  double M_Epoch_J2000[3][3];
  double TM[3][3],BM[3][3];
  double vsky[3],dvsky[3],UVW[3],dUVW[3];

  double p2[3],*d2,p1mp2[3],r1mr2[3],d1md2[3],nv[3],nv1[3],nv2[3];
  double dc1[3],dc2[3],c1[3],c2[3],drp[3];
  double dvnorm;
  double d1n2,d2n1;
  double dmin,tmin,tmin1,tmin2;
  double vrel[3],vrelmag;
  double dthres;
  double ting;

  ////////////////////////////////////////////////////
  //GLOBAL DEFINITIONS
  ////////////////////////////////////////////////////
  pxform_c("J2000","GALACTIC",0,TM);

  ////////////////////////////////////////////////////
  //READING DATA
  ////////////////////////////////////////////////////

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //READING WANDERERS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sprintf(Filename,"scratch/wanderer-%s.csv",WANDERER);
  FILE *fc=fopen(Filename,"r");

  if(fc==NULL){
    fprintf(stderr,"You must first propagate the wanderers\n");
    exit(1);
  }
  //HEADER
  fscanf(fc,"%s",line);

  //NOMINAL SOLUTION
  fscanf(fc,"%s",line);
  parseLine(line,fields,&nfields);

  for(i=Wanderer::XGAL,j=0;i<=Wanderer::VZGAL;i++) posbody[j++]=atof(fields[i]);
  tbody=atof(fields[Wanderer::TING]);
  direction=tbody/fabs(tbody);

  VPRINT(stdout,"TIME (year): %e\n",tbody/YEAR);
  VPRINT(stdout,"POS.GAL. (km): %s\n",vec2strn(posbody,3,"%.17e "));
  VPRINT(stdout,"VEL.UVW. (km/s): %s\n",vec2strn(posbody+3,3,"%.17e "));

  //SKEW LINE VECTORS
  d1=posbody+3;
  vscl_c(1E3/PARSEC,posbody,p1);
  VPRINT(stdout,"Object skew line: p=(%s), d=(%s)\n",
	 vec2strn(p1,3,"%.17e,"),
	 vec2strn(d1,3,"%.17e,"));

  ////////////////////////////////////////////////////
  //READING GAIA DATABASE
  ////////////////////////////////////////////////////
  fc=fopen("db/src/AstroRV.csv","r");

  sprintf(Filename,"scratch/encounters-%s.csv",WANDERER);
  FILE* fe=fopen(Filename,"w");

  sprintf(Filename,"scratch/candidates-%s.csv",WANDERER);
  FILE* fg=fopen(Filename,"w");

  //READING HEADER
  fscanf(fc,"%s",head);

  fprintf(fe,"n,postarx,postary,postarz,velstarx,velstary,velstarz,d,dmin,tmin,vrelx,vrely,vrelz,vrel,%s\n",head);
  fprintf(fg,"n,postarx,postary,postarz,velstarx,velstary,velstarz,d,dmin,tmin,vrelx,vrely,vrelz,vrel,%s\n",head);

  //COMPUTING LMA MINIMUM DISTANCE TO STARS
  printHeader(stdout,"COMPUTING MINIMUM DISTANCE TO CATALOG STARS");
  int Nfreq=10000;
  n=0;
  k=0;

  int qdist;
  int ndist=0;
  int nclos=0;

  sprintf(Filename,"scratch/thresholds-%s.csv",WANDERER);
  FILE* fth=fopen(Filename,"w");

  int Ndir=0;
  telaps=0.0;
  elapsedTime();
  while(fscanf(fc,"%s",line)==1){
    strcpy(aline,line);
    //fprintf(stdout,"%s\n",aline);

    if((n%Nfreq)==0){
      fprintf(stdout,"Analysing encounter of star %d...\n",n);
    }

    n++;

    //PARSE FIELDS
    parseLine(line,fields,&nfields);

    //DEBUGGING
    //if(!(strcmp(fields[Stars::TYCHO2_ID],"6995-570-1")==0)) continue;
    //if(!(strcmp(fields[Stars::HIP],"27913")==0) and VERBOSE) continue;
    //if(!(strcmp(fields[Stars::HIP],"62512")==0)) continue;
    //if(!(strcmp(fields[Stars::HIP],"64532")==0)) continue;
    //if(n<243520) continue;
    //if(n<14053 && VERBOSE) continue;
    //if(!(strcmp(fields[Stars::HIP],"43667")==0)) continue;
    //if(!(strcmp(fields[Stars::HIP],"103749")==0)) continue;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //ID
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VPRINT(stdout,"Star %d, HIP %s, TYC2 %s, HD %s, NAME %s:\n",n,
	   fields[Stars::HIP],fields[Stars::TYCHO2_ID],
	   fields[Stars::HENRYDRAPERID],fields[Stars::NAME_SIMBAD]);

    if(nfields!=47){
      fprintf(stderr,"Star %d, HIP %s, TYC2 %s, HD %s, NAME %s:, nfields = %d\n",n,
	      fields[Stars::HIP],fields[Stars::TYCHO2_ID],
	      fields[Stars::HENRYDRAPERID],fields[Stars::NAME_SIMBAD],nfields);
      exit(0);
    }
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //Primary-
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //CHECK IF ASTROMETRIC DATA IS AVAILABLE
    int qhip=0;
    ra=atof(fields[Stars::RA]);dra=atof(fields[Stars::RA_ERROR]);
    dec=atof(fields[Stars::DEC]);ddec=atof(fields[Stars::DEC_ERROR]);
    par=atof(fields[Stars::PARALLAX]);
    dpar=atof(fields[Stars::PARALLAX_ERROR]);	     
    mura=atof(fields[Stars::PMRA]);dmura=atof(fields[Stars::PMRA_ERROR]);
    mudec=atof(fields[Stars::PMDEC]);dmudec=atof(fields[Stars::PMDEC_ERROR]);
    VPRINT(stdout,"\tData for %d %s %s: %e, %e, %e, %e, %e, %e\n",
	   n,fields[Stars::HIP],fields[Stars::TYCHO2_ID],ra,dec,par,dpar,mura,mudec);

    //OTHER
    vr=atof(fields[Stars::RV]);
    dvr=atof(fields[Stars::E_RV]);
    gmag=atof(fields[Stars::PHOT_G_MEAN_MAG]);
    l=atof(fields[Stars::L]);
    b=atof(fields[Stars::B]);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //FILTER STARS
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //NULL PARALLAX
    if(par==0){
      VPRINT(stderr,"\tNull parallax for star %d %s %s. Skipping\n",
	     n,fields[Stars::HIP],fields[Stars::TYCHO2_ID]);
      continue;
    }
    if(fabs(par/dpar)<pardpar_min){
      VPRINT(stderr,"\tParallax is very uncertain for star %d %s %s. Skipping\n",
	     n,fields[Stars::HIP],fields[Stars::TYCHO2_ID]);
      continue;
    }
    if(fabs(mura/dmura)<muradmura_min){
      VPRINT(stderr,"\tProper motion in RA is very uncertain for star %d %s %s. Skipping\n",
	     n,fields[Stars::HIP],fields[Stars::TYCHO2_ID]);
      continue;
    }
    if(fabs(mudec/dmudec)<mudecdmudec_min){
      VPRINT(stderr,"\tProper motion in DEC is very uncertain for star %d %s %s. Skipping\n",
	     n,fields[Stars::HIP],fields[Stars::TYCHO2_ID]);
      continue;
    }
    if(fabs(vr/dvr)<vrdvr_min){
      VPRINT(stderr,"\tRadial velocity is very uncertain for star %d %s %s. Skipping\n",
	     n,fields[Stars::HIP],fields[Stars::TYCHO2_ID]);
      continue;
    }

    //COORDINATES AT EPOCH
    VPRINT(stdout,"\tRA(epoch) = %.17lf +/- %.3lf mas\n",ra,dra);
    VPRINT(stdout,"\tDEC(epoch) = %.17lf +/- %.3lf mas\n",dec,ddec);
    VPRINT(stdout,"\tRA(epoch) = %s, DEC(epoch) = %s\n",dec2sex(ra/15.0),dec2sex(dec));
    VPRINT(stdout,"\tParallax = %.17lf +/- %.3lf mas\n",par,dpar);
    VPRINT(stdout,"\tmu_RA(epoch) = %.17lf +/- %.3lf mas\n",mura,dmura);
    VPRINT(stdout,"\tmu_DEC(epoch) = %.17lf +/- %.3lf mas\n",mudec,dmudec);
    VPRINT(stdout,"\tv_r = %.17lf +/- %.3lf km/s\n",vr,dvr);
    VPRINT(stdout,"\tmag_g = %.3lf\n",gmag);
    VPRINT(stdout,"\tl = %.17lf, b = %.17lf\n",l,b);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //SECONDARY
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VPRINT(stdout,"\tSecondary properties:\n",gMag);
    
    //DISTANCE
    d=AU/tan(par/(60*60*1000.0)*DEG)/PARSEC;
    dd=d*dpar/par;
    VPRINT(stdout,"\td(pc) = %.17lf +/- %.3lf \n",d,dd);

    //VELOCITY WITH RESPECT TO SKY
    vsky[0]=vr;//RADIAL, km/s
    vsky[1]=KC2*mura/par;//RA, km/s
    vsky[2]=KC2*mudec/par;//DEC, km/s
    VPRINT(stdout,"\tVelocity in the sky: %s\n",vec2str(vsky,"%.5f "));

    //ABSOLUTE MAGNITUDE
    gMag=gmag-5*log10(d/10);
    VPRINT(stdout,"\tAbsolute magnitude = %.3lf\n",gMag);

    //POSITION TO STAR RESPECT TO J2000
    radrec_c(d,ra*DEG,dec*DEG,postar);
    VPRINT(stdout,"\tPosition in true epoch = %s\n",vec2str(postar,"%.17e "));

    //TRANSFORM TO GALACTIC TO CHECK
    mxv_c(TM,postar,postar);
    VPRINT(stdout,"\tPosition in Galactic = %s\n",vec2str(postar,"%.17e,"));

    //GALACTIC COORDINATES
    recrad_c(postar,&tmp,&l,&b);
    VPRINT(stdout,"\tl = %.17lf, b = %.17lf\n",l*RAD,b*RAD);

    //CALCULATE UVW
    calcUVW(ra,dec,par,dpar,mura,dmura,mudec,dmudec,vr,dvr,UVW,dUVW);
    VPRINT(stdout,"\tVelocity w.r.t. LSR: %s +/- %s\n",
	   vec2str(UVW,"%.5f,"),
	   vec2str(dUVW,"%.5f "));

    double vsun[]={USUN,VSUN+VCIRC,WSUN};
    vadd_c(UVW,vsun,vsun);
    VPRINT(stdout,"\tGalactocentric velocity: %s\n",vec2str(vsun,"%.17f,"));

    double vmag=vnorm_c(UVW);
    
    if(vmag>=vgc_max){
      VPRINT(stderr,"\tThe star %d %s %s is going too fast. Skipping\n",
	     n,fields[Stars::HIP],fields[Stars::TYCHO2_ID]);
      continue;
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //CALCULATE MINIMUM DISTANCE
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //BACKWARD STAR UNTIL TING
    vscl_c(tbody,UVW,p2);
    VPRINT(stdout,"\tDisplacement (km): %s\n",vec2str(p2,"%.17e,"));
    vscl_c(1e3/PARSEC,p2,p2);
    VPRINT(stdout,"\tDisplacement (pc): %s\n",vec2str(p2,"%.17e,"));
    vadd_c(postar,p2,p2);
    d2=UVW;
    
    VPRINT(stdout,"\tPosition particle: %s\n",vec2str(p1,"%.5lf,"));
    VPRINT(stdout,"\tVelocity particle: %s\n",vec2str(d1,"%.5lf,"));

    VPRINT(stdout,"\tPosition star today: %s\n",vec2str(postar,"%.17lf,"));
    VPRINT(stdout,"\tPosition star: %s\n",vec2str(p2,"%.17lf,"));
    VPRINT(stdout,"\tVelocity star: %s\n",vec2str(d2,"%.5lf,"));

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //CALCULATE PERISTAR
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //r1-r2
    vsub_c(p1,p2,p1mp2);
    //v1-v2
    vsub_c(d1,d2,d1md2);

    //Compute time at minimum
    vscl_c(1e3/KC1,d1md2,d1md2);
    dvnorm=vnorm_c(d1md2);
    tmin=-vdot_c(p1mp2,d1md2)/(dvnorm*dvnorm);
    VPRINT(stdout,"\tTime for minimum distance dynamic = %.17e\n",tmin);

    //Compute minimum distance
    vscl_c(tmin,d1md2,d1md2);
    vadd_c(p1mp2,d1md2,r1mr2);
    dmin=vnorm_c(r1mr2);
    VPRINT(stdout,"\tMinimum distance dynamic = %.17e\n",dmin);
    //if(dmin<10) getchar();
    
    //RELATIVE VELOCITY AT MINIMUM DISTANCE
    vsub_c(d2,d1,vrel);
    vrelmag=vnorm_c(vrel);

    //STORE INFORMATION
    fprintf(fe,"%d,",n);
    fprintf(fe,"%s%s",vec2str(postar,"%.17e,"),vec2str(UVW,"%.17e,"));
    fprintf(fe,"%.17e,%.17e,%.17e,",d,dmin,tmin);
    fprintf(fe,"%s%.17e,",vec2str(vrel,"%.17e,"),vrelmag);
    fprintf(fe,"%s",aline);
    fprintf(fe,"\n");

    //CONDITION FOR CANDIDATES
    dthres=MAX(dmax1,dfactor*d);

    VPRINT(stdout,"\tDistance threshold (tmin = %e, tRet = %e, d = %e):%e\n",tmin,tRet,d,dthres);
    if(direction*tmin>0){
      //Stars fullfiling that minimum encounter is in the direction of integration
      Ndir++;
      fprintf(fth,"%e %e %e\n",tmin,d,dthres);
      if(dmin<=dthres && vrelmag<50.0){
	fprintf(fg,"%d,",n);
	fprintf(fg,"%s%s",vec2str(postar,"%.17e,"),vec2str(UVW,"%.17e,"));
	fprintf(fg,"%.17e,%.17e,%.17e,",d,dmin,tmin);
	fprintf(fg,"%s%.17e,",vec2str(vrel,"%.17e,"),vrelmag);
	fprintf(fg,"%s",aline);
	fprintf(fg,"\n");
	Ncand++;
      }
    }
    k++;
    //if(VERBOSE) break;
    telaps+=elapsedTime();
  }
  telaps/=k;
  fprintf(stdout,"Average time per star: %f ms\n",telaps/1e-3);
  fclose(fc);
  fclose(fe);
  fclose(fg);
  fclose(fth);
  
  Nstars=n;
  Naccept=k;
  fprintf(stdout,"Total number of stars: %d\n",Nstars);
  fprintf(stdout,"Accepted stars: %d\n",Naccept);
  fprintf(stdout,"Stars with potential encounters: %d\n",Ndir);
  fprintf(stdout,"Candidate stars (fulfilling thresholds): %d\n",Ncand);

  printHeader(stdout,"DONE.",'!');

  telaps=elapsedTime(0);
  fprintf(stdout,"Total elapsed time = %.5f s (%.5f min)\n",telaps,telaps/60.0);
  return 0;
}
