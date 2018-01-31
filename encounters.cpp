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

#define VERBOSE 0 //Verbosity level
#define OSTREAM stdout //Stream where the output is redirected
#define VSTREAM stderr //Stream where the error output is redirected

int main(int argc,char* argv[])
{
  /*
    Example: ./encounters.exe

    Function: 

    This program perform two different tasks:

    1) Compute the LMA minimum distance and time to all stars in the
       AstroRV catalogue...

    2) Select the progenitor candidates.

    Input: 
    * wanderer.csv

    Output: 

    * encounters-<Wanderer>.csv: all the columns of the input catalog (AstroRV)
      plus additional information computed from the LMA approximation.

    * candidates-<Wanderer>.csv: list of objects fulfilling certain
      selection criteria that classify them as close encounters candidates.

  */

  ////////////////////////////////////////////////////
  //CONFIGURATION
  ////////////////////////////////////////////////////
  #include <iwander.conf>
  #include <encounters.conf>
  printHeader(OSTREAM,"ENCOUNTER CANDIDATES",'*');

  ////////////////////////////////////////////////////
  //INITIALIZE iWANDER
  ////////////////////////////////////////////////////
  initWander();
  
  ////////////////////////////////////////////////////
  //VARIABLES DECLARATION
  ////////////////////////////////////////////////////
  //COUNTERS
  int i,j,k,n,nfields;
  int Nfreq=10000;

  int Nstars_total=0,Nstars_cand=0;
  int Nstars_noastro=0,Nstars_null=0,Nstars_fast=0,Nstars_nothresh=0,Nstars_nodir=0;

  //COORDINATES
  double vsky[3],dvsky[3],UVW[3],dUVW[3];
  double p2[3],p1mp2[3],r1mr2[3],d1md2[3],nv[3],nv1[3],nv2[3];
  double dc1[3],dc2[3],c1[3],c2[3],drp[3];
  double *d2;

  //ASTROMETRIC POSITIONS
  double mura,dmura,mudec,dmudec;
  double vr,dvr;
  double ra,dra,dec,ddec;
  double raep,decep;
  double lep,bep,l,b;
  double par,dpar,d,dd;

  //ALLOCATION
  double M_Epoch_J2000[3][3];
  double TM[3][3],BM[3][3];
  pxform_c("J2000","GALACTIC",0,TM);
  char **fields=charMatrixAllocate(MAXCOLS,MAXTEXT);

  //MISC
  double tmp,telaps;
  int qastro;
  double posbody[6],tbody,direction;
  double p1[3],*d1;

  double tstar,dt;
  double postar[6];
  double gmag,gMag;

  double dvnorm;
  double d1n2,d2n1;
  double dmin,tmin,tmin1,tmin2;
  double vrel[3],vrelmag;
  double dthres;
  double ting;

  //FILES
  FILE *fc;
  FILE *fe;
  FILE *fg;
  FILE *fw;

  ////////////////////////////////////////////////////
  //OUTPUT FILE
  ////////////////////////////////////////////////////
  //Openning the AstroRV database
  fc=fopen("db/src/AstroRV.csv","r");
  fscanf(fc,"%s",SLINE);

  //Openning the output files
  sprintf(FILENAME,"scratch/encounters-%s.csv",Wanderer);
  fe=fopen(FILENAME,"w");
  sprintf(FILENAME,"scratch/candidates-%s.csv",Wanderer);
  fg=fopen(FILENAME,"w");

  ////////////////////////////////////////////////////
  //READING DATA
  ////////////////////////////////////////////////////
  printHeader(OSTREAM,"READING WANDERERS INFORMATION",'-');
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //READING WANDERERS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sprintf(FILENAME,"scratch/wanderer-%s.csv",Wanderer);
  fw=fopen(FILENAME,"r");
  if(fc==NULL){
    fprintf(stderr,"You must first propagate the wanderers\n");
    exit(1);
  }
  //HEADER
  fscanf(fw,"%s",LINE);

  //NOMINAL SOLUTION
  print0(OSTREAM,"\tReading nominal solution\n");
  fscanf(fw,"%s",LINE);
  parseLine(LINE,fields,&NFIELDS);

  for(i=Wanderer::XGAL,j=0;i<=Wanderer::VZGAL;i++) posbody[j++]=atof(fields[i]);
  tbody=atof(fields[Wanderer::TING]);
  direction=tbody/fabs(tbody);

  print1(VSTREAM,"\tIngress time (year): %e\n",tbody/YEAR);
  print1(VSTREAM,"\tPosition in the Galaxy (km): %s\n",
	 vec2strn(posbody,3,"%.17e "));
  print1(VSTREAM,"\tVelocity in the Galaxy (km/s): %s\n",
	 vec2strn(posbody+3,3,"%.17e "));

  //SKEW LINE VECTORS
  d1=posbody+3;
  vscl_c(1E3/PARSEC,posbody,p1);
  print1(VSTREAM,"\tObject skew line:\n\t\tp (pc)=(%s)\n\t\td(km/s)=(%s)\n",
	 vec2strn(p1,3,"%.17e "),
	 vec2strn(d1,3,"%f "));
  fclose(fw);

  ////////////////////////////////////////////////////
  //PROCEDING TO COMPUTE LMA CONDITIONS
  ////////////////////////////////////////////////////
  printHeader(OSTREAM,"COMPUTING MINIMUM DISTANCE TO CATALOG STARS",'-');

  //Header of files
  fprintf(fe,"n,postarx,postary,postarz,velstarx,velstary,velstarz,d,dmin,tmin,vrelx,vrely,vrelz,vrel,qastro,%s\n",SLINE);
  fprintf(fg,"n,postarx,postary,postarz,velstarx,velstary,velstarz,d,dmin,tmin,vrelx,vrely,vrelz,vrel,qastro,%s\n",SLINE);

  elapsedTime();
  n=0;
  k=0;
  while(fscanf(fc,"%s",LINE)==1){
    elapsedTime();

    //Save line
    strcpy(SLINE,LINE);

    //Show encounter stars
    if((n%Nfreq)==0){
      print0(OSTREAM,"\tAnalysing encounter of star %d...\n",n);
    }
    n++;

    //Parse fields
    parseLine(LINE,FIELDS,&NFIELDS);

    //Reading id
    print1(VSTREAM,"\tStar %d, HIP %s, TYC2 %s, HD %s, NAME %s:\n",n,
	   FIELDS[Stars::HIP],FIELDS[Stars::TYCHO2_ID],
	   FIELDS[Stars::HENRYDRAPERID],FIELDS[Stars::NAME_SIMBAD]);

    if(NFIELDS!=47){
      fprintf(stderr,"\t\tStar %d, HIP %s, TYC2 %s, HD %s, NAME %s:, nfields = %d\n",n,
	      FIELDS[Stars::HIP],FIELDS[Stars::TYCHO2_ID],
	      FIELDS[Stars::HENRYDRAPERID],FIELDS[Stars::NAME_SIMBAD],NFIELDS);
      exit(1);
    }
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //Primary
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //CHECK IF ASTROMETRIC DATA IS AVAILABLE
    ra=atof(FIELDS[Stars::RA]);dra=atof(FIELDS[Stars::RA_ERROR]);
    dec=atof(FIELDS[Stars::DEC]);ddec=atof(FIELDS[Stars::DEC_ERROR]);
    par=atof(FIELDS[Stars::PARALLAX]);
    dpar=atof(FIELDS[Stars::PARALLAX_ERROR]);	     
    mura=atof(FIELDS[Stars::PMRA]);dmura=atof(FIELDS[Stars::PMRA_ERROR]);
    mudec=atof(FIELDS[Stars::PMDEC]);dmudec=atof(FIELDS[Stars::PMDEC_ERROR]);
    //OTHER
    vr=atof(FIELDS[Stars::RV]);
    dvr=atof(FIELDS[Stars::E_RV]);
    gmag=atof(FIELDS[Stars::PHOT_G_MEAN_MAG]);
    l=atof(FIELDS[Stars::L]);
    b=atof(FIELDS[Stars::B]);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //FILTER STARS
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(par<=0){
      print2(VSTREAM,"\tNull parallax for star %d %s %s. Skipping\n",
	     n,FIELDS[Stars::HIP],FIELDS[Stars::TYCHO2_ID]);
      Nstars_null++;
      continue;
    }
    if(fabs(par/dpar)<pardpar_min){
      print2(VSTREAM,"\tParallax is very uncertain for star %d %s %s. Skipping\n",
	     n,FIELDS[Stars::HIP],FIELDS[Stars::TYCHO2_ID]);
      Nstars_noastro++;
      continue;
    }
    if(fabs(mura/dmura)<muradmura_min){
      print2(VSTREAM,"\tProper motion in RA is very uncertain for star %d %s %s. Skipping\n",
	     n,FIELDS[Stars::HIP],FIELDS[Stars::TYCHO2_ID]);
      Nstars_noastro++;
      continue;
    }
    if(fabs(mudec/dmudec)<mudecdmudec_min){
      print2(VSTREAM,"\tProper motion in DEC is very uncertain for star %d %s %s. Skipping\n",
	     n,FIELDS[Stars::HIP],FIELDS[Stars::TYCHO2_ID]);
      Nstars_noastro++;
      continue;
    }
    if(fabs(vr/dvr)<vrdvr_min){
      print2(VSTREAM,"\tRadial velocity is very uncertain for star %d %s %s. Skipping\n",
	     n,FIELDS[Stars::HIP],FIELDS[Stars::TYCHO2_ID]);
      Nstars_noastro++;
      continue;
    }
    
    //Astromerty quality factor
    qastro=(int)(MIN(fabs(par/dpar),MIN(fabs(mura/dmura),MIN(fabs(mudec/dmudec),fabs(vr/dvr)))));

    //COORDINATES AT EPOCH
    print1(VSTREAM,"\t\tPrimary properties:\n",gMag);
    print1(VSTREAM,"\t\t\tRA(epoch) = %.17lf +/- %.3lf mas\n",ra,dra);
    print1(VSTREAM,"\t\t\tDEC(epoch) = %.17lf +/- %.3lf mas\n",dec,ddec);
    print1(VSTREAM,"\t\t\tRA(epoch) = %s, DEC(epoch) = %s\n",dec2sex(ra/15.0),dec2sex(dec));
    print1(VSTREAM,"\t\t\tParallax = %.17lf +/- %.3lf mas\n",par,dpar);
    print1(VSTREAM,"\t\t\tmu_RA(epoch) = %.17lf +/- %.3lf mas\n",mura,dmura);
    print1(VSTREAM,"\t\t\tmu_DEC(epoch) = %.17lf +/- %.3lf mas\n",mudec,dmudec);
    print1(VSTREAM,"\t\t\tv_r = %.17lf +/- %.3lf km/s\n",vr,dvr);
    print1(VSTREAM,"\t\t\tmag_g = %.3lf\n",gmag);
    print1(VSTREAM,"\t\t\tl = %.17lf, b = %.17lf\n",l,b);
    print1(VSTREAM,"\t\t\tqastro = %.3lf\n",qastro);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //SECONDARY
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print1(VSTREAM,"\t\tSecondary properties:\n",gMag);
    
    //DISTANCE
    d=AU/tan(par/(60*60*1000.0)*DEG)/PARSEC;
    dd=d*dpar/par;
    print1(VSTREAM,"\t\t\td(pc) = %.17lf +/- %.3lf \n",d,dd);

    //VELOCITY WITH RESPECT TO SKY
    vsky[0]=vr;//RADIAL, km/s
    vsky[1]=KC2*mura/par;//RA, km/s
    vsky[2]=KC2*mudec/par;//DEC, km/s
    print1(VSTREAM,"\t\t\tVelocity in the sky: %s\n",vec2str(vsky,"%.5f "));

    //ABSOLUTE MAGNITUDE
    gMag=gmag-5*log10(d/10);
    print1(VSTREAM,"\t\t\tAbsolute magnitude = %.3lf\n",gMag);

    //POSITION TO STAR RESPECT TO J2000
    radrec_c(d,ra*DEG,dec*DEG,postar);
    print1(VSTREAM,"\t\t\tPosition in true epoch = %s\n",vec2str(postar,"%.17e "));

    //TRANSFORM TO GALACTIC TO CHECK
    mxv_c(TM,postar,postar);
    print1(VSTREAM,"\t\t\tPosition in Galactic = %s\n",vec2str(postar,"%.17e,"));

    //GALACTIC COORDINATES
    recrad_c(postar,&tmp,&l,&b);
    print1(VSTREAM,"\t\t\tl = %.17lf, b = %.17lf\n",l*RAD,b*RAD);

    //CALCULATE UVW
    calcUVW(ra,dec,par,dpar,mura,dmura,mudec,dmudec,vr,dvr,UVW,dUVW);
    print1(VSTREAM,"\t\t\tVelocity w.r.t. LSR: %s +/- %s\n",
	   vec2str(UVW,"%.5f,"),
	   vec2str(dUVW,"%.5f "));

    double vsun[]={USUN,VSUN+VCIRC,WSUN};
    vadd_c(UVW,vsun,vsun);
    print1(VSTREAM,"\t\t\tGalactocentric velocity: %s\n",vec2str(vsun,"%.17f,"));

    double vmag=vnorm_c(UVW);
    
    if(vmag>=vgc_max){
      print2(VSTREAM,"\t\t***The star %d %s %s is going too fast. Skipping***\n",
	     n,FIELDS[Stars::HIP],FIELDS[Stars::TYCHO2_ID]);
      Nstars_fast++;
      continue;
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //CALCULATE MINIMUM DISTANCE
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //BACKWARD STAR UNTIL TING
    vscl_c(tbody,UVW,p2);
    print1(VSTREAM,"\t\tDisplacement (km): %s\n",vec2str(p2,"%.17e,"));
    vscl_c(1e3/PARSEC,p2,p2);
    print1(VSTREAM,"\t\tDisplacement (pc): %s\n",vec2str(p2,"%.17e,"));
    vadd_c(postar,p2,p2);
    d2=UVW;
    
    print1(VSTREAM,"\t\tPosition particle: %s\n",vec2str(p1,"%.5lf,"));
    print1(VSTREAM,"\t\tVelocity particle: %s\n",vec2str(d1,"%.5lf,"));

    print1(VSTREAM,"\t\tPosition star today: %s\n",vec2str(postar,"%.17lf,"));
    print1(VSTREAM,"\t\tPosition star: %s\n",vec2str(p2,"%.17lf,"));
    print1(VSTREAM,"\t\tVelocity star: %s\n",vec2str(d2,"%.5lf,"));

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
    print1(VSTREAM,"\t\tTime for minimum distance dynamic = %.17e\n",tmin);

    //Compute minimum distance
    vscl_c(tmin,d1md2,d1md2);
    vadd_c(p1mp2,d1md2,r1mr2);
    dmin=vnorm_c(r1mr2);
    print1(VSTREAM,"\t\tMinimum distance dynamic = %.17e\n",dmin);
    
    //RELATIVE VELOCITY AT MINIMUM DISTANCE
    vsub_c(d2,d1,vrel);
    vrelmag=vnorm_c(vrel);
    print1(VSTREAM,"\t\tRelative velocity at minimum = %.17e\n",vrelmag);

    //STORE INFORMATION
    fprintf(fe,"%d,",n);
    fprintf(fe,"%s%s",vec2str(postar,"%.17e,"),vec2str(UVW,"%.17e,"));
    fprintf(fe,"%.17e,%.17e,%.17e,",d,dmin,tmin);
    fprintf(fe,"%s%.17e,",vec2str(vrel,"%.17e,"),vrelmag);
    fprintf(fe,"%d,",qastro);
    fprintf(fe,"%s",SLINE);
    fprintf(fe,"\n");

    //CONDITION FOR CANDIDATES
    dthres=MAX(dmax1,dfactor*d);

    print1(VSTREAM,"\t\tDistance threshold (tmin = %e, TRet = %e, d = %e):%e\n",tmin,TRet,d,dthres);

    if(direction*tmin>0){
      //Stars fullfiling that minimum encounter is in the direction of integration
      print2(VSTREAM,"\t\tThe star is in the right direction\n");
      if(dmin<=dthres && vrelmag<50.0){
	print2(VSTREAM,"\t\t***The star is accepted***\n");
	fprintf(fg,"%d,",n);
	fprintf(fg,"%s%s",vec2str(postar,"%.17e,"),vec2str(UVW,"%.17e,"));
	fprintf(fg,"%.17e,%.17e,%.17e,",d,dmin,tmin);
	fprintf(fg,"%s%.17e,",vec2str(vrel,"%.17e,"),vrelmag);
	fprintf(fg,"%d,",qastro);
	fprintf(fg,"%s",SLINE);
	fprintf(fg,"\n");
	Nstars_cand++;
      }else{
	print2(VSTREAM,"\t\t***The star is beyond the threshold. Skipping***\n");
	Nstars_nothresh++;
      }
    }else{
      print2(VSTREAM,"\t\t***The star is not in the right direction. Skipping***\n");
      Nstars_nodir++;
      continue;
    }
    k++;
    TELAPS+=elapsedTime();
  }
  printHeader(OSTREAM,"RESULTS SUMMARY",'-');

  TELAPS/=k;
  print0(OSTREAM,"\tAverage time per star: %f ms\n",TELAPS/1e-3);
  fclose(fc);
  fclose(fe);
  fclose(fg);
  
  //Summary
  print0(OSTREAM,"\tTotal number of stars: %d\n",Nstars_total);
  print0(OSTREAM,"\tStars rejected by astrometry quality: %d\n",Nstars_noastro);
  print0(OSTREAM,"\tStars rejected by null quantity: %d\n",Nstars_null);
  print0(OSTREAM,"\tStars too fast: %d\n",Nstars_fast);
  print0(OSTREAM,"\tStars beyond the distance threshold: %d\n",Nstars_nothresh);
  print0(OSTREAM,"\tStars in the wrong direction: %d\n",Nstars_nodir);
  print0(OSTREAM,"\tAccepted stars: %d\n",Nstars_cand);

  TELAPS=elapsedTime(0);
  print0(OSTREAM,"Total elapsed time = %.5f s (%.5f min)\n",TELAPS,TELAPS/60.0);
  return 0;
}
