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

    Input: None

    Output: 

    * encounters.csv
      Cols:
          0:i
	  1-8:Initial elements, q,e,i,W,w,Mo,to,mu
	  9-15:Asymptotic elements, q,e,i,W,w,Mo,to,mu
	  16:ting (time of ingress, seconds)
	  17-22:Position Ecliptic J2000
	  23-28:Position J2000
	  29-34:Position Galactic J2000
	  35:RA(h) (ingress)
	  36:DEC(deg)
	  37:l(deg)
	  38:b(deg)

    * candidates.csv
      Cols:
          0: n
          1-6: postar (present)
	  7-9: pos.body periastro
	  10-12: pos.star periastro
	  13: dmin (pc)
	  14: tmin1 (yr)
	  14: tmin2 (yr)
	  15-17: relative velocity (km/s)
	  18: relative speed (km/s)
  */

  ////////////////////////////////////////////////////
  //CONFIGURATION
  ////////////////////////////////////////////////////
  #include <encounters.conf>

  ////////////////////////////////////////////////////
  //INITIALIZE iWANDER
  ////////////////////////////////////////////////////
  initWander();
  
  ////////////////////////////////////////////////////
  //GENERAL VARIABLES
  ////////////////////////////////////////////////////
  double tmp;
  char ctmp[100],line[10000],aline[10000],head[10000];
  double posbody[6],tbody;
  char **fields=charMatrixAllocate(MAXCOLS,MAXTEXT);
  int i,j,k,n,nfields;
  double p1[3],*d1;

  int Nstars,Naccept;

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

  double *p2,*d2,p1mp2[3],r1mr2[3],d1md2[3],nv[3],nv1[3],nv2[3];
  double dc1[3],dc2[3],c1[3],c2[3],drp[3];
  double dvnorm;
  double d1n2,d2n1;
  double dmin,tmin,tmin1,tmin2;
  double vrel[3],vrelmag;

  double dthres;

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
  FILE *fc=fopen("wanderer.csv","r");
  //HEADER
  fscanf(fc,"%s",line);

  //NOMINAL SOLUTION
  fscanf(fc,"%s",line);
  parseLine(line,fields,&nfields);

  for(i=Wanderer::XGAL,j=0;i<=Wanderer::VZGAL;i++) posbody[j++]=atof(fields[i]);
  tbody=atof(fields[Wanderer::TING]);

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
  FILE* fe=fopen("encounters.csv","w");
  FILE* fg=fopen("candidates.csv","w");

  //READING HEADER
  fscanf(fc,"%s",head);

  fprintf(fe,"n,postarx,postary,postarz,velstarx,velstary,velstarz,posbodyperix,posbodyperiy,posbodyperiz,postarperix,postarperiy,postarperiz,d,dmin,tmin,vrelx,vrely,vrelz,vrel\n");
  fprintf(fg,"n,postarx,postary,postarz,velstarx,velstary,velstarz,posbodyperix,posbodyperiy,posbodyperiz,postarperix,postarperiy,postarperiz,d,dmin,tmin,vrelx,vrely,vrelz,vrel,%s\n",head);

  //COMPUTING LMA MINIMUM DISTANCE TO STARS
  int Nfreq=10000;
  n=0;
  k=0;
  while(fscanf(fc,"%s",line)==1){
    strcpy(aline,line);

    if((n%Nfreq)==0){
      fprintf(stdout,"Analysing encounter of star %d...\n",n);
    }

    n++;

    //PARSE FIELDS
    parseLine(line,fields,&nfields);
    //if(!(strcmp(fields[Stars::HIP],"62512")==0)) continue;
    //if(!(strcmp(fields[Stars::HIP],"64532")==0)) continue;
    //if(n<243520) continue;
    if(n<14053 && VERBOSE) continue;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //ID
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VPRINT(stdout,"Star %d, HIP %s, TYC2 %s, HD %s, NAME %s:\n",n,
	   fields[Stars::HIP],fields[Stars::TYCHO2_ID],
	   fields[Stars::HENRYDRAPERID_TYC],fields[Stars::NAME_SIMBAD]);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //PRIMARY
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //CHECK IF ASTROMETRIC DATA IS AVAILABLE
    if(fields[Stars::RA][0]!='N'){
      ra=atof(fields[Stars::RA]);dra=atof(fields[Stars::RA_ERROR]);
      dec=atof(fields[Stars::DEC]);ddec=atof(fields[Stars::DEC_ERROR]);
      par=atof(fields[Stars::PARALLAX]);
      dpar=atof(fields[Stars::PARALLAX_ERROR]);	     
      mura=atof(fields[Stars::PMRA]);dmura=atof(fields[Stars::PMRA_ERROR]);
      mudec=atof(fields[Stars::PMDEC]);dmudec=atof(fields[Stars::PMDEC_ERROR]);
      VPRINT(stdout,"\tGaia data for %d %s %s: %e, %e, %e, %e, %e, %e\n",
	      n,fields[Stars::HIP],fields[Stars::TYCHO2_ID],ra,dec,par,dpar,mura,mudec);
    }else if(fields[Stars::RA_HIP][0]!='N'){
      ra=atof(fields[Stars::RA_HIP2]);dra=atof(fields[Stars::RA_ERROR_HIP]);
      dec=atof(fields[Stars::DEC_HIP2]);ddec=atof(fields[Stars::DEC_ERROR_HIP]);
      par=atof(fields[Stars::PARALLAX_HIP2]);
      dpar=atof(fields[Stars::PARALLAX_ERROR_HIP2]);	     
      mura=atof(fields[Stars::PMRA_HIP2]);dmura=atof(fields[Stars::PMRA_ERROR_HIP2]);
      mudec=atof(fields[Stars::PMDEC_HIP2]);dmudec=atof(fields[Stars::PMDEC_ERROR_HIP2]);
      VPRINT(stdout,"\tHipparcos data for %d %s %s: %e, %e, %e, %e, %e, %e\n",
	      n,fields[Stars::HIP],fields[Stars::TYCHO2_ID],ra,dec,par,dpar,mura,mudec);
    }else{
      fprintf(stderr,"\tAstrometric data not available for star %d %s %s\n",
	      n,fields[Stars::HIP],fields[Stars::TYCHO2_ID]);
      continue;
    }      
    //OTHER
    vr=atof(fields[Stars::RV]);
    dvr=atof(fields[Stars::ERV]);
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

    double vmag=vnorm_c(UVW);
    
    if(vmag>=vgc_max){
      fprintf(stderr,"\tThe star %d %s %s is going too fast. Skipping\n",
	      n,fields[Stars::HIP],fields[Stars::TYCHO2_ID]);
      continue;
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //CALCULATE MINIMUM DISTANCE
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //MINIMUM DISTANCE
    p2=postar;
    d2=UVW;

    VPRINT(stdout,"\tPosition particle: %s\n",vec2str(p1,"%.5lf,"));
    VPRINT(stdout,"\tVelocity particle: %s\n",vec2str(d1,"%.5lf,"));
    VPRINT(stdout,"\tPosition star: %s\n",vec2str(p2,"%.5lf,"));
    VPRINT(stdout,"\tVelocity star: %s\n",vec2str(d2,"%.5lf,"));

    vsub_c(p1,p2,p1mp2);
    ucrss_c(d1,d2,nv);
    VPRINT(stdout,"\tNormal vector to skew lines: %s\n",vec2str(nv,"%.5f "));
    dmin=fabs(vdot_c(nv,p1mp2));
    VPRINT(stdout,"\tMinimum distance = %.17e\n",dmin);

    //POINT OF MINIMUM DISTANCE
    ucrss_c(d1,nv,nv1);
    ucrss_c(d2,nv,nv2);
    d1n2=vdot_c(d1,nv2);
    d2n1=vdot_c(d2,nv1);
    vscl_c(-vdot_c(p1mp2,nv2)/d1n2,d1,dc1);
    vadd_c(p1,dc1,c1);
    vscl_c(+vdot_c(p1mp2,nv1)/d2n1,d2,dc2);
    vadd_c(p2,dc2,c2);
    VPRINT(stdout,"\tPoint of peristar (body): %s\n",vec2str(c1,"%.5f "));
    VPRINT(stdout,"\tPoint of peristar (star): %s\n",vec2str(c2,"%.5f "));
    vsub_c(c1,p1,drp);
    VPRINT(stdout,"\tDistance traveleded by body until crossing point: %s\n",vec2str(drp,"%.5f "));

    //TIME OF MINIMUM DISTANCE
    tmin1=(c1[0]-p1[0])/(d1[0]*1e3/KC1);
    tmin2=(c2[0]-p2[0])/(d2[0]*1e3/KC1);
    VPRINT(stdout,"\tTime for minimum distance 1 = %.17e\n",tmin1);
    VPRINT(stdout,"\tTime for minimum distance 2 = %.17e\n",tmin2);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //CALCULATE PERISTAR
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vsub_c(d1,d2,d1md2);
    vscl_c(1e3/KC1,d1md2,d1md2);
    dvnorm=vnorm_c(d1md2);
    tmin=-vdot_c(p1mp2,d1md2)/(dvnorm*dvnorm);
    VPRINT(stdout,"\tTime for minimum distance dynamic = %.17e\n",tmin);
    vscl_c(tmin,d1md2,d1md2);
    vadd_c(p1mp2,d1md2,r1mr2);
    dmin=vnorm_c(r1mr2);
    VPRINT(stdout,"\tMinimum distance dynamic = %.17e\n",dmin);
    
    //RELATIVE VELOCITY AT MINIMUM DISTANCE
    vsub_c(d2,d1,vrel);
    vrelmag=vnorm_c(vrel);

    //STORE INFORMATION
    fprintf(fe,"%d,",n);
    fprintf(fe,"%s%s",vec2str(p2,"%.5e,"),vec2str(UVW,"%.5e,"));
    fprintf(fe,"%s%s",vec2str(c1,"%.5e,"),vec2str(c2,"%.5e,"));
    fprintf(fe,"%.5e,%.5e,%.5e,",d,dmin,tmin);
    fprintf(fe,"%s%.5e",vec2str(vrel,"%.5e,"),vrelmag);
    fprintf(fe,"\n");

    //CONDITION FOR CANDIDATES
    dthres=d<10*dmax1?0.1*d:dmax1;
    VPRINT(stdout,"\tDistance threshold (d = %e):%e\n",d,dthres);
    if(tmin<0 && dmin<=dthres && vrelmag<50.0){
      fprintf(fg,"%d,",n);
      fprintf(fg,"%s%s",vec2str(p2,"%.5e,"),vec2str(UVW,"%.5e,"));
      fprintf(fg,"%s%s",vec2str(c1,"%.5e,"),vec2str(c2,"%.5e,"));
      fprintf(fg,"%.5e,%.5e,%.5e,",d,dmin,tmin);
      fprintf(fg,"%s%.5e,",vec2str(vrel,"%.5e,"),vrelmag);
      fprintf(fg,"%s",aline);
      fprintf(fg,"\n");
    }
    k++;
    if(VERBOSE) break;
  }
  fclose(fc);
  fclose(fe);
  fclose(fg);
  
  Nstars=n+1;
  Naccept=k;
  fprintf(stdout,"Total number of stars: %d\n",Nstars);
  fprintf(stdout,"Accepted stars: %d\n",Naccept);
  return 0;
}