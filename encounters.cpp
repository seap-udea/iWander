#include <iwander.cpp>
using namespace std;
#define VERBOSE 0

int main(int argc,char* argv[])
{
  /*
    Example: ./wanderer.exe 10 1
    Where:
      10: Number of particles
      1: Use diagonal covariance

    Function: 

    This program perform three different tasks:

    1) Calculate the time t_asymp when a single conic approximation is
       good enough to predict the future position of the interstellar
       object.

    2) Calculate the time t_ingress at which the object was at a half
       of the truncation tidal radius of the Solar System, ie. 100,000
       AU.

    3) Predict the position and velocity of the surrogate objects at
       t_ingress.

    Input: None

    Output: 

    * wanderer.csv
      Rows: 1 is for nominal solution, the rest is for random particle

      Cols:
          0:i
	  1-6:Initial elements, q,e,i,W,w,Mo,to,mu
	  7:tdb (terminal)
	  8-13:Position Ecliptic J2000
	  14-19:Position J2000
	  20-25:Position Galactic J2000
	  26:RA(h) (terminal)
	  27:DEC(deg)
	  28:l(deg)
	  29:b(deg)
	  30:d(AU)
	  31-36:Asymptotic elements, q,e,i,W,w,Mo,to,mu
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

  int Nstars;

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

  ////////////////////////////////////////////////////
  //GLOBAL DEFINITIONS
  ////////////////////////////////////////////////////
  pxform_c("J2000","GALACTIC",0,TM);

  ////////////////////////////////////////////////////
  //READING DATA
  ////////////////////////////////////////////////////

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //READING WANDERER SOLUTION
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FILE *fc=fopen("wanderer.csv","r");
  fscanf(fc,"%s",line);
  fscanf(fc,"%s",line);
  parseLine(line,fields,&nfields);

  for(i=Wanderer::XGAL,j=0;i<=Wanderer::VZGAL;i++) posbody[j++]=atof(fields[i]);
  tbody=atof(fields[Wanderer::TING]);

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

  //READING HEADER
  fscanf(fc,"%s",head);

  /*
    Columns:
    0: n
    1-6: pstar,dstar
    7-9: pos.body periastro
    10-12: pos.star periastro
    13: dmin (pc)
    14: tmin1 (yr)
    14: tmin2 (yr)
    15-17: relative velocity (km/s)
    18: relative speed (km/s)
  */
  FILE* fe=fopen("encounters.csv","w");
  FILE* fg=fopen("candidates.csv","w");

  fprintf(fe,"n,postarx,postary,postarz,velstarx,velstary,velstarz,posbodyperix,posbodyperiy,posbodyperiz,postarperix,postarperiy,postarperiz,dmin,tmin,vrelx,vrely,vrelz,vrel\n");
  fprintf(fg,"n,postarx,postary,postarz,velstarx,velstary,velstarz,posbodyperix,posbodyperiy,posbodyperiz,postarperix,postarperiy,postarperiz,dmin,tmin,vrelx,vrely,vrelz,vrel,%s\n",head);

  int Nfreq=10000;
  n=0;
  while(fscanf(fc,"%s",line)==1){
    strcpy(aline,line);

    if((n%Nfreq)==0){
      fprintf(stdout,"Analysing encounter of star %d...\n",n);
    }

    n++;

    //PARSE FIELDS
    parseLine(line,fields,&nfields);

    /*
      if(n<230611) continue;
      fprintf(stdout,"%s\n",aline);
    */
    
    /*
    for(i=0;i<nfields;i++)
      fprintf(stdout,"Field %d: %s\n",i,fields[i]);
    i=Stars::PMRA_PMDEC_CORR_TYC;
    fprintf(stdout,"%d, %s\n",i,fields[i]);
    exit(0);
    */

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //ID
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VPRINT(stdout,"Star %d, HIP %s, TYC2 %s, HD %s, NAME %s:\n",n,
	   fields[Stars::HIP],fields[Stars::TYCHO2_ID],fields[Stars::HENRYDRAPERID_TYC],fields[Stars::NAME_SIMBAD]);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //PRIMARY
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //EPOCH
    /*
    tstar=atof(fields[REF_EPOCH]);
    VPRINT(stdout,"\tEpoch (year): = %lf\n",tstar);
    sprintf(ctmp,"01/01/%d 00:00:00.000",(int)tstar);
    str2et_c(ctmp,&tstar);
    deltet_c(tstar,"et",&dt);
    tstar-=dt;
    VPRINT(stdout,"\tEpoch: tdb = %lf\n",tstar);
    */

    //CHECK IF ASTROMETRIC DATA IS AVAILABLE
    if(fields[Stars::RA][0]!='N'){
      ra=atof(fields[Stars::RA]);dra=atof(fields[Stars::RA_ERROR]);
      dec=atof(fields[Stars::DEC]);ddec=atof(fields[Stars::DEC_ERROR]);
      par=atof(fields[Stars::PARALLAX]);
      dpar=atof(fields[Stars::PARALLAX_ERROR]);	     
      mura=atof(fields[Stars::PMRA]);dmura=atof(fields[Stars::PMRA_ERROR]);
      mudec=atof(fields[Stars::PMDEC]);dmudec=atof(fields[Stars::PMDEC_ERROR]);
      fprintf(stdout,"\tGaia data for %d %s %s: %e, %e, %e, %e, %e, %e\n",
	      n,fields[Stars::HIP],fields[Stars::TYCHO2_ID],ra,dec,par,dpar,mura,mudec);
    }else if(fields[Stars::RA_HIP][0]!='N'){
      ra=atof(fields[Stars::RA_HIP]);dra=atof(fields[Stars::RA_ERROR_HIP]);
      dec=atof(fields[Stars::DEC_HIP]);ddec=atof(fields[Stars::DEC_ERROR_HIP]);
      par=atof(fields[Stars::PARALLAX_HIP]);
      dpar=atof(fields[Stars::PARALLAX_ERROR_HIP]);	     
      mura=atof(fields[Stars::PMRA_HIP]);dmura=atof(fields[Stars::PMRA_ERROR_HIP]);
      mudec=atof(fields[Stars::PMDEC_HIP]);dmudec=atof(fields[Stars::PMDEC_ERROR_HIP]);
      fprintf(stdout,"\tHipparcos data for %d %s %s: %e, %e, %e, %e, %e, %e\n",
	      n,fields[Stars::HIP],fields[Stars::TYCHO2_ID],ra,dec,par,dpar,mura,mudec);
    }else{
      ra=atof(fields[Stars::RA_TYC]);dra=atof(fields[Stars::RA_ERROR_TYC]);
      dec=atof(fields[Stars::DEC_TYC]);ddec=atof(fields[Stars::DEC_ERROR_TYC]);
      par=atof(fields[Stars::PARALLAX_TYC]);
      dpar=atof(fields[Stars::PARALLAX_ERROR_TYC]);	     
      mura=atof(fields[Stars::PMRA_TYC]);dmura=atof(fields[Stars::PMRA_ERROR_TYC]);
      mudec=atof(fields[Stars::PMDEC_TYC]);dmudec=atof(fields[Stars::PMDEC_ERROR_TYC]);
      fprintf(stdout,"\tTycho data for %d %s %s: %e, %e, %e, %e, %e, %e\n",
	      n,fields[Stars::HIP],fields[Stars::TYCHO2_ID],ra,dec,par,dpar,mura,mudec);
    }      
      
    //COORDINATES AT EPOCH
    VPRINT(stdout,"\tRA(epoch) = %.17lf +/- %.3lf mas\n",ra,dra);
    VPRINT(stdout,"\tDEC(epoch) = %.17lf +/- %.3lf mas\n",dec,ddec);
    VPRINT(stdout,"\tRA(epoch) = %s, DEC(epoch) = %s\n",dec2sex(ra/15.0),dec2sex(dec));

    //PARALLAX
    VPRINT(stdout,"\tParallax = %.17lf +/- %.3lf mas\n",par,dpar);

    //PROPER MOTION
    
    VPRINT(stdout,"\tmu_RA(epoch) = %.17lf +/- %.3lf mas\n",mura,dmura);
    VPRINT(stdout,"\tmu_DEC(epoch) = %.17lf +/- %.3lf mas\n",mudec,dmudec);

    //RADIAL VELOCITY
    vr=atof(fields[Stars::RV]);
    dvr=atof(fields[Stars::ERV]);
    VPRINT(stdout,"\tv_r = %.17lf +/- %.3lf km/s\n",vr,dvr);

    //GMAG
    gmag=atof(fields[Stars::PHOT_G_MEAN_MAG]);
    VPRINT(stdout,"\tmag_g = %.3lf\n",gmag);
    
    //READ GALACTIC COORDINATES
    l=atof(fields[Stars::L]);
    b=atof(fields[Stars::B]);
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
    vsky[1]=KC2*mura*cos(dec*DEG)/par;//RA, km/s
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
    fprintf(fe,"%s %s ",vec2str(p2,"%.5e,"),vec2str(UVW,"%.5e,"));
    fprintf(fe,"%s %s ",vec2str(c1,"%.5e,"),vec2str(c2,"%.5e,"));
    fprintf(fe,"%.5e,%.5e,",dmin,tmin);
    fprintf(fe,"%s %.5e",vec2str(vrel,"%.5e,"),vrelmag);
    fprintf(fe,"\n");

    //CONDITION FOR CANDIDATES
    if(tmin<0 && dmin<=5.0){
      fprintf(fg,"%d,",n);
      fprintf(fg,"%s %s ",vec2str(p2,"%.5e,"),vec2str(UVW,"%.5e,"));
      fprintf(fg,"%s %s ",vec2str(c1,"%.5e,"),vec2str(c2,"%.5e,"));
      fprintf(fg,"%.5e,%.5e,",dmin,tmin);
      fprintf(fg,"%s %.5e,",vec2str(vrel,"%.5e,"),vrelmag);
      fprintf(fg,"%s ",aline);
      fprintf(fg,"\n");
    }
    
    //break;
    //if(n>10) break;
  }
  /*
  fprintf(stdout,"Star %d, HIP %s, TYC2 %s, HD %s, NAME %s:\n",n,
	  fields[Stars::HIP],fields[Stars::TYCHO2_ID],fields[Stars::HENRYDRAPERID_TYC],fields[Stars::NAME_SIMBAD]);
  */
  fclose(fc);
  fclose(fe);
  fclose(fg);
  
  Nstars=n+1;
  fprintf(stdout,"Number of stars: %d\n",Nstars);
  return 0;
}

/*
hip,tycho2_id,solution_id,source_id,random_index,ref_epoch,ra,ra_error,dec,dec_error,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,ra_dec_corr,ra_parallax_corr,ra_pmra_corr,ra_pmdec_corr,dec_parallax_corr,dec_pmra_corr,dec_pmdec_corr,parallax_pmra_corr,parallax_pmdec_corr,pmra_pmdec_corr,astrometric_n_obs_al,astrometric_n_obs_ac,astrometric_n_good_obs_al,astrometric_n_good_obs_ac,astrometric_n_bad_obs_al,astrometric_n_bad_obs_ac,astrometric_delta_q,astrometric_excess_noise,astrometric_excess_noise_sig,astrometric_primary_flag,astrometric_relegation_factor,astrometric_weight_al,astrometric_weight_ac,astrometric_priors_used,matched_observations,duplicated_source,scan_direction_strength_k1,scan_direction_strength_k2,scan_direction_strength_k3,scan_direction_strength_k4,scan_direction_mean_k1,scan_direction_mean_k2,scan_direction_mean_k3,scan_direction_mean_k4,phot_g_n_obs,phot_g_mean_flux,phot_g_mean_flux_error,phot_g_mean_mag,phot_variable_flag,l,b,ecl_lon,ecl_lat,RAJ2000,DEJ2000,RV,eRV,CAT
NULL,55-72-1,1635378410781933568,16870631694208,1081909,2015.0,45.1127787196,0.206981018629,0.38084351483400003,0.150942984191,2.09081229243,0.222205844313,-1.57292770474,1.7331941971200002,-11.6616159547,0.982994152333,-0.18713556,0.41238126,-0.28435326,0.09948922,-0.13931505,0.23891407,-0.35899627,0.09374801,-0.277202,-0.9042231,79,79,79,76,0,3,NULL,0.253709333143,88.6260569708,True,2.6545267000000003,13.709007999999999,2.5737182999999998e-05,5,10,False,0.39078984,0.3336923,0.40038670000000004,0.9007094999999999,-88.302666,14.786093,-47.974358,27.022827000000003,87,967144.0942899999,601.801663362,10.5610421033,NOT_AVAILABLE,176.665031972,-48.5568508384,42.7641399122,-16.0048855323,45.112829999999995,0.38092,2.061,1.073,RAVE-DR5.tsv
*/
