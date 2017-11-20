#include <gravray.cpp>
using namespace std;

#define VERBOSE 1

enum {HIP,TYCHO2_ID,SOLUTION_ID,SOURCE_ID,RANDOM_INDEX,REF_EPOCH,RA,RA_ERROR,DEC,DEC_ERROR,PARALLAX,PARALLAX_ERROR,PMRA,PMRA_ERROR,PMDEC,PMDEC_ERROR,RA_DEC_CORR,RA_PARALLAX_CORR,RA_PMRA_CORR,RA_PMDEC_CORR,DEC_PARALLAX_CORR,DEC_PMRA_CORR,DEC_PMDEC_CORR,PARALLAX_PMRA_CORR,PARALLAX_PMDEC_CORR,PMRA_PMDEC_CORR,ASTROMETRIC_N_OBS_AL,ASTROMETRIC_N_OBS_AC,ASTROMETRIC_N_GOOD_OBS_AL,ASTROMETRIC_N_GOOD_OBS_AC,ASTROMETRIC_N_BAD_OBS_AL,ASTROMETRIC_N_BAD_OBS_AC,ASTROMETRIC_DELTA_Q,ASTROMETRIC_EXCESS_NOISE,ASTROMETRIC_EXCESS_NOISE_SIG,ASTROMETRIC_PRIMARY_FLAG,ASTROMETRIC_RELEGATION_FACTOR,ASTROMETRIC_WEIGHT_AL,ASTROMETRIC_WEIGHT_AC,ASTROMETRIC_PRIORS_USED,MATCHED_OBSERVATIONS,DUPLICATED_SOURCE,SCAN_DIRECTION_STRENGTH_K1,SCAN_DIRECTION_STRENGTH_K2,SCAN_DIRECTION_STRENGTH_K3,SCAN_DIRECTION_STRENGTH_K4,SCAN_DIRECTION_MEAN_K1,SCAN_DIRECTION_MEAN_K2,SCAN_DIRECTION_MEAN_K3,SCAN_DIRECTION_MEAN_K4,PHOT_G_N_OBS,PHOT_G_MEAN_FLUX,PHOT_G_MEAN_FLUX_ERROR,PHOT_G_MEAN_MAG,PHOT_VARIABLE_FLAG,L,B,ECL_LON,ECL_LAT,RAJ2000,DEJ2000,RV,ERV,CAT};

int main(int argc,char* argv[])
{
  /*
    Example: ./select.exe 

    Input: 
    * cloud.data
    * RVGaia.data
    Output: 
  */
  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  int Npart=1;
  if(argc>1){
    Npart=atoi(argv[1]);
  }

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initSpice();
  
  ////////////////////////////////////////////////////
  //GENERAL VARIABLES
  ////////////////////////////////////////////////////
  double tmp;
  char ctmp[100],line[10000],aline[10000];

  ////////////////////////////////////////////////////
  //READING DATA
  ////////////////////////////////////////////////////

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //READING NOMINAL SOLUTION
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FILE *fc=fopen("cloud.data","r");
  double posbody[6],tbody;
  for(int i=0;i<=0;i++) fscanf(fc,"%lf",&tmp);
  fscanf(fc,"%lf",&tbody);
  VPRINT(stdout,"Epoch body: %.17e\n",tbody);
  for(int i=2;i<=39;i++) fscanf(fc,"%lf",&tmp);
  for(int j=0;j<6;j++) fscanf(fc,"%lf",&posbody[j]);//GALACTIC
  VPRINT(stdout,"UVW: %s\n",vec2strn(posbody+3,3,"%.17e "));
  fclose(fc);

  //SKEW LINE VECTORS
  double p1[3],*d1;
  d1=posbody+3;
  vscl_c(1E3/PARSEC,posbody,p1);
  VPRINT(stdout,"Object skew line: p=(%s), d=(%s)\n",
	 vec2strn(p1,3,"%.17e,"),
	 vec2strn(d1,3,"%.17e,"));

  ////////////////////////////////////////////////////
  //READING GAIA DATABASE
  ////////////////////////////////////////////////////
  fc=fopen("../RVGaia/DB/RVGaia.csv","r");

  //READING HEADER
  fscanf(fc,"%s",line);

  //LOOP OVER STARS
  int n=0;
  int Nstars;
  int nfields;

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

  //GALACTIC TRANSFORMATION MATRIX
  pxform_c("J2000","GALACTIC",0,TM);

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

  fprintf(fg,"n,postarx,postary,postarz,velstarx,velstary,velstarz,posbodyperix,posbodyperiy,posbodyperiz,postarperix,postarperiy,postarperiz,dmin,tmin,vrelx,vrely,vrelz,vrel,hip,tycho2_id,solution_id,source_id,random_index,ref_epoch,ra,ra_error,dec,dec_error,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,ra_dec_corr,ra_parallax_corr,ra_pmra_corr,ra_pmdec_corr,dec_parallax_corr,dec_pmra_corr,dec_pmdec_corr,parallax_pmra_corr,parallax_pmdec_corr,pmra_pmdec_corr,astrometric_n_obs_al,astrometric_n_obs_ac,astrometric_n_good_obs_al,astrometric_n_good_obs_ac,astrometric_n_bad_obs_al,astrometric_n_bad_obs_ac,astrometric_delta_q,astrometric_excess_noise,astrometric_excess_noise_sig,astrometric_primary_flag,astrometric_relegation_factor,astrometric_weight_al,astrometric_weight_ac,astrometric_priors_used,matched_observations,duplicated_source,scan_direction_strength_k1,scan_direction_strength_k2,scan_direction_strength_k3,scan_direction_strength_k4,scan_direction_mean_k1,scan_direction_mean_k2,scan_direction_mean_k3,scan_direction_mean_k4,phot_g_n_obs,phot_g_mean_flux,phot_g_mean_flux_error,phot_g_mean_mag,phot_variable_flag,l,b,ecl_lon,ecl_lat,RAJ2000,DEJ2000,RV,eRV,CAT\n");

  char **fields=(char**)malloc(MAXCOLS*sizeof(char*));
  for(int i=0;i<MAXCOLS;i++) fields[i]=(char*)malloc(MAXTEXT*sizeof(char));

  int Nfreq=10000;
  while(fscanf(fc,"%s",line)==1){

    strcpy(aline,line);

    //SHOW
    if((n%Nfreq)==0){
      fprintf(stdout,"Analysing encounter of star %d...\n",n);
    }

    //PARSE FIELDS
    parseLine(line,fields,&nfields);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //ID
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VPRINT(stdout,"Star %d, HIP %s, TYC2 %s:\n",n,fields[HIP],fields[TYCHO2_ID]);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //PRIMARY
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //EPOCH
    tstar=atof(fields[REF_EPOCH]);
    VPRINT(stdout,"\tEpoch (year): = %lf\n",tstar);
    sprintf(ctmp,"01/01/%d 00:00:00.000",(int)tstar);
    str2et_c(ctmp,&tstar);
    deltet_c(tstar,"et",&dt);
    tstar-=dt;
    VPRINT(stdout,"\tEpoch: tdb = %lf\n",tstar);

    //COORDINATES AT EPOCH
    ra=atof(fields[RA]);dra=atof(fields[RA_ERROR]);
    dec=atof(fields[DEC]);ddec=atof(fields[DEC_ERROR]);
    VPRINT(stdout,"\tRA(epoch) = %.17lf +/- %.3lf mas\n",ra,dra);
    VPRINT(stdout,"\tDEC(epoch) = %.17lf +/- %.3lf mas\n",dec,ddec);
    VPRINT(stdout,"\tRA(epoch) = %s, DEC(epoch) = %s\n",dec2sex(ra/15.0),dec2sex(dec));

    //PARALLAX
    par=atof(fields[PARALLAX]);
    dpar=atof(fields[PARALLAX_ERROR]);	     
    VPRINT(stdout,"\tParallax = %.17lf +/- %.3lf mas\n",par,dpar);

    //PROPER MOTION
    mura=atof(fields[PMRA]);dmura=atof(fields[PMRA_ERROR]);
    mudec=atof(fields[PMDEC]);dmudec=atof(fields[PMDEC_ERROR]);
    
    VPRINT(stdout,"\tmu_RA(epoch) = %.17lf +/- %.3lf mas\n",mura,dmura);
    VPRINT(stdout,"\tmu_DEC(epoch) = %.17lf +/- %.3lf mas\n",mudec,dmudec);

    //RADIAL VELOCITY
    vr=atof(fields[RV]);
    dvr=atof(fields[ERV]);
    VPRINT(stdout,"\tv_r = %.17lf +/- %.3lf km/s\n",vr,dvr);

    //GMAG
    gmag=atof(fields[PHOT_G_MEAN_MAG]);
    VPRINT(stdout,"\tmag_g = %.3lf\n",gmag);
    
    //READ GALACTIC COORDINATES
    l=atof(fields[L]);
    b=atof(fields[B]);
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
    exit(0);

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
    /*
    if(strcmp(fields[HIP],"40170.0")==0 || n==214593){
      printf("Este es: HIP=%s, tmin=%lf, dmin=%lf\n",fields[HIP],tmin,dmin);
      printf("Este es: p1=[%s], d1=[%s], p2=[%s], d2=[%s]\n",
	     vec2str(p1,"%.3lf,"),
	     vec2str(d1,"%.3lf,"),
	     vec2str(p2,"%.3lf,"),
	     vec2str(d2,"%.3lf,"));
      printf("Este es: p1-p2=[%s], d1-d2=[%s], dfut=[%s], r1-r2=[%s]\n",
	     vec2str(p1mp2,"%.3lf,"),
	     vec2str(d1md2,"%.3lf,"),
	     vec2str(dfut,"%.3lf,"),
	     vec2str(r1mr2,"%.3lf,"));
      exit(0);
    }
    */
    
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
    
    //if(n>10) break;
    n++;
  }
  fclose(fc);
  fclose(fe);
  fclose(fg);
  
  Nstars=n+1;
  VPRINT(stdout,"Number of stars: %d\n",Nstars);
  return 0;
}

/*
hip,tycho2_id,solution_id,source_id,random_index,ref_epoch,ra,ra_error,dec,dec_error,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,ra_dec_corr,ra_parallax_corr,ra_pmra_corr,ra_pmdec_corr,dec_parallax_corr,dec_pmra_corr,dec_pmdec_corr,parallax_pmra_corr,parallax_pmdec_corr,pmra_pmdec_corr,astrometric_n_obs_al,astrometric_n_obs_ac,astrometric_n_good_obs_al,astrometric_n_good_obs_ac,astrometric_n_bad_obs_al,astrometric_n_bad_obs_ac,astrometric_delta_q,astrometric_excess_noise,astrometric_excess_noise_sig,astrometric_primary_flag,astrometric_relegation_factor,astrometric_weight_al,astrometric_weight_ac,astrometric_priors_used,matched_observations,duplicated_source,scan_direction_strength_k1,scan_direction_strength_k2,scan_direction_strength_k3,scan_direction_strength_k4,scan_direction_mean_k1,scan_direction_mean_k2,scan_direction_mean_k3,scan_direction_mean_k4,phot_g_n_obs,phot_g_mean_flux,phot_g_mean_flux_error,phot_g_mean_mag,phot_variable_flag,l,b,ecl_lon,ecl_lat,RAJ2000,DEJ2000,RV,eRV,CAT
NULL,55-72-1,1635378410781933568,16870631694208,1081909,2015.0,45.1127787196,0.206981018629,0.38084351483400003,0.150942984191,2.09081229243,0.222205844313,-1.57292770474,1.7331941971200002,-11.6616159547,0.982994152333,-0.18713556,0.41238126,-0.28435326,0.09948922,-0.13931505,0.23891407,-0.35899627,0.09374801,-0.277202,-0.9042231,79,79,79,76,0,3,NULL,0.253709333143,88.6260569708,True,2.6545267000000003,13.709007999999999,2.5737182999999998e-05,5,10,False,0.39078984,0.3336923,0.40038670000000004,0.9007094999999999,-88.302666,14.786093,-47.974358,27.022827000000003,87,967144.0942899999,601.801663362,10.5610421033,NOT_AVAILABLE,176.665031972,-48.5568508384,42.7641399122,-16.0048855323,45.112829999999995,0.38092,2.061,1.073,RAVE-DR5.tsv
*/
