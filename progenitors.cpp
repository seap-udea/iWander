#include <gravray.cpp>
using namespace std;

#define VERBOSE 1

//FIELDS OF THE POTENTIAL.CSV FILE
enum {DYNTMIN,DYNDMIN,DYNVREL,N,POSTARX,POSTARY,POSTARZ,VELSTARX,VELSTARY,VELSTARZ,POSBODYPERIX,POSBODYPERIY,POSBODYPERIZ,POSTARPERIX,POSTARPERIY,POSTARPERIZ,DMIN,TMIN,VRELX,VRELY,VRELZ,VREL,HIP,TYCHO2_ID,SOLUTION_ID,SOURCE_ID,RANDOM_INDEX,REF_EPOCH,RA,RA_ERROR,DEC,DEC_ERROR,PARALLAX,PARALLAX_ERROR,PMRA,PMRA_ERROR,PMDEC,PMDEC_ERROR,RA_DEC_CORR,RA_PARALLAX_CORR,RA_PMRA_CORR,RA_PMDEC_CORR,DEC_PARALLAX_CORR,DEC_PMRA_CORR,DEC_PMDEC_CORR,PARALLAX_PMRA_CORR,PARALLAX_PMDEC_CORR,PMRA_PMDEC_CORR,ASTROMETRIC_N_OBS_AL,ASTROMETRIC_N_OBS_AC,ASTROMETRIC_N_GOOD_OBS_AL,ASTROMETRIC_N_GOOD_OBS_AC,ASTROMETRIC_N_BAD_OBS_AL,ASTROMETRIC_N_BAD_OBS_AC,ASTROMETRIC_DELTA_Q,ASTROMETRIC_EXCESS_NOISE,ASTROMETRIC_EXCESS_NOISE_SIG,ASTROMETRIC_PRIMARY_FLAG,ASTROMETRIC_RELEGATION_FACTOR,ASTROMETRIC_WEIGHT_AL,ASTROMETRIC_WEIGHT_AC,ASTROMETRIC_PRIORS_USED,MATCHED_OBSERVATIONS,DUPLICATED_SOURCE,SCAN_DIRECTION_STRENGTH_K1,SCAN_DIRECTION_STRENGTH_K2,SCAN_DIRECTION_STRENGTH_K3,SCAN_DIRECTION_STRENGTH_K4,SCAN_DIRECTION_MEAN_K1,SCAN_DIRECTION_MEAN_K2,SCAN_DIRECTION_MEAN_K3,SCAN_DIRECTION_MEAN_K4,PHOT_G_N_OBS,PHOT_G_MEAN_FLUX,PHOT_G_MEAN_FLUX_ERROR,PHOT_G_MEAN_MAG,PHOT_VARIABLE_FLAG,L,B,ECL_LON,ECL_LAT,RAJ2000,DEJ2000,RV,ERV,CAT};

int main(int argc,char* argv[])
{
  /*
    Example: ./probability.exe 10

    Where: 
    10: Number of particles

    Input: 
    * cloud.data: test particles compatible with A/2017U1 orbit.
    * potential.csv: candidates to close encounters.

    Output:
    *
  */

  ////////////////////////////////////////////////////
  //INPUT PARAMETERS
  ////////////////////////////////////////////////////
  int Npart=1;
  if(argc>1){
    Npart=atoi(argv[1]);
  }
  VPRINT(stdout,"Analysing %d test particles against candidate stars\n",Npart);
  VPRINT(stdout,"\n");

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initSpice();
  
  ////////////////////////////////////////////////////
  //GENERAL VARIABLES
  ////////////////////////////////////////////////////
  double tmp;
  char ctmp[100],line[MAXLINE];
  int Ntimes,Ntimesp,Nobs,nsys,nsysp; 
  int ip,n;
  int nfields;
  double params[10],mparams[23];
  double hstep,duration;
  //MATRICES WITH INTEGRATIONS
  double *xnom0,*xnoms0,*xIntp0,*xIntc0,**xIntp,**xIntc;
  double *xInt0,**xInt;
  double *tsp,*ts;
  //INITIAL CONDITIONS
  double *dxIntdt,*x,*xg,*xp1,*xp2,*xpmin,*dx,*x0;
  double dmin,tmin,ftmin,dyn_tmin,dyn_dmin,dyn_vrel;
  double G;
  double t;
  int it;
  double Pprob,Psur,fvel,fdist;

  //     0  1   2   3    4     5 
  double ra,dec,par,mura,mudec,vr;
  double dra,ddec,dpar,dmura,dmudec,dvr;
  double **cov,**obs,*mobs;

  double UVW[3];
  SpiceDouble M_J2000_Galactic[3][3];
  double d,l,b,h;
  double mint,maxt;
  double ps[23];

  ////////////////////////////////////////////////////
  //UNITS
  ////////////////////////////////////////////////////
  UM=MSUN;
  UL=PARSEC;
  UT=YEAR;
  UV=UL/UT;
  G=GCONST/(UL*UL*UL/(UM*UT*UT));
  
  VPRINT(stdout,"Units:\n\tUM = %.5e kg=%.5e Msun\n\tUL = %.17e m = %.17e pc\n\tUT = %.5e s = %.5e yr\n\tUV = %.5e m/s = %.5e km/s\n\tG = %.5e\n",
	 UM,UM/MSUN,UL,UL/PARSEC,UT,UT/YEAR,UV,UV/1e3,G);
  VPRINT(stdout,"\n");

  ////////////////////////////////////////////////////
  //GLOBAL ALLOCATION
  ////////////////////////////////////////////////////
  nsysp=6*Npart;
  nsys=6;

  Ntimesp=10000;
  Ntimes=100;
  Nobs=10;

  x=(double*)malloc(6*sizeof(double));//LSR STATE VECTOR
  xg=(double*)malloc(6*sizeof(double));//GC STATE VECTOR
  dx=(double*)malloc(6*sizeof(double));//LSR STATE VECTOR
  xpmin=(double*)malloc(6*sizeof(double));//GC STATE VECTOR

  xIntp0=(double*)malloc(nsysp*sizeof(double));
  xIntc0=(double*)malloc(nsysp*sizeof(double));
  xInt0=(double*)malloc(nsys*sizeof(double));
  xnom0=(double*)malloc(nsysp*sizeof(double));
  xnoms0=(double*)malloc(nsysp*sizeof(double));
  
  xInt=(double**)malloc(Ntimes*sizeof(double*));
  for(int j=0;j<Ntimes;j++) xInt[j]=(double*)malloc(nsys*sizeof(double));

  xIntp=(double**)malloc(Ntimesp*sizeof(double*));
  for(int j=0;j<Ntimesp;j++) xIntp[j]=(double*)malloc(nsysp*sizeof(double));

  xIntc=(double**)malloc(Ntimesp*sizeof(double*));
  for(int j=0;j<Ntimesp;j++) xIntc[j]=(double*)malloc(nsysp*sizeof(double));
  
  char **fields=(char**)malloc(MAXCOLS*sizeof(char*));
  for(int i=0;i<MAXCOLS;i++) fields[i]=(char*)malloc(MAXTEXT*sizeof(char));

  tsp=(double*)malloc(Ntimesp*sizeof(double));
  ts=(double*)malloc(Ntimesp*sizeof(double));

  mobs=(double*)malloc(6*sizeof(double));

  obs=(double**)malloc(Nobs*sizeof(double*));
  for(int i=0;i<Nobs;i++) obs[i]=(double*)malloc(6*sizeof(double));

  cov=(double**)malloc(6*sizeof(double*));
  for(int i=0;i<6;i++) cov[i]=(double*)malloc(6*sizeof(double));

  pxform_c("J2000","GALACTIC",0,M_J2000_Galactic);

  ////////////////////////////////////////////////////
  //GLOBAL PROPERTIES FOR INTEGRATION
  ////////////////////////////////////////////////////
  ip=1;
  params[ip++]=G*MDISK*MSUN/UM;
  params[ip++]=ADISK*PARSEC/UL;
  params[ip++]=BDISK*PARSEC/UL;
  params[ip++]=G*MBULGE*MSUN/UM;
  params[ip++]=0.0;
  params[ip++]=BBULGE*PARSEC/UL;
  params[ip++]=G*MHALO*MSUN/UM;
  params[ip++]=0.0;
  params[ip++]=BHALO*PARSEC/UL;

  ////////////////////////////////////////////////////
  //READ PARTICLES POSITION
  ////////////////////////////////////////////////////
  FILE *fc;
  if((fc=fopen("cloud-int.csv","r"))==NULL){
    fprintf(stdout,"Houston we've got a problem\n");
    exit(0);
  }
  fgets(line,MAXLINE,fc);//HEADER
  int i=0;
  while(fgets(line,MAXLINE,fc)!=NULL){
    parseLine(line,fields,&nfields);
    tsp[i]=atof(fields[0]);
    n=1;
    for(int j=0;j<Npart;j++){
      ip=6*j;
      x=xIntp0+ip;
      for(int k=0;k<6;k++) x[k]=atof(fields[n++]);
      x=xIntc0+ip;
      for(int k=0;k<6;k++) x[k]=atof(fields[n++]);
    }
    //READ ONLY INITIAL CONDITIONS
    break;
  }
  fclose(fc);
  copyVec(xnoms0,xIntp0,6);

  ////////////////////////////////////////////////////
  //READING POTENTIAL OBJECTS
  ////////////////////////////////////////////////////
  params[0]=nsys;
  int Nsur=10;

  //CHOOSE H FOR PROBABILITY CALCULATIONS
  double hprob=1.0*PARSEC/UL;//PC
  double sigma=wNormalization(hprob);

  n=0;
  fc=fopen("potential.csv","r");
  fgets(line,MAXLINE,fc);//HEADER

  while(fgets(line,MAXLINE,fc)!=NULL){
    
    //PARSE FIELDS
    parseLine(line,fields,&nfields);
    n++;
    //if(strcmp(fields[TYCHO2_ID],"7774-308-1")!=0) continue;

    fprintf(stdout,"Star %d,%s,%s:\n",n,fields[HIP],fields[TYCHO2_ID]);

    //ESTIMATED TIME OF ENCOUNTER
    tmin=atof(fields[DYNTMIN]);
    dmin=atof(fields[DYNDMIN]);
    double tmin0=tmin;
    double tmins=tmin0;

    mint=1.2*tmin;
    maxt=0.8*tmin;

    fprintf(stdout,"\tEncounter estimated time, tmin = %.6e\n",tmin);
    fprintf(stdout,"\t\tTesting range = [%.6e,%.6e]\n",mint,maxt);
    fprintf(stdout,"\tEncounter estimated distance, dmin = %.6e\n",dmin);

    if(!(fabs(tmin0)>0)){
      fprintf(stdout,"This star is too close\n");
      getchar();
      continue;
    }
    
    //INFORMATION REQUIRED
    mobs[0]=ra=atof(fields[RA]);
    dra=atof(fields[RA_ERROR])*MAS;

    mobs[1]=dec=atof(fields[DEC]);
    ddec=atof(fields[DEC_ERROR])*MAS;

    mobs[2]=par=atof(fields[PARALLAX]);
    dpar=atof(fields[PARALLAX_ERROR]);

    mobs[3]=mura=atof(fields[PMRA]);
    dmura=atof(fields[PMRA_ERROR]);

    mobs[4]=mudec=atof(fields[PMDEC]);
    dmudec=atof(fields[PMDEC_ERROR]);

    mobs[5]=vr=atof(fields[RV]);
    dvr=atof(fields[ERV]);

    //COVARIANCE MATRIX
    /*RA*/cov[0][0]=dra*dra;
    cov[0][1]=atof(fields[RA_DEC_CORR])*dra*ddec;
    cov[0][2]=atof(fields[RA_PARALLAX_CORR])*dra*dpar;
    cov[0][3]=atof(fields[RA_PMRA_CORR])*dra*dmura;
    cov[0][4]=atof(fields[RA_PMDEC_CORR])*dra*dmudec;
    cov[0][5]=0.0;
    /*DEC*/cov[1][1]=ddec*ddec;
    cov[1][0]=cov[0][1];
    cov[1][2]=atof(fields[DEC_PARALLAX_CORR])*ddec*dpar;
    cov[1][3]=atof(fields[DEC_PMRA_CORR])*ddec*dmura;
    cov[1][4]=atof(fields[DEC_PMDEC_CORR])*ddec*dmudec;
    cov[1][5]=0.0;
    /*PAR*/cov[2][2]=dpar*dpar;
    cov[2][0]=cov[0][2];
    cov[2][1]=cov[1][2];
    cov[2][3]=atof(fields[PARALLAX_PMRA_CORR])*dpar*dmura;
    cov[2][4]=atof(fields[PARALLAX_PMDEC_CORR])*dpar*dmudec;
    cov[2][5]=0.0;
    /*MURA*/cov[3][3]=dmura*dmura;
    cov[3][0]=cov[0][3];
    cov[3][1]=cov[1][3];
    cov[3][2]=cov[2][3];
    cov[3][4]=atof(fields[PMRA_PMDEC_CORR])*dmura*dmudec;
    cov[3][5]=0.0;
    /*MUDEC*/cov[4][4]=dmudec*dmudec;
    cov[4][0]=cov[0][4];
    cov[4][1]=cov[1][4];
    cov[4][2]=cov[2][4];
    cov[4][3]=cov[3][4];
    cov[4][5]=0.0;
    /*RV*/cov[5][5]=dvr*dvr;
    cov[5][0]=cov[0][5];
    cov[5][1]=cov[1][5];
    cov[5][2]=cov[2][5];
    cov[5][3]=cov[3][5];
    cov[5][4]=cov[4][5];
    
    VPRINT(stdout,"\tStellar properties: %s\n",vec2strn(mobs,6,"%.5e "));
    VPRINT(stdout,"\t\tErrors:");
    for(int i=0;i<6;i++) VPRINT(stdout,"%.6e ",sqrt(cov[i][i]));
    VPRINT(stdout,"\n");
    VPRINT(stdout,"\tGalactic coordinates: l = %lf, b = %lf\n",
	   atof(fields[L]),atof(fields[B]));

    VPRINT(stdout,"\tStar Covariance Matrix:\n");
    for(int i=0;i<6;i++)
      fprintf(stdout,"\t\t|%s|\n",vec2strn(cov[i],6,"%-+15.3e"));

    generateMultivariate(cov,mobs,obs,6,Nobs);

    VPRINT(stdout,"\tSurrogate random properties:\n");
    for(int i=Nobs;i-->0;){
      VPRINT(stdout,"\t\tObservation %d: %s\n",i,vec2strn(obs[i],6,"%.10e "));
    }

    //CALCULATE PROBABILITIES
    Pprob=0;
    for(int i=0;i<Nsur;i++){

      VPRINT(stdout,"\tSurrogate %d:\n",i);
      copyVec(xnom0,xnoms0,6);

      //GET KEY PROPERTIES OF SURROGATE
      ra=obs[i][0];
      dec=obs[i][1];
      par=obs[i][2];
      mura=obs[i][3];
      mudec=obs[i][4];
      vr=obs[i][5];

      //INITIAL POSITION RELATIVE TO SUN
      d=AU/tan(par/(60*60*1000.0)*DEG)/PARSEC;
      radrec_c(d,ra*DEG,dec*DEG,xg);
      mxv_c(M_J2000_Galactic,xg,x);
      recrad_c(x,&tmp,&l,&b);
      calcUVW(ra,dec,par,dpar,mura,dmura,mudec,dmudec,vr,dvr,x+3,dx+3);
      //INITIAL POSITION RELATIVE TO GALACTIC CENTER
      vscl_c(PARSEC/1e3,x,x);
      LSR2GC(x,xg);
      vscl_c(1e3/UL,xg,xg);//SET UNITS
      vscl_c(1e3/UV,xg+3,xg+3);
      //INITIAL POLAR COORDINATES OF SURROGATE
      cart2polar(xg,xInt0,1.0);

      VPRINT(stdout,"\t\ttmin0: %.6e\n",tmin0);
      VPRINT(stdout,"\t\tObservations: %s\n",vec2strn(obs[i],6,"%.5e "));
      VPRINT(stdout,"\t\tDistance: %e pc\n",d);
      VPRINT(stdout,"\t\tGalactic coordinates: l = %lf, b = %lf\n",l*RAD,b*RAD);
      VPRINT(stdout,"\t\tInitial position cartesian: %s\n",vec2strn(x,6,"%.5e "));
      VPRINT(stdout,"\t\tInitial position cylindrical: %s\n",vec2strn(xnom0,6,"%.5e "));

      //CALCULATE MINIMUM DISTANCE AND TIME OF NOMINAL SOLUTION TO SURROGATE
      params[0]=6;

      try{
	minDistance2(xInt0,xnom0,tmin0,&dmin,&tmin,params);
      }catch(int e){
	fprintf(stdout,"¡No minimum!\n");
	getchar();
	continue;
      }
      fprintf(stdout,"\t\tMinimum distance at: t = %.6e, d = %.6e\n",tmin,dmin);
      VPRINT(stdout,"\t\tFinal position star: %s\n",vec2strn(xInt0,6,"%.5e "));
      VPRINT(stdout,"\t\tFinal position nominal: %s\n",vec2strn(xnom0,6,"%.5e "));

      //PROPAGATE ALL TEST PARTICLE
      params[0]=6*Npart;
      h=fabs(tmin)/100;
      VPRINT(stdout,"\t\tInitial conditions for all particles: %s\n",vec2strn(xIntp0,6*Npart,"%.5e "));
      integrateEoM(0,xIntp0,h,2,tmin,6*Npart,EoMGalactic,params,ts,xIntp);
      VPRINT(stdout,"\t\tIntegration result for all particles: %s\n",vec2strn(xIntp[1],6*Npart,"%.5e "));

      double D,Dmax=0,*xt1,*xt2,vrel;
      
      //COMPUTE SPH-LIKE PROBABILITY
      double pd,pv;
      Psur=0.0;
      fvel=0.0;
      fprintf(stdout,"\t\tComparing test particle position with star position\n");
      for(int j=0;j<Npart;j++){
	vsubg_c(xInt0,xIntp[1]+6*j,6,dx);
	D=vnorm_c(dx);
	vrel=vnorm_c(dx+3);
	fprintf(stdout,"\t\t\tDistance to test particle %d: v=%.6e,vrel=%.6e\n",j,D,vrel*UV/1e3);
	//CONTRIBUTION TO P FROM DISTANCE
	pd=sigma*wFunction(D,&hprob);

	//CONTRIBUTION TO P FROM VELOCITY
	pv=1.0;
	
	fprintf(stdout,"\t\t\t\tDistance probability: %.6e\n",pd);
	fprintf(stdout,"\t\t\t\tVelocity probability: %.6e\n",pv);

	Psur+=pd*pv;
	//COMPUTE CORRECTION FOR RELATIVE STELLAR VELOCITY
      }

      //COMPUTE CORRECTION FOR STELLAR DISTANCE
      fdist=1/(d*d);
      Psur*=fdist;
      
      //SURROGATE PROBABILITY
      fprintf(stdout,"\t\tSurrogate probability: %.6e\n",Psur);

      //ACCUMULATE
      Pprob+=Psur;
      //getchar();
    }
    fprintf(stdout,"Probability for star: %.6e\n",Pprob);
    getchar();
    //break;
  }
  fclose(fc);
    
  return 0;
}
