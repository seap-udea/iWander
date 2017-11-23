#include <iwander.cpp>
using namespace std;

#define VERBOSE 1

int main(int argc,char* argv[])
{
  /*
    Example: ./probability.exe

    Function: calculate the IOP for a list of candidates.

    Input:
    * potential.csv

    Output: 

  */

  ////////////////////////////////////////////////////
  //CONFIGURATION
  ////////////////////////////////////////////////////
  #include <probability.conf>

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initWander();
  
  ////////////////////////////////////////////////////
  //GENERAL VARIABLES
  ////////////////////////////////////////////////////
  double tmp;
  char ctmp[100],line[MAXLINE],values[MAXLINE];
  int Ntimes,Ntimesp,Nobs,nsys,nsysp; 
  int ip,n;
  int nfields;
  double params[10],mparams[23];
  double hstep,duration;
  //MATRICES WITH INTEGRATIONS
  double *xnom0,*xnoms0,*xIntp0,*xIntc0,**xIntp,**xFullp,**xFullc,**xIntc;
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

  int RAf,RA_ERRORf,DECf,DEC_ERRORf;
  int PMRAf,PMRA_ERRORf;
    
  enum{
    RA,
    RA_ERROR,
    RA_DEC_CORR,
    RA_PARALLAX_CORR,
    RA_PMRA_CORR,
    RA_PMDEC_CORR,
    DEC,
    DEC_ERROR,
    DEC_PARALLAX_CORR,
    DEC_PMRA_CORR,
    DEC_PMDEC_CORR,
    PMRA,
    PMRA_ERROR,
    PMRA_PMDEC_CORR,
    PMDEC,
    PMDEC_ERROR,
    PARALLAX,
    PARALLAX_ERROR,
    PARALLAX_PMRA_CORR,
    PARALLAX_PMDEC_CORR
  };
  int FC[]={
    Candidates::RA,
    Candidates::RA_ERROR,
    Candidates::RA_DEC_CORR,
    Candidates::RA_PARALLAX_CORR,
    Candidates::RA_PMRA_CORR,
    Candidates::RA_PMDEC_CORR,
    Candidates::DEC,
    Candidates::DEC_ERROR,
    Candidates::DEC_PARALLAX_CORR,
    Candidates::DEC_PMRA_CORR,
    Candidates::DEC_PMDEC_CORR,
    Candidates::PMRA,
    Candidates::PMRA_ERROR,
    Candidates::PMRA_PMDEC_CORR,
    Candidates::PMDEC,
    Candidates::PMDEC_ERROR,
    Candidates::PARALLAX,
    Candidates::PARALLAX_ERROR,
    Candidates::PARALLAX_PMRA_CORR,
    Candidates::PARALLAX_PMDEC_CORR
  };

  ////////////////////////////////////////////////////
  //UNITS
  ////////////////////////////////////////////////////
  VPRINT(stdout,"Units:\n\tUM = %.5e kg=%.5e Msun\n\tUL = %.17e m = %.17e pc\n\tUT = %.5e s = %.5e yr\n\tUV = %.5e m/s = %.5e km/s\n\tG = %.5e\n",
	 UM,UM/MSUN,UL,UL/PARSEC,UT,UT/YEAR,UV,UV/1e3,GGLOBAL);
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

  xFullp=(double**)malloc(Ntimesp*sizeof(double*));
  for(int j=0;j<Ntimesp;j++) xFullp[j]=(double*)malloc(nsysp*sizeof(double));

  xFullc=(double**)malloc(Ntimesp*sizeof(double*));
  for(int j=0;j<Ntimesp;j++) xFullc[j]=(double*)malloc(nsysp*sizeof(double));

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
  params[ip++]=GGLOBAL*MDISK*MSUN/UM;
  params[ip++]=ADISK*PARSEC/UL;
  params[ip++]=BDISK*PARSEC/UL;
  params[ip++]=GGLOBAL*MBULGE*MSUN/UM;
  params[ip++]=0.0;
  params[ip++]=BBULGE*PARSEC/UL;
  params[ip++]=GGLOBAL*MHALO*MSUN/UM;
  params[ip++]=0.0;
  params[ip++]=BHALO*PARSEC/UL;

  ////////////////////////////////////////////////////
  //READ PARTICLES POSITION
  ////////////////////////////////////////////////////
  FILE *fc;
  fc=fopen("wanderer.csv","r");
  fgets(line,MAXLINE,fc);//HEADER

  int i=0;
  fprintf(stdout,"Reading initial conditions of test particles:\n");
  while(fgets(line,MAXLINE,fc)!=NULL){
    parseLine(line,fields,&nfields);
    tsp[i]=0.0;
    ip=6*i;

    n=Wanderer::XGAL;
    for(int k=0;k<6;k++) x[k]=atof(fields[n++]);

    LSR2GC(x,xg);
    vscl_c(1e3/UL,xg,xg);//SET UNITS
    vscl_c(1e3/UV,xg+3,xg+3);
    //CONVERT TO CYLINDRICAL GALACTIC COORDINATES
    cart2polar(xg,xIntp0+ip,1.0);
    copyVec(xIntc0+ip,xIntp0+ip,6);

    fprintf(stdout,"\tInitial condition for particle %d: %s\n",i,vec2strn(xIntp0+ip,6,"%e "));
    i++;
    if(i==Npart) break;
  }
  fclose(fc);
  copyVec(xnoms0,xIntp0,6);
  fprintf(stdout,"Initial condition nominal test particle: %s\n",vec2strn(xIntp0,6,"%e "));
  
  ////////////////////////////////////////////////////
  //READ PREINTEGRATED TEST PARTICLES
  ////////////////////////////////////////////////////
  fprintf(stdout,"Reading preintegrated %d test particles:\n",Npart);
  fc=fopen("cloud.csv","r");
  fgets(line,MAXLINE,fc);//HEADER
  i=0;
  while(fgets(line,MAXLINE,fc)!=NULL){
    parseLine(line,fields,&nfields);
    tsp[i]=atof(fields[0]);
    n=1;
    for(int j=0;j<Npart;j++){
      ip=6*j;
      x=xFullp[i]+ip;
      for(int k=0;k<6;k++) x[k]=atof(fields[n++]);
      x=xFullc[i]+ip;
      for(int k=0;k<6;k++) x[k]=atof(fields[n++]);
    }
    i++;
  }
  fclose(fc);
  fprintf(stdout,"Initial condition nominal test particle (preintegrated): %s\n",vec2strn(xFullp[0],6,"%e "));
  
  ////////////////////////////////////////////////////
  //READING POTENTIAL OBJECTS
  ////////////////////////////////////////////////////
  params[0]=nsys;

  //CHOOSE H FOR PROBABILITY CALCULATIONS
  double hprob=1.0*PARSEC/UL;//PC
  double sigma=wNormalization(hprob);

  n=0;
  //fc=fopen("potential.csv","r");
  fc=fopen("candidates.csv","r");
  fgets(line,MAXLINE,fc);//HEADER

  FILE *fp=fopen("progenitors.csv","w");
  fprintf(fp,"Pprob,nomdmin,nomtmin,mindmin,maxdmin,velrelmin,velrelmax,%s",line);

  int qinterrupt=0;
  while(fgets(line,MAXLINE,fc)!=NULL){
    
    //LOCAL MINIMA
    double mindmin=1e100;
    double maxdmin=0.0;
    double velrelmin,velrelmax;

    //PARSE FIELDS
    strcpy(values,line);
    parseLine(line,fields,&nfields);
    n++;

    //if(strcmp(fields[TYCHO2_ID],"7774-308-1")!=0) continue;

    /*
    if(strcmp(fields[Candidates::HIP],"103749")!=0) continue;
    else qinterrupt=1;
    //*/

    fprintf(stdout,"Star %d,%s,%s:\n",n,fields[Potential::HIP],fields[Potential::TYCHO2_ID]);

    //ESTIMATED TIME OF ENCOUNTER
    /*
    tmin=atof(fields[Potential::DYNTMIN]);
    dmin=atof(fields[Potential::DYNDMIN]);
    */
    tmin=atof(fields[Candidates::TMIN]);
    dmin=atof(fields[Candidates::DMIN]);

    double tmin0=tmin;
    double tmins=tmin0;
    double dmin0=dmin;
    double nomtmin,nomdmin;

    mint=1.2*tmin;
    maxt=0.8*tmin;

    fprintf(stdout,"\tEncounter estimated time, tmin = %.6e\n",tmin);
    fprintf(stdout,"\t\tTesting range = [%.6e,%.6e]\n",mint,maxt);
    fprintf(stdout,"\tEncounter estimated distance, dmin = %.6e\n",dmin);
    
    if(!(fabs(tmin0)>0)){
      fprintf(stdout,"This star is too close\n");
      //getchar();
      continue;
    }
    
    if(fields[Candidates::RA][0]=='N'){
      fprintf(stdout,"*************** USING DATA FROM HIPPARCOS ****************\n");
      int FChip[]={
	Candidates::RA_HIP,
	Candidates::RA_ERROR_HIP,
	Candidates::RA_DEC_CORR_HIP,
	Candidates::RA_PARALLAX_CORR_HIP,
	Candidates::RA_PMRA_CORR_HIP,
	Candidates::RA_PMDEC_CORR_HIP,
	Candidates::DEC_HIP,
	Candidates::DEC_ERROR_HIP,
	Candidates::DEC_PARALLAX_CORR_HIP,
	Candidates::DEC_PMRA_CORR_HIP,
	Candidates::DEC_PMDEC_CORR_HIP,
	Candidates::PMRA_HIP,
	Candidates::PMRA_ERROR_HIP,
	Candidates::PMRA_PMDEC_CORR_HIP,
	Candidates::PMDEC_HIP,
	Candidates::PMDEC_ERROR_HIP,
	Candidates::PARALLAX_HIP,
	Candidates::PARALLAX_ERROR_HIP,
	Candidates::PARALLAX_PMRA_CORR_HIP,
	Candidates::PARALLAX_PMDEC_CORR
      };
      copyVecInt(FC,FChip,20);
    }
    /*
    else{
      continue;
    }
    */

    //INFORMATION REQUIRED
    mobs[0]=ra=atof(fields[FC[RA]]);
    dra=atof(fields[FC[RA_ERROR]])*MAS;

    mobs[1]=dec=atof(fields[FC[DEC]]);
    ddec=atof(fields[FC[DEC_ERROR]])*MAS;

    mobs[2]=par=atof(fields[FC[PARALLAX]]);
    dpar=atof(fields[FC[PARALLAX_ERROR]]);

    mobs[3]=mura=atof(fields[FC[PMRA]]);
    dmura=atof(fields[FC[PMRA_ERROR]]);

    mobs[4]=mudec=atof(fields[FC[PMDEC]]);
    dmudec=atof(fields[FC[PMDEC_ERROR]]);

    mobs[5]=vr=atof(fields[Candidates::RV]);
    dvr=atof(fields[Candidates::ERV]);

    //COVARIANCE MATRIX
    /*RA*/cov[0][0]=dra*dra;
    cov[0][1]=atof(fields[FC[RA_DEC_CORR]])*dra*ddec;
    cov[0][2]=atof(fields[FC[RA_PARALLAX_CORR]])*dra*dpar;
    cov[0][3]=atof(fields[FC[RA_PMRA_CORR]])*dra*dmura;
    cov[0][4]=atof(fields[FC[RA_PMDEC_CORR]])*dra*dmudec;
    cov[0][5]=0.0;
    /*DEC*/cov[1][1]=ddec*ddec;
    cov[1][0]=cov[0][1];
    cov[1][2]=atof(fields[FC[DEC_PARALLAX_CORR]])*ddec*dpar;
    cov[1][3]=atof(fields[FC[DEC_PMRA_CORR]])*ddec*dmura;
    cov[1][4]=atof(fields[FC[DEC_PMDEC_CORR]])*ddec*dmudec;
    cov[1][5]=0.0;
    /*PAR*/cov[2][2]=dpar*dpar;
    cov[2][0]=cov[0][2];
    cov[2][1]=cov[1][2];
    cov[2][3]=atof(fields[FC[PARALLAX_PMRA_CORR]])*dpar*dmura;
    cov[2][4]=atof(fields[FC[PARALLAX_PMDEC_CORR]])*dpar*dmudec;
    cov[2][5]=0.0;
    /*MURA*/cov[3][3]=dmura*dmura;
    cov[3][0]=cov[0][3];
    cov[3][1]=cov[1][3];
    cov[3][2]=cov[2][3];
    cov[3][4]=atof(fields[FC[PMRA_PMDEC_CORR]])*dmura*dmudec;
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
	   atof(fields[Candidates::L]),atof(fields[Candidates::B]));

    VPRINT(stdout,"\tStar Covariance Matrix:\n");
    for(int i=0;i<6;i++)
      fprintf(stdout,"\t\t|%s|\n",vec2strn(cov[i],6,"%-+15.3e"));

    generateMultivariate(cov,mobs,obs,6,Nobs);

    VPRINT(stdout,"\tSurrogate random properties:\n");
    for(int i=Nobs;i-->0;){
      VPRINT(stdout,"\t\tObservation %d: %s\n",i,vec2strn(obs[i],6,"%.10e "));
    }

    //COMPUTE THE NOMINAL CONDITIONS
    d=AU/tan(par/(60*60*1000.0)*DEG)/PARSEC;
    radrec_c(d,ra*DEG,dec*DEG,xg);
    mxv_c(M_J2000_Galactic,xg,x);
    recrad_c(x,&tmp,&l,&b);
    calcUVW(ra,dec,par,dpar,mura,dmura,mudec,dmudec,vr,dvr,x+3,dx+3);
    vscl_c(PARSEC/1e3,x,x);
    LSR2GC(x,xg);
    vscl_c(1e3/UL,xg,xg);//SET UNITS
    vscl_c(1e3/UV,xg+3,xg+3);
    cart2polar(xg,xInt0,1.0);
    params[0]=6;
    copyVec(xnom0,xnoms0,6);

    //DISCRETE METHOD
    //minDistanceDiscrete(xInt0,xFullc,tsp,Ntimesp,1.1*tmin0,params,&dmin,&tmin);

    try{
      minDistance2(xInt0,xnom0,tmin0,&dmin,&tmin,params);
    }catch(int e){
      fprintf(stdout,"¡No minimum!\n");
      //continue;
    }
    fprintf(stdout,"\t\tMinimum distance from nominal to nominal (nom. t=%.6e, d=%.6e): t = %.6e, d = %.6e\n",tmin0,dmin0,tmin,dmin);
    nomdmin=dmin;
    nomtmin=tmin;

    /*
    int it=findTime(tmin,tsp,Ntimesp);
    fprintf(stdout,"\t\tTime %e correspond to interval %d: [%e,%e]\n",tmin,it,tsp[it],tsp[it+1]);
    double *xp1=xFullc[it];
    double *xp2=xFullc[it+1];
    fprintf(stdout,"\t\tNominal particle position 1 (%e):%s\n",tsp[it],vec2strn(xp1,6,"%.5e "));
    fprintf(stdout,"\t\tNominal particle position 2 (%e):%s\n",tsp[it+1],vec2strn(xp2,6,"%.5e "));
    */

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
	//getchar();
	continue;
      }
      fprintf(stdout,"\t\tMinimum distance at (nom. d=%.6e, t=%.6e): t = %.6e, d = %.6e\n",tmin0,dmin0,tmin,dmin);
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

	//FIND MINIMUM APPROXIMATION
	if(D<mindmin){
	  mindmin=D;
	  velrelmin=vrel;
	}
	if(D>maxdmin){
	  maxdmin=D;
	  velrelmax=vrel;
	}

	fprintf(stdout,"\t\t\tDistance to test particle %d: d=%.6e,vrel=%.6e\n",j,D,vrel*UV/1e3);
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
    fprintf(stdout,"Probability for star %d: %.6e\n",n,Pprob);
    fprintf(stdout,"Minimum distance (lma.dmin=%e,nom.dmin=%e): min.dmin = %e, max.dmin = %e\n",dmin0,nomdmin,mindmin,maxdmin);

    fprintf(fp,"%e,%e,%e,%e,%e,%e,%e,",Pprob,nomdmin,nomtmin,mindmin,maxdmin,velrelmin*UV/1e3,velrelmax*UV/1e3);
    fprintf(fp,"%s",values);
    fflush(fp);

    //if(VERBOSE) getchar();
    if(qinterrupt) break;
  }
  fclose(fc);
  fclose(fp);

  return 0;
}
