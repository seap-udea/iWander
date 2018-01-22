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
# Compute encounter probabilities
#########################################################
*/
#include <iwander.cpp>
using namespace std;

#define VERBOSE 1

struct vinfpar {
  double xmin,xmax;
  gsl_interp_accel *a;
  gsl_spline *s;
};
double vinfPosterior(double x,void *params)
{
  struct vinfpar *pars=(struct vinfpar*)params;
  double p;
  if(x<pars->xmin) x=1.01*pars->xmin;
  if(x>pars->xmax){
    p=12*pow10(-3*x);
    /*
    VPRINT(stdout,"Extrapolation for %e: %e\n",x,p);
    p=gsl_spline_eval(pars->s,pars->xmax,pars->a);
    VPRINT(stdout,"Compare with maximum: %e\n",p);
    */
  }else{
    p=gsl_spline_eval(pars->s,x,pars->a);
  }
  return p;
}
double vinfProbability(double v1,double v2,void *params)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  gsl_function F={.function=vinfPosterior,.params=params};
  double vint,vint_error;
  gsl_integration_qags(&F,v1,v2,0,1e-5,1000,w,&vint,&vint_error);
  return vint;
}

int main(int argc,char* argv[])
{
  /*
    Example: ./probability.exe

    Function: calculate the IOP for a list of candidates.

    Input:
    * wanderer.csv
    * candidates.csv

    Output: 
    * progenitors.csv

      Cols:
          0: IOP for this candidate, Pprob
	  1: Average position probability, Psmed
	  2: Average velocity probability, Pvmed
	  3: Probability distance factor, fdist
	  4-6: Nominal minimum time, minimum distance, relative velocity
	  7,8: Minimum and maximum tmin
	  9,10: Minimum and maximum dmin
	  11,12: Minimum and maximum vrel
	  13...: Same as candidates.csv
  */
  ////////////////////////////////////////////////////
  //CONFIGURATION
  ////////////////////////////////////////////////////
  #include <iwander.conf>
  #include <probability.conf>

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initWander();

  ////////////////////////////////////////////////////
  //LOAD COMMAND LINE OPTIONS
  ////////////////////////////////////////////////////
  int ipart=-1;
  char basename[100];
  sprintf(basename,"scratch/candidates-%s.csv",WANDERER);
  if(argc>1){
    strcpy(basename,argv[1]);
    if(argc>2) ipart=atoi(argv[2]);
    else ipart=1;
    fprintf(stdout,"Analysing '%s' part %d\n",basename,ipart);
  }else{
    fprintf(stdout,"Analysing whole candidate set with basename '%s'\n",basename);
  }

  //REPORT THAT THE COMPUTATION HAS STARTED
  if(ipart<0){
    sprintf(Filename,"touch log/start-%s",WANDERER);
    system(Filename);
  }else{
    sprintf(Filename,"touch log/start-%s.%05d",WANDERER,ipart);
    system(Filename);
  }
  
  ////////////////////////////////////////////////////
  //GENERAL VARIABLES
  ////////////////////////////////////////////////////
  double nomtmin,nomdmin,nomvrel,telaps;
  double D,Dmax=0,*xt1,*xt2,vrel,Pvmed,Psmed;
  double tmp,ting;
  char ctmp[100],line[MAXLINE],values[MAXLINE];
  int Ntimes,Ntimesp,Nobs,nsys,nsysp,Ntelaps; 
  int ip,n;
  int nfields;
  double params[10],mparams[23];
  double hstep;
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

  int Nsur_acc;
    
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
  nsysp=6*Ntest;
  nsys=6;

  Ntimesp=10000;
  Ntimes=100;
  Nobs=Nsur;

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
  fprintf(stdout,"Params = %s\n",vec2strn(params,10,"%e "));

  ////////////////////////////////////////////////////
  //READ PARTICLES POSITION
  ////////////////////////////////////////////////////
  printHeader(stdout,"READING SURROGATE OBJECTS");
  FILE *fc;

  sprintf(Filename,"scratch/wanderer-%s.csv",WANDERER);
  fc=fopen(Filename,"r");

  fgets(line,MAXLINE,fc);//HEADER

  int i=0;
  VPRINT(stdout,"Reading initial conditions of test particles:\n");
  while(fgets(line,MAXLINE,fc)!=NULL){
    parseLine(line,fields,&nfields);
    tsp[i]=0.0;
    ip=6*i;

    ting=atof(fields[Wanderer::TING])/UT;

    n=Wanderer::XGAL;
    for(int k=0;k<6;k++) x[k]=atof(fields[n++]);

    VPRINT(stdout,"\t\tHeliocentric (km,km/s): %s\n",vec2strn(x,6,"%e "));
    
    LSR2GC(x,xg);

    VPRINT(stdout,"\t\tGalactocentric (km,km/s): %s\n",vec2strn(xg,6,"%e "));

    vscl_c(1e3/UL,xg,xg);//SET UNITS
    vscl_c(1e3/UV,xg+3,xg+3);


    VPRINT(stdout,"\t\tGalactocentric (UL,UV): %s\n",vec2strn(xg,6,"%e "));

    //CONVERT TO CYLINDRICAL GALACTIC COORDINATES
    cart2polar(xg,xIntp0+ip,1.0);
    copyVec(xIntc0+ip,xIntp0+ip,6);

    VPRINT(stdout,"\tInitial condition for particle %d: %s\n",i,vec2strn(xIntp0+ip,6,"%e "));
    i++;
    if(i==Ntest) break;
  }
  Ntest=i;
  fprintf(stdout,"Ntest = %d\n",Ntest);
  fclose(fc);
  copyVec(xnoms0,xIntp0,6);
  VPRINT(stdout,"Initial condition nominal test particle: %s\n",vec2strn(xIntp0,6,"%e "));

  ////////////////////////////////////////////////////
  //READ POSTERIOR EJECTION VELOCITY DISTRIBUTION
  ////////////////////////////////////////////////////
  FILE *fv=fopen("db/ejection-posterior.data","r");
  double ul,um,ut,uv,vcan,vradiuscan;
  double pvs[MAXCOLS],vinfs[MAXCOLS];
  int nv=0;
  while(fgets(line,MAXLINE,fv)!=NULL){
    parseLine(line,fields,&nfields," ");
    vinfs[nv]=atof(fields[0]);
    pvs[nv]=atof(fields[1]);
    nv++;
  }
  gsl_interp_accel *acc=gsl_interp_accel_alloc();
  gsl_spline *spline=gsl_spline_alloc(gsl_interp_cspline,nv);
  gsl_spline_init(spline,vinfs,pvs,nv);
  struct vinfpar vpar={.xmin=0,.xmax=vinfs[nv-1],.a=acc,.s=spline};
  VPRINT(stdout,"P(vinf=%e) = %e\n",0.2,vinfPosterior(0.2,&vpar));
  VPRINT(stdout,"Integral = %e +/- %e\n",vinfProbability(0,10,&vpar));

  ////////////////////////////////////////////////////
  //READING POTENTIAL OBJECTS
  ////////////////////////////////////////////////////
  params[0]=nsys;

  //CHOOSE H FOR PROBABILITY CALCULATIONS
  double hprob=1.0*PARSEC/UL;//PC
  double sigma=wNormalization(hprob);

  n=0;

  if(ipart<0){
    strcpy(Filename,basename);
  }else{
    sprintf(Filename,"%s.%05d",basename,ipart);
  }
  if((fc=fopen(Filename,"r"))==NULL){
    fprintf(stderr,"No candidates file '%s'\n",Filename);
  }
  fprintf(stdout,"Reading input filename %s...\n",Filename);
  fgets(line,MAXLINE,fc);//HEADER

  if(ipart<0){
    sprintf(Filename,"scratch/progenitors-%s.csv",WANDERER);
  }else{
    sprintf(Filename,"scratch/progenitors-%s.csv.%05d",WANDERER,ipart);
  }
  fprintf(stdout,"Generationg output filename %s...\n",Filename);
  FILE *fp=fopen(Filename,"w");

  fprintf(fp,"Pprob,Psurmed,Pvelmed,Pdist,nomtmin,nomdmin,nomvrel,tminl,tminmed,tminu,dminl,dminmed,dminu,vrell,vrelmed,vrelu,%s",line);

  int qinterrupt=0;

  printHeader(stdout,"CALCULATING PROBABILITIES");
  
  telaps=0.0;
  Ntelaps=0.0;
  elapsedTime();

  double dmins[Ntest*Nsur],tmis[Ntest*Nsur],vrels[Ntest*Nsur];
  while(fgets(line,MAXLINE,fc)!=NULL){
    
    //LOCAL MINIMA
    int nt;
    double mindmin=1e100,maxdmin=-1e100;
    double minvrel=1e100,maxvrel=-1e100;
    double mintmin=1e100,maxtmin=-1e100;

    //PARSE FIELDS
    strcpy(values,line);
    parseLine(line,fields,&nfields);
    n++;

    //if(strcmp(fields[TYCHO2_ID],"7774-308-1")!=0) continue;
    /*
    if(strcmp(fields[Candidates::HIP],"103749")!=0) continue;
    else qinterrupt=1;
    //*/
    /*
    if(strcmp(fields[Candidates::HIP],hip_single)!=0) continue;
    else qinterrupt=1;
    //*/
    //*
    if(qsingle){
      if(strlen(hip_single)==0)
	if(strcmp(fields[Candidates::TYCHO2_ID],tyc_single)!=0) continue;
	else qinterrupt=1;
      else
	if(strcmp(fields[Candidates::HIP],hip_single)!=0) continue;
	else qinterrupt=1;
    }
    //*/
    //if(n<239) continue;
    fprintf(stdout,"Computing probabilities for candidate star %d (%s,%s)...\n",
	    n,fields[Candidates::HIP],fields[Candidates::TYCHO2_ID]);

    //ESTIMATED TIME OF ENCOUNTER
    /*
    tmin=atof(fields[Candidates::DYNTMIN]);
    dmin=atof(fields[Candidates::DYNDMIN]);
    */
    tmin=atof(fields[Candidates::TMIN]);
    dmin=atof(fields[Candidates::DMIN]);

    fprintf(stdout,"\tEstimated nominal time and distance: dmin=%e, tmin=%e\n",dmin,tmin);

    double tmin0=tmin;
    double tmins=tmin0;
    double dmin0=dmin;

    mint=1.2*tmin;
    maxt=0.8*tmin;

    VPRINT(stdout,"\tEncounter estimated time, tmin = %.6e\n",tmin);
    VPRINT(stdout,"\t\tTesting range = [%.6e,%.6e]\n",mint,maxt);
    VPRINT(stdout,"\tEncounter estimated distance, dmin = %.6e\n",dmin);
    if(!(fabs(tmin0)>0)){
      fprintf(stdout,"\tThis star is too close\n");
      continue;
    }

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
    dvr=atof(fields[Candidates::E_RV]);
    //fprintf(stdout,"dvr = %e\n",dvr);

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
      VPRINT(stdout,"\t\t|%s|\n",vec2strn(cov[i],6,"%-+15.3e"));

    generateMultivariate(cov,mobs,obs,6,Nobs);

    //FIRST ONE IS ALWAYS THE NOMINAL ONE
    obs[0][0]=ra;
    obs[0][1]=dec;
    obs[0][2]=par;
    obs[0][3]=mura;
    obs[0][4]=mudec;
    obs[0][5]=vr;

    VPRINT(stdout,"\tSurrogate random properties:\n");
    for(int i=Nobs;i-->0;){
      VPRINT(stdout,"\t\tObservation %d: %s\n",i,vec2strn(obs[i],6,"%.10e "));
    }

    //COMPUTE THE NOMINAL CONDITIONS
    d=AU/tan(par/(60*60*1000.0)*DEG)/PARSEC;
    VPRINT(stdout,"\t\tDistance: %e\n",d);
    radrec_c(d,ra*DEG,dec*DEG,xg);
    mxv_c(M_J2000_Galactic,xg,x);
    recrad_c(x,&tmp,&l,&b);
    calcUVW(ra,dec,par,dpar,mura,dmura,mudec,dmudec,vr,dvr,x+3,dx+3);
    VPRINT(stdout,"\t\tConditions at present (heliocentric): %s\n",vec2strn(x,6,"%e "));
    vscl_c(PARSEC/1e3,x,x);
    LSR2GC(x,xg);
    vscl_c(1e3/UL,xg,xg);//SET UNITS
    vscl_c(1e3/UV,xg+3,xg+3);
    cart2polar(xg,xInt0,1.0);

    //INTEGRATE BACK TO TIME OF INGRESS
    //*
    hstep=fabs(ting)/10;
    polar2cart(xInt0,x);
    VPRINT(stdout,"\t\tConditions at present: %s\n",vec2strn(xInt0,6,"%e "));
    vscl_c(ting*UT*UV/PARSEC,xg+3,dx);
    vadd_c(xg,dx,xg);
    VPRINT(stdout,"\t\tExpected final conditions (ting = %e, hstep = %e): %s\n",ting,hstep,vec2strn(xg,6,"%e "));
    params[0]=6;

    try{
      integrateEoM(0,xInt0,hstep,2,ting,6,EoMGalactic,params,ts,xInt);
    }catch(int e){
      fprintf(stdout,"\t\t****No suitable integration for star %d***\n",n);
      //getchar();
      continue;
    }
    VPRINT(stdout,"\t\tTimes: %s\n",vec2strn(ts,2,"%e "));
    copyVec(xInt0,xInt[1],6);
    polar2cart(xInt0,xg);
    VPRINT(stdout,"\t\tInitial conditions: %s\n",vec2strn(xg,6,"%e "));
    VPRINT(stdout,"\t\tInitial conditions (polar): %s\n",vec2strn(xInt0,6,"%e "));
    //exit(0);
    //*/

    //INITIAL CONDITION FOR TEST PARTICLES
    copyVec(xnom0,xnoms0,6);

    //DISCRETE METHOD
    try{
      minDistance2(xInt0,xnom0,tmin0,&dmin,&tmin,params);
    }catch(int e){
      fprintf(stdout,"\tNo minimum!\n");
      continue;
    }
    //COMPUTE THE RELATIVE VELOCITY WITH XINT0 AND XNOM0
    nomdmin=dmin;
    nomtmin=tmin;
    polar2cart(xInt0,x);
    polar2cart(xnom0,xg);
    vsubg_c(x,xg,6,dx);
    nomvrel=vnorm_c(dx+3);
    VPRINT(stdout,"\t\tMinimum distance from nominal to nominal (nom. t=%.6e, d=%.6e): t = %.6e, d = %.6e, dv = %.6e\n",tmin0,dmin0,nomtmin,nomdmin,nomvrel*UV/1e3);

    //FILTER OBJECTS TOO FAR FROM THE PRESENT
    if(fabs(tmin)>tRet){
      fprintf(stdout,"\tThis star has gone too far (tmin = %e). Excluding it\n",tmin);
      continue;
    }

    //CALCULATE PROBABILITIES
    Pprob=0;
    Psmed=0;
    Nsur_acc=0;
    nt=0;
    for(int i=0;i<Nsur;i++){
      //if(i<35) continue;
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

      //INTEGRATE BACK TO TIME OF INGRESS
      //*
      hstep=fabs(ting)/10;
      polar2cart(xInt0,x);
      VPRINT(stdout,"\t\tConditions at present: %s\n",vec2strn(xg,6,"%e "));
      vscl_c(ting*UT*UV/PARSEC,xg+3,dx);
      vadd_c(xg,dx,xg);
      VPRINT(stdout,"\t\tExpected final conditions: %s\n",vec2strn(xg,6,"%e "));
      params[0]=6;
      try{
	integrateEoM(0,xInt0,hstep,2,ting,6,EoMGalactic,params,ts,xInt);
      }catch(int e){
	fprintf(stdout,"\t\t****No suitable integration for star %d, surrgate %d***\n",n,i);
	//getchar();
	continue;
      }
      copyVec(xInt0,xInt[1],6);
      polar2cart(xInt0,xg);
      VPRINT(stdout,"\t\tInitial conditions: %s\n",vec2strn(xg,6,"%e "));
      //*/

      VPRINT(stdout,"\t\ttmin0: %.6e\n",tmin0);
      VPRINT(stdout,"\t\tObservations: %s\n",vec2strn(obs[i],6,"%.5e "));
      VPRINT(stdout,"\t\tDistance: %e pc\n",d);
      VPRINT(stdout,"\t\tGalactic coordinates: l = %lf, b = %lf\n",l*RAD,b*RAD);
      VPRINT(stdout,"\t\tInitial position cartesian: %s\n",vec2strn(xg,6,"%.5e "));
      VPRINT(stdout,"\t\tInitial position cylindrical: %s\n",vec2strn(xInt0,6,"%.5e "));

      //CALCULATE MINIMUM DISTANCE AND TIME OF NOMINAL SOLUTION TO SURROGATE
      params[0]=6;

      try{
	minDistance2(xInt0,xnom0,tmin0,&dmin,&tmin,params);
      }catch(int e){
	VPRINT(stdout,"No minimum!\n");
	continue;
      }
      VPRINT(stdout,"\t\tMinimum distance at (nom. d=%.6e, t=%.6e): t = %.6e, d = %.6e\n",tmin0,dmin0,tmin,dmin);
      VPRINT(stdout,"\t\tFinal position star: %s\n",vec2strn(xInt0,6,"%.5e "));
      VPRINT(stdout,"\t\tFinal position nominal: %s\n",vec2strn(xnom0,6,"%.5e "));

      if(fabs(tmin)>tRetMax){
	VPRINT(stdout,"\t\t\tThis surrogate has gone too far. Excluding\n");
	continue;
      }
      if(fabs(tmin)==0){
	VPRINT(stdout,"\t\t\tThis surrogate does not got too far\n");
	continue;
      }

      //PROPAGATE ALL TEST PARTICLE
      params[0]=6*Ntest;
      h=fabs(tmin)/100;
      VPRINT(stdout,"\t\tInitial conditions for all particles: %s\n",vec2strn(xIntp0,6*Ntest,"%.5e "));
      try{
	integrateEoM(0,xIntp0,h,2,tmin,6*Ntest,EoMGalactic,params,ts,xIntp);
      }catch(int e){
	fprintf(stdout,"\t\t****No suitable integration for star %d, surrogate %d***\n",n,i);
	//getchar();
      }
      VPRINT(stdout,"\t\tIntegration result for all particles: %s\n",vec2strn(xIntp[1],6*Ntest,"%.5e "));

      //COMPUTE THE SIZE OF THE CLOUD
      double vradius,rinter,vinter;
      cloudProperties(xIntp[1],Ntest,&hprob,&vradius,&rinter,&vinter);
      sigma=wNormalization2(hprob);
      VPRINT(stdout,"\t\tSize of cloud at = %e Myr\n",tmin);
      VPRINT(stdout,"\t\thprob = %e pc\n",hprob);
      VPRINT(stdout,"\t\tVel.radius = %e km/s\n",vradius*UV/1e3);
      VPRINT(stdout,"\t\tAverage distance = %e pc\n",rinter);
      VPRINT(stdout,"\t\tAverage velocity difference = %e km/s\n",vinter*UV/1e3);
      VPRINT(stdout,"\t\tDensity normalization = %e\n",sigma);
      //VPRINT(fcl,"%e %e %e %e %e\n",tmin,hprob,rinter,vradius*UV/1e3,vinter*UV/1e3);
      //getchar();
      
      //COMPUTE SPH-LIKE PROBABILITY
      double pd,pv;
      Psur=0.0;
      fvel=0.0;
      VPRINT(stdout,"\t\tComparing test particle position with star position\n");
      Pvmed=0.0;
      for(int j=0;j<Ntest;j++){

	polar2cart(xInt0,x);
	polar2cart(xIntp[1]+6*j,xg);
	vsubg_c(x,xg,6,dx);

	D=vnorm_c(dx);
	vrel=vnorm_c(dx+3);

	//FIND MINIMUM APPROXIMATION
	mindmin=MIN(D,mindmin);
	maxdmin=MAX(D,maxdmin);
	mintmin=MIN(tmin,mintmin);
	maxtmin=MAX(tmin,maxtmin);
	minvrel=MIN(vrel,minvrel);
	maxvrel=MAX(vrel,maxvrel);
	dmins[nt]=D;
	tmis[nt]=tmin;
	vrels[nt]=vrel;
	nt++;

	VPRINT(stdout,"\t\t\tDistance to test particle %d (hprob = %e): d=%.6e,vrel=%.6e\n",j,hprob,D,vrel*UV/1e3);

	//CONTRIBUTION TO P FROM DISTANCE
	pd=sigma*wFunction2(D,&hprob)*deltaV;
	VPRINT(stdout,"\t\t\tdV (variable) = %e\n",rinter*rinter*rinter);
	VPRINT(stdout,"\t\t\tdV (fixed) = %e\n",deltaV);

	//CONTRIBUTION TO P FROM VELOCITY
	//Compute the actual scale of the velocity
	ul=1*AU;
	um=0.5*MSUN;
	ut=sqrt(ul*ul*ul/(GCONST*um));
	uv=ul/ut;
	VPRINT(stdout,"\t\t\tUnits of velocity: %e\n",uv);
	vcan=vrel*UV/uv;
	vradiuscan=vradius*UV/uv;
	VPRINT(stdout,"\t\t\tVelocity %e km/s in canonic units: %e\n",vrel*UV/1e3,vcan);
	VPRINT(stdout,"\t\t\tVelocity radius %e km/s in canonic units: %e\n",vradius*UV/1e3,vradiuscan);
	//pv=vinfProbability(vcan-vradius,vcan+vradius,&vpar);
	pv=vinfPosterior(vcan,&vpar)*vradiuscan;
	Pvmed+=pv;
	
	VPRINT(stdout,"\t\t\t\tDistance probability: %.6e\n",pd);
	VPRINT(stdout,"\t\t\t\tVelocity probability: %.6e\n",pv);
	//getchar();
	
	Psur+=pd*pv;
	//COMPUTE CORRECTION FOR RELATIVE STELLAR VELOCITY
      }
      Pvmed/=Ntest;
      Psur/=Ntest;
      Psmed+=Psur;

      //COMPUTE CORRECTION FOR STELLAR DISTANCE
      fdist=(RTRUNC*AU/PARSEC)*(RTRUNC*AU/PARSEC)/(d*d);
      VPRINT(stdout,"\t\tProbability for distance d = %e, RT = %e: %e\n",
	      d,(RTRUNC*AU/PARSEC),fdist);
      Psur*=fdist;
      
      //SURROGATE PROBABILITY
      VPRINT(stdout,"\t\tSurrogate probability: %.6e\n",Psur);

      //ACCUMULATE
      Pprob+=Psur;
      Nsur_acc++;
      //getchar();
    }
    Psmed/=Nsur;
    Pprob/=Nsur;
    fprintf(stdout,"\tNumber of accepted surrogates: %d/%d\n",Nsur_acc,Nsur);
    fprintf(stdout,"\tPsmed = %e, Pvmed = %e, fdist = %e, Pprob = %e\n",Psmed,Pvmed,fdist,Pprob);

    double dminmin,dminmax,dminmed,dminl,dminu;
    double tminmin,tminmax,tminmed,tminl,tminu;
    double vrelmin,vrelmax,vrelmed,vrell,vrelu;
    quantilesVector(dmins,nt,&dminmin,&dminmax,&dminmed,&dminl,&dminu);
    quantilesVector(tmis,nt,&tminmin,&tminmax,&tminmed,&tminl,&tminu);
    quantilesVector(vrels,nt,&vrelmin,&vrelmax,&vrelmed,&vrell,&vrelu);
    fprintf(stdout,"\tMinimum distance: %e\n",dminmin);
    fprintf(stdout,"\tMaximum distance: %e\n",dminmax);
    fprintf(stdout,"\tMedian distance: %e\n",dminmed);
    fprintf(stdout,"\t95%% CL: %e,%e\n",dminl,dminu);
    /*
    FILE*fmin=fopen("dmins.dat","w");
    VPRINT(fmin,"%s\n",vec2strn(dmins,nt,"%e\n"));
    fclose(fmin);
    */

    VPRINT(stdout,"Probability for star %d: %.6e\n",n,Pprob);
    VPRINT(stdout,"Minimum distance (lma.dmin=%e,nom.dmin=%e): min.dmin = %e, max.dmin = %e\n",dmin0,nomdmin,mindmin,maxdmin);

    
    fprintf(fp,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,",
	    Pprob,Psmed,Pvmed,fdist,
	    nomtmin,nomdmin,nomvrel*UV/1e3,
	    tminl,tminmed,tminu,
	    dminl,dminmed,dminu,
	    vrell*UV/1e3,vrelmed*UV/1e3,vrelu*UV/1e3);
    
    fprintf(fp,"%s",values);
    fflush(fp);

    Ntelaps++;
    telaps+=elapsedTime();
    //if(VERBOSE) getchar();
    if(qinterrupt) break;
  }
  telaps/=Ntelaps;
  fprintf(stdout,"Average time per candidate: %f s\n",telaps);
  fclose(fc);
  fclose(fp);
  //fclose(fcl);
  VPRINT(stdout,"DONE.\n");
  
  //REPORT THAT THE COMPUTATION HAS BEEN DONE
  if(ipart<0){
    sprintf(Filename,"rm log/start-%s",WANDERER);
    system(Filename);
    sprintf(Filename,"touch log/done-%s",WANDERER);
    system(Filename);
  }else{
    sprintf(Filename,"rm log/start-%s.%05d",WANDERER,ipart);
    system(Filename);
    sprintf(Filename,"touch log/done-%s.%05d",WANDERER,ipart);
    system(Filename);
  }

  telaps=elapsedTime(0);
  fprintf(stdout,"Total elapsed time = %.5f (%.5f min)\n",telaps,telaps/60.0);
  return 0;
}
