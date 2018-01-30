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
# Calculate for a single object
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

  if(x<pars->xmin) x=pars->xmin;
  if(x>pars->xmax) x=pars->xmax;

  double p=gsl_spline_eval(pars->s,x,pars->a);
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
    * potential.csv

    Output: 

  */

  ////////////////////////////////////////////////////
  //CONFIGURATION
  ////////////////////////////////////////////////////
  #include <iwander.conf>
  #include <single.conf>

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initWander();
  
  ////////////////////////////////////////////////////
  //GENERAL VARIABLES
  ////////////////////////////////////////////////////
  double D,Dmax=0,*xt1,*xt2,vrel;
  double tmp,ting;
  char ctmp[100],line[MAXLINE],values[MAXLINE];
  int Ntimes,Ntimesp,Nobs,nsys,nsysp; 
  int ip,n;
  int nfields;
  double params[10],mparams[23];
  double hstep;
  //MATRICES WITH INTEGRATIONS
  double *xnom0,*xnoms0,*xIntp0,*xIntc0,**xIntp,**xFullp,**xFullc,**xIntc;
  double *xInt0,**xInt,*xTraj0,**xTraj1,**xTraj2;
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
  Nobs=Nsur;

  x=(double*)malloc(6*sizeof(double));//LSR STATE VECTOR
  xg=(double*)malloc(6*sizeof(double));//GC STATE VECTOR
  dx=(double*)malloc(6*sizeof(double));//LSR STATE VECTOR
  xpmin=(double*)malloc(6*sizeof(double));//GC STATE VECTOR

  xIntp0=(double*)malloc(nsysp*sizeof(double));
  xIntc0=(double*)malloc(nsysp*sizeof(double));
  xInt0=(double*)malloc(nsys*sizeof(double));
  xTraj0=(double*)malloc(nsys*sizeof(double));
  xnom0=(double*)malloc(nsysp*sizeof(double));
  xnoms0=(double*)malloc(nsysp*sizeof(double));
  
  xInt=(double**)malloc(Ntimes*sizeof(double*));
  for(int j=0;j<Ntimes;j++) xInt[j]=(double*)malloc(nsys*sizeof(double));

  xTraj1=(double**)malloc(Ntimesp*sizeof(double*));
  for(int j=0;j<Ntimesp;j++) xTraj1[j]=(double*)malloc(nsys*sizeof(double));
  xTraj2=(double**)malloc(Ntimesp*sizeof(double*));
  for(int j=0;j<Ntimesp;j++) xTraj2[j]=(double*)malloc(nsys*sizeof(double));

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

    ting=atof(fields[Wanderer::TING])/UT;
      
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
  fprintf(stdout,"Time of ingress: %e\n",ting);

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
  fprintf(stdout,"P(vinf=%e) = %e\n",0.2,vinfPosterior(0.2,&vpar));
  fprintf(stdout,"Integral = %e +/- %e\n",vinfProbability(0,10,&vpar));

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

  //SAVE SURROGATES
  FILE *fss=fopen("surrogate.csv","w");
  fprintf(fss,"ra,dec,par,mura,mudec,vr,d,dummy\n");

  //SAVE CLOUD AND SURROGATES
  FILE *fso=fopen("surrogatenom.dat","w");

  //SAVE TIMES AND DISTANCES AT MINIMA 
  FILE *fsm=fopen("surrogatemin.csv","w");
  fprintf(fsm,"tmin,dmin,vrel\n");

  FILE *fst;

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

    /*
    if(strcmp(fields[Candidates::HIP],hip_single)!=0) continue;
    else qinterrupt=1;
    //*/
    //*
    if(strcmp(fields[Candidates::TYCHO2_ID],tyc_single)!=0) continue;
    else qinterrupt=1;
    //*/

    fprintf(stdout,"Star %d,%s,%s:\n",n,fields[Candidates::HIP],fields[Candidates::TYCHO2_ID]);

    //ESTIMATED TIME OF ENCOUNTER
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

    //GENERATE THE REST
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
    radrec_c(d,ra*DEG,dec*DEG,xg);
    mxv_c(M_J2000_Galactic,xg,x);
    recrad_c(x,&tmp,&l,&b);
    calcUVW(ra,dec,par,dpar,mura,dmura,mudec,dmudec,vr,dvr,x+3,dx+3);
    vscl_c(PARSEC/1e3,x,x);
    LSR2GC(x,xg);
    vscl_c(1e3/UL,xg,xg);//SET UNITS
    vscl_c(1e3/UV,xg+3,xg+3);
    cart2polar(xg,xInt0,1.0);

    //INTEGRATE STAR BACK TO TIME OF INGRESS
    //*
    hstep=fabs(ting)/10;
    polar2cart(xInt0,x);
    fprintf(stdout,"\t\tConditions at present: %s\n",vec2strn(xg,6,"%e "));
    vscl_c(ting*UT*UV/PARSEC,xg+3,dx);
    vadd_c(xg,dx,xg);
    fprintf(stdout,"\t\tExpected final conditions: %s\n",vec2strn(xg,6,"%e "));
    params[0]=6;
    integrateEoM(0,xInt0,hstep,2,ting,6,EoMGalactic,params,ts,xInt);
    copyVec(xInt0,xInt[1],6);
    polar2cart(xInt0,xg);
    fprintf(stdout,"\t\tInitial conditions: %s\n",vec2strn(xg,6,"%e "));
    //*/

    //INITIAL CONDITION FOR TEST PARTICLES
    copyVec(xnom0,xnoms0,6);

    params[0]=6;
    try{
      minDistance2(xInt0,xnom0,tmin0,&dmin,&tmin,params);
    }catch(int e){
      fprintf(stdout,"¡No minimum!\n");
      //continue;
    }
    fprintf(stdout,"\t\tMinimum distance from nominal to nominal (nom. t=%.6e, d=%.6e): t = %.6e, d = %.6e\n",tmin0,dmin0,tmin,dmin);
    nomdmin=dmin;
    nomtmin=tmin;

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
      
      //INTEGRATE BACK TO TIME OF INGRESS
      //*
      hstep=fabs(ting)/10;
      polar2cart(xInt0,x);
      fprintf(stdout,"\t\tConditions at present: %s\n",vec2strn(xg,6,"%e "));
      vscl_c(ting*UT*UV/PARSEC,xg+3,dx);
      vadd_c(xg,dx,xg);
      fprintf(stdout,"\t\tExpected final conditions: %s\n",vec2strn(xg,6,"%e "));
      params[0]=6;
      integrateEoM(0,xInt0,hstep,2,ting,6,EoMGalactic,params,ts,xInt);
      copyVec(xInt0,xInt[1],6);
      polar2cart(xInt0,xg);
      fprintf(stdout,"\t\tInitial conditions: %s\n",vec2strn(xg,6,"%e "));
      //*/

      VPRINT(stdout,"\t\ttmin0: %.6e\n",tmin0);
      VPRINT(stdout,"\t\tObservations: %s\n",vec2strn(obs[i],6,"%.5e "));
      VPRINT(stdout,"\t\tDistance: %e pc\n",d);
      VPRINT(stdout,"\t\tGalactic coordinates: l = %lf, b = %lf\n",l*RAD,b*RAD);
      VPRINT(stdout,"\t\tInitial position cartesian: %s\n",vec2strn(xg,6,"%.5e "));
      VPRINT(stdout,"\t\tInitial position cylindrical: %s\n",vec2strn(xInt0,6,"%.5e "));

      //SAVE SURROGATE
      fprintf(fss,"%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,",
	      ra,dec,par,mura,mudec,vr,d);
      
      //CALCULATE MINIMUM DISTANCE AND TIME OF NOMINAL SOLUTION TO SURROGATE
      params[0]=6;

      copyVec(xTraj0,xInt0,6);
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

      //COMPUTE THE SIZE OF THE CLOUD
      double vradius,rinter,vinter;
      cloudProperties(xIntp[1],Npart,&hprob,&vradius,&rinter,&vinter);
      sigma=wNormalization2(hprob);
      fprintf(stdout,"\t\thprob = %e pc\n",hprob);
      fprintf(stdout,"\t\tVel.radius = %e km/s\n",vradius*UV/1e3);
      fprintf(stdout,"\t\tAverage distance = %e pc\n",rinter);
      fprintf(stdout,"\t\tAverage velocity difference = %e km/s\n",vinter*UV/1e3);
      fprintf(stdout,"\t\tDensity normalization = %e\n",sigma);
      //getchar();

      //SAVE POSITIONS
      fprintf(fso,"%s",vec2strn(xInt0,6,"%e "));
      fprintf(fso,"%s\n",vec2strn(xIntp[1],6*Npart,"%e "));

      //COMPUTE SPH-LIKE PROBABILITY
      double pd,pv;
      Psur=0.0;
      fvel=0.0;
      fprintf(stdout,"\t\tComparing test particle position with star position\n");
      for(int j=0;j<Npart;j++){
	
	polar2cart(xInt0,x);
	polar2cart(xIntp[1]+6*j,xg);
	vsubg_c(x,xg,6,dx);

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

	fprintf(stdout,"\t\t\tDistance to test particle %d (hprob = %e): d=%.6e,vrel=%.6e\n",j,hprob,D,vrel*UV/1e3);
	//CONTRIBUTION TO P FROM DISTANCE
	//pd=sigma*wFunction2(D,&hprob)*rinter*rinter*rinter;
	pd=sigma*wFunction2(D,&hprob)*deltaV;
	fprintf(stdout,"\t\t\tdV (variable) = %e\n",rinter*rinter*rinter);
	fprintf(stdout,"\t\t\tdV (fixed) = %e\n",deltaV);

	//CONTRIBUTION TO P FROM VELOCITY
	//Compute the actual scale of the velocity
	ul=1*AU;
	um=0.5*MSUN;
	ut=sqrt(ul*ul*ul/(GCONST*um));
	uv=ul/ut;
	fprintf(stdout,"\t\t\tUnits of velocity: %e\n",uv);
	vcan=vrel*UV/uv;
	vradiuscan=vradius*UV/uv;
	fprintf(stdout,"\t\t\tVelocity %e km/s in canonic units: %e\n",vrel*UV/1e3,vcan);
	fprintf(stdout,"\t\t\tVelocity radius %e km/s in canonic units: %e\n",vradius*UV/1e3,vradiuscan);
	//pv=vinfProbability(vcan-vinter*UV/uv,vcan+vinter*UV/uv,&vpar);
	pv=vinfPosterior(vcan,&vpar)*vradiuscan;

	fprintf(stdout,"\t\t\t\tDistance probability: %.6e\n",pd);
	fprintf(stdout,"\t\t\t\tVelocity probability: %.6e\n",pv);
	fprintf(fsm,"%e,%e,%e\n",tmin,D,vrel*UV/1e3);
	//getchar();

	Psur+=pd*pv;
	//COMPUTE CORRECTION FOR RELATIVE STELLAR VELOCITY
      }
      Psur/=Npart;

      //SAVE TRAJECTORIES
      char fname[100];
      sprintf(fname,"scratch/surrogatetraj-%02d.dat",i);
      fst=fopen(fname,"w");

      //NOMINAL PARTICLE TRAJECTORY
      params[0]=6;
      h=fabs(tmin)/100;
      Ntimesp=1000;
      integrateEoM(0,xIntp0,h,Ntimesp,2*tmin,6,EoMGalactic,params,ts,xTraj1);
      integrateEoM(0,xTraj0,h,Ntimesp,2*tmin,6,EoMGalactic,params,ts,xTraj2);
      for(int it=0;it<Ntimesp;it++){
	//PRECISE
	polar2cart(xTraj1[it],x);
	polar2cart(xTraj2[it],xg);
	vsub_c(xg,x,dx);
	D=vnorm_c(dx);
	fprintf(fst,"%e %e ",ts[it],D);
	fprintf(fst,"%s",vec2strn(x,6,"%e "));
	fprintf(fst,"%s",vec2strn(xg,6,"%e "));
	//USING LMAX
	polar2cart(xIntp0,x);
	polar2cart(xTraj0,xg);
	vscl_c(UV*ts[it]*YEAR/PARSEC,x+3,dx);
	vadd_c(x,dx,x);
	fprintf(fst,"%s",vec2strn(x,3,"%e "));
	vscl_c(UV*ts[it]*YEAR/PARSEC,xg+3,dx);
	vadd_c(xg,dx,xg);
	fprintf(fst,"%s",vec2strn(xg,3,"%e "));
	//DISTANCE BY LMA
	vsub_c(xg,x,dx);
	D=vnorm_c(dx);
	fprintf(fst,"%e ",D);
	fprintf(fst,"\n");
      }
      fclose(fst);

      //COMPUTE CORRECTION FOR STELLAR DISTANCE
      fdist=(RTRUNC*AU/PARSEC)*(RTRUNC*AU/PARSEC)/(d*d);
      fprintf(stdout,"\t\tProbability for distance d = %e, RT = %e: %e\n",
	      d,(RTRUNC*AU/PARSEC),fdist);
      Psur*=fdist;
      
      //SURROGATE PROBABILITY
      fprintf(stdout,"\t\tSurrogate probability: %.6e\n",Psur);
      //getchar();
      
      //ACCUMULATE
      Pprob+=Psur;
      //getchar();
      fprintf(fss,"\n");
    }
    Pprob/=Nsur;
    fprintf(stdout,"Probability for star %d: %.6e\n",n,Pprob);
    fprintf(stdout,"Minimum distance (lma.dmin=%e,nom.dmin=%e): min.dmin = %e, max.dmin = %e\n",dmin0,nomdmin,mindmin,maxdmin);

    fprintf(fp,"%e,%e,%e,%e,%e,%e,%e,",Pprob,nomdmin,nomtmin,mindmin,maxdmin,velrelmin*UV/1e3,velrelmax*UV/1e3);
    fprintf(fp,"%s",values);
    fflush(fp);

    if(VERBOSE) getchar();
    if(qinterrupt) break;
  }
  fclose(fc);
  fclose(fp);
  fclose(fss);
  fclose(fso);
  fclose(fsm);

  return 0;
}
