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

//SET TO 1 IF YOU WANT A VERBOSE RUN
#define VERBOSE 2
#define VSTREAM stderr
//#define VSTREAM stdout

int main(int argc,char* argv[])
{
  /*
    Example: ./probability.exe [candidates-Oumuamu.csv 1]

    Where:
      candidadtes-Oumuamua.csv: base name of the candidate file.
      1: number of candidate file (candidadtes-Oumuamua.csv.00001)

    Function: calculate the IOP probability for a list of candidates.

    Input:
    * wanderer.csv
    * candidates.csv

    Output: 

    * progenitors-<wanderer>.csv (alternatively progenitors-<wanderer>.csv.00001)

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
  #include <probability2.conf>
  printHeader(stdout,"COMPUTING PROBABILITIES",'*');

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initWander();

  ////////////////////////////////////////////////////
  //VARIABLES DECLARATION
  ////////////////////////////////////////////////////
  //COUNTERS
  int i,j,n,nt;
  int ip,jp;

  //INPUT
  int ipart=-1;
  char basename[100];
  char suffix[100];

  //GENERAL VARIABLES
  double dmin=0.0,tmin=0.0;
  double tmin0,dmin0;
  double mint=0.0,maxt=0.0;
  double nomtmin=0.0,nomdmin=0.0,nomvrel=0.0;
  double dminmin,dminmax,dminl,dminu,dminmed;

  double ni=0.0,nsum=0.0;
  double P=0.0,Pspeed=0.0,Psum=0.0,Psump=0.0,Psumv=0.0;

  double Pp,Ppos,Pv,Pvel,Pposvel=0.0;
  double Pdist=0.0;
  double IOP=0.0;

  double Pvmed=0.0,Psmed=0.0,Pprob=0.0,Psur=0.0;
  double fvel=0.0,fdist=0.0;
  double hstep;
  double mindmin=1e100,maxdmin=-1e100;
  double minvrel=1e100,maxvrel=-1e100;
  double mintmin=1e100,maxtmin=-1e100;
  double phi,theta;

  //INTEGER
  int nsys=6;
  int nsysp=6*Ntest;
  int nfullp=6*(Ntest+Nsur);
  int Ntimes=100;
  int Ntimesp=10000;

  //PARAMETERS
  double params[10],mparams[23];

  //COUNTERS
  int Ncand=0;
  int Nstar_close=0,Nstar_far=0,Nstar_zero=0;
  int Nstar_prop=0,Nstar_nosur=0,Nstar_nomin=0,Nstar_noint=0,Nstar_nointall=0;
  int Nsur_suc=0,Nsur_acc=0; 

  //VELOCITY VARIABLES
  int nv=0;
  double ul,um,ut,uv,dv;
  double pvs[MAXCOLS],vinfs[MAXCOLS];
  gsl_interp_accel *acc;
  gsl_spline *spline;
  struct vinfpar vpar;
  double r90,v90;

  //Units of ejection-speed program
  ul=1*AU;
  um=0.5*MSUN;
  ut=sqrt(ul*ul*ul/(GCONST*um));
  uv=ul/ut;

  //SWITCHES
  int qinterrupt=0;
  
  //DENSITY VARIABLES
  double hprob=0.5;
  double sigma=wNormalization(hprob);
  double D=0.0;
  double Dmin;
  double deltaOmegav,deltaOmegaing;

  //VECTORS AND MATRICES WITH INTEGRATIONS
  double *xdum;
  double *x=vectorAllocate(6);
  double *xg=vectorAllocate(6);
  double *xc=vectorAllocate(6);
  double *dx=vectorAllocate(6);
  double *mobs=vectorAllocate(6);

  double *xIntp0=vectorAllocate(nsysp);
  double *xFullp0=vectorAllocate(nfullp);
  double *xFullpe=vectorAllocate(nfullp);
  double *xIntc0=vectorAllocate(nsysp);
  double *xInt0=vectorAllocate(nsys);
  double *xnom0=vectorAllocate(nsysp);
  double *xnoms0=vectorAllocate(nsysp);

  double *tsp=vectorAllocate(Ntimesp);
  double *ts=vectorAllocate(Ntimesp);

  double *xTraj0=vectorAllocate(nsys);
  double **xTraj1=matrixAllocate(Ntimesp,nsys);
  double **xTraj2=matrixAllocate(Ntimesp,nsys);
  
  double **xInt=matrixAllocate(Ntimes,nsys);
  double **xIntp=matrixAllocate(Ntimesp,nsysp);
  double **xFullp=matrixAllocate(Ntimesp,nfullp);
  double **xIntc=matrixAllocate(Ntimesp,nsysp);
  double **obs=matrixAllocate(Nsur,6);
  double **cov=matrixAllocate(6,6);
  char **Fields=charMatrixAllocate(MAXCOLS,MAXTEXT);

  double **vrelvec=matrixAllocate(Ntest,3);
  double **v1=matrixAllocate(Ntest,3);
  double *Ds=vectorAllocate(Ntest*Nsur);
  double *tmis=vectorAllocate(Ntest*Nsur);
  double *vrels=vectorAllocate(Ntest*Nsur);

  //STELLAR PROPERTIES
  //     0  1   2   3    4     5 
  double ra,dec,par,mura,mudec,vr;
  double dra,ddec,dpar,dmura,dmudec,dvr;
  double d,l,b,h;
  double UVW[3];
  
  //SPICE VARIABLES
  SpiceDouble M_J2000_Galactic[3][3];
  pxform_c("J2000","GALACTIC",0,M_J2000_Galactic);

  //CLOUD PROPERTIES
  double rinter,vinter,vradius;

  //ENUMERATIONS
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

  //FIELDS OF STELLAR PROPERTIES
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

  //FILES
  FILE *fc;
  FILE *fv;
  FILE *fp;
  FILE *fso;
  FILE *fst;

  //DEPRECATED
  double *dxIntdt,*xp1,*xp2,*xpmin,*x0;
  double Dmax=0,*xt1,*xt2,vrel;
  double tmp,ting;
  char ctmp[100];
  double ftmin;
  double t;
  int it;

  ////////////////////////////////////////////////////
  //LOAD COMMAND LINE OPTIONS
  ////////////////////////////////////////////////////
  ipart=-1;
  //By default basename is candidates-<wanderer>.csv
  sprintf(basename,"%s/candidates-%s.csv",SCR_DIR,WANDERER);
  sprintf(suffix,"%s",WANDERER);
  if(argc>1){
    //The base name is provided
    strcpy(basename,argv[1]);
    //If the number of the file is provided get it
    if(argc>2) ipart=atoi(argv[2]);
    //If not assume 00001
    else ipart=1;
    //Suffix
    sprintf(suffix,"%s.%05d",WANDERER,ipart);
  }
  print0(stdout,"\tAnalysing candidates in '%s'\n",basename);
  if(ipart>0)
    print0(stdout,"\this is the %d part of a parallel analysis\n",ipart);

  ////////////////////////////////////////////////////
  //SIGNAL THE START OF THE COMPUTATION
  ////////////////////////////////////////////////////
  sprintf(Filename,"date +%%s > %s/start-%s",LOG_DIR,suffix);
  system(Filename);

  ////////////////////////////////////////////////////
  //GLOBAL PROPERTIES FOR INTEGRATION IN GALACTIC POT.
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
  printHeader(stdout,"READING SURROGATE OBJECTS",'-');

  //Read file
  sprintf(Filename,"scratch/wanderer-%s.csv",WANDERER);
  fc=fopen(Filename,"r");

  //Get header
  fgets(Line,MAXLINE,fc);//HEADER
  print0(stdout,"\tReading initial conditions of test particles\n");
  print1(VSTREAM,"Reading initial conditions of test particles\n");

  i=0;
  while(fgets(Line,MAXLINE,fc)!=NULL){

    print2(VSTREAM,"\tParticle %d:\n",i);
    //Parse line
    parseLine(Line,Fields,&Nfields);

    //Initial time
    tsp[i]=0.0;

    //Position in the initial conditions vector
    ip=6*i;

    //Ingress time
    ting=atof(Fields[Wanderer::TING])/UT;
    print2(VSTREAM,"\t\tIngress time: %e\n",ting*UT/YEAR);

    //Reading heliocentric ingress position
    j=Wanderer::XECL;
    for(int k=0;k<6;k++) x[k]=atof(Fields[j++]);
    copyVec(vrelvec[i],x,3);
    print2(VSTREAM,"\t\tHeliocentric ingress (km,km/s): %s\n",vec2strn(x,6,"%e "));

    //Reading heliocentric galactic position
    j=Wanderer::XGAL;
    for(int k=0;k<6;k++) x[k]=atof(Fields[j++]);
    print2(VSTREAM,"\t\tHeliocentric galactic (km,km/s): %s\n",vec2strn(x,6,"%e "));
    
    //Converting from heliocentric to galactocentric
    LSR2GC(x,xg);
    print2(VSTREAM,"\t\tGalactocentric (km,km/s): %s\n",vec2strn(xg,6,"%e "));

    //Convert from kilometers to program units
    vscl_c(1e3/UL,xg,xg);
    vscl_c(1e3/UV,xg+3,xg+3);
    print2(VSTREAM,"\t\tGalactocentric (UL,UV): %s\n",vec2strn(xg,6,"%e "));

    //Convert to cylindrical
    cart2polar(xg,xIntp0+ip);
    print2(VSTREAM,"\t\tInitial condition for particle %d: %s\n",i,
	   vec2strn(xIntp0+ip,6,"%e "));

    i++;

    //Break at a given number of test particles (see configuration file)
    if(i==Ntest) break;
  }
  print0(stdout,"\tNtest = %d\n",Ntest);
  fclose(fc);

  //Compute
  deltaOmegaing=solidAngle(vrelvec,Ntest,&vrel,&dv);
  print2(VSTREAM,"\tIngress solid angle:%e\n",deltaOmegaing);
  print2(VSTREAM,"\tAverage distance (AU): %e\n",vrel*1e3/AU);
  print2(VSTREAM,"\tDistance dispesion (AU): %e\n",dv*1e3/AU);

  //Backup nominal initial conditions
  copyVec(xnoms0,xIntp0,6);

  //Initial conditions
  copyVec(xFullp0,xIntp0,nsysp);

  print2(VSTREAM,"Initial condition nominal test particle: %s\n",
	 vec2strn(xIntp0,6,"%e "));

  ////////////////////////////////////////////////////
  //READ POSTERIOR EJECTION VELOCITY DISTRIBUTION
  ////////////////////////////////////////////////////
  printHeader(stdout,"READING VELOCITY DISTRIBUTION",'-');
  fv=fopen("db/ejection-posterior.data","r");
  nv=0;
  while(fgets(Line,MAXLINE,fv)!=NULL){
    parseLine(Line,Fields,&Nfields," ");

    //Velocity
    vinfs[nv]=atof(Fields[0]);

    //Probability
    pvs[nv]=atof(Fields[1]);

    nv++;
  }
  //Interpolating
  acc=gsl_interp_accel_alloc();
  spline=gsl_spline_alloc(gsl_interp_cspline,nv);
  gsl_spline_init(spline,vinfs,pvs,nv);

  //Initializing vpar struct
  vpar.xmin=0.0;
  vpar.xmax=vinfs[nv-1];
  vpar.a=acc;
  vpar.s=spline;

  print2(VSTREAM,"\tP(vinf=%e) = %e\n",0.2,vinfPosterior(0.2,&vpar));
  print0(stdout,"\tIntegral = %e +/- %e\n",vinfProbability(0,10,&vpar));

  ////////////////////////////////////////////////////
  //READING POTENTIAL OBJECTS
  ////////////////////////////////////////////////////
  printHeader(stdout,"READING CANDIDATES FILE AND OPENNING OUTPUT FILE",'-');

  //Input candidates
  if(ipart<0){
    strcpy(Filename,basename);
  }else{
    sprintf(Filename,"%s.%05d",basename,ipart);
  }
  if((fc=fopen(Filename,"r"))==NULL){
    fprintf(stderr,"No candidates file '%s'\n",Filename);
  }
  print0(stdout,"\tReading input filename %s\n",Filename);
  fgets(Line,MAXLINE,fc);

  //Output progenitors
  if(ipart<0){
    sprintf(Filename,"scratch/progenitors-%s.csv",WANDERER);
  }else{
    sprintf(Filename,"scratch/progenitors-%s.csv.%05d",WANDERER,ipart);
  }
  if(qsingle){
    sprintf(Filename,"scratch/progenitors-single.csv",WANDERER);
  }
  print2(VSTREAM,"\tGenerationg output filename %s...\n",Filename);
  fp=fopen(Filename,"w");

  //Header of progenitors
  fprintf(fp,"nomtmin,nomdmin,nomvrel,dminl,dminmed,dminu,");
  //fprintf(fp,"tminl,tminmed,tminu,vrell,vrelmed,vrelu,");
  fprintf(fp,"Ppos,Pvel,Pposvel,Pdist,IOP,");
  fprintf(fp,"%s",Line);

  //If qsingle save all particles position
  if(qsingle){
    if(strlen(hip_single)==0) sprintf(basename,"TYC%s",tyc_single);
    else sprintf(basename,"HIP%s",hip_single);
    sprintf(Filename,"scratch/particles-%s-%s.dat",WANDERER,basename);
    fso=fopen(Filename,"w");
  }

  ////////////////////////////////////////////////////
  //START PROBABILITY CALCULATION
  ////////////////////////////////////////////////////
  printHeader(stdout,"CALCULATING PROBABILITIES",'-');

  //TIME
  Telaps=0.0;
  Nelaps=0.0;
  elapsedTime();

  params[0]=nsys;

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //START READING OF STARS
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  n=0;
  while(fgets(Line,MAXLINE,fc)!=NULL){
    
    //PARSE FIELDS
    strcpy(Values,Line);
    parseLine(Line,Fields,&Nfields);
    n++;

    //LOCAL MINIMA
    mindmin=1e100;maxdmin=-1e100;
    minvrel=1e100;maxvrel=-1e100;
    mintmin=1e100;maxtmin=-1e100;
    
    //IF QSINGLE
    if(qsingle){
      if(strlen(hip_single)==0)
	if(strcmp(Fields[Candidates::TYCHO2_ID],tyc_single)!=0) continue;
	else qinterrupt=1;
      else
	if(strcmp(Fields[Candidates::HIP],hip_single)!=0) continue;
	else qinterrupt=1;
    }
    
    print0(stdout,"Computing probabilities for candidate star %d (%s,%s)...\n",
	    n,Fields[Candidates::HIP],Fields[Candidates::TYCHO2_ID]);
    print1(VSTREAM,"Computing probabilities for candidate star %d (%s,%s)...\n",
	   n,Fields[Candidates::HIP],Fields[Candidates::TYCHO2_ID]);

    //LMA tmin and dmin
    tmin=atof(Fields[Candidates::TMIN]);
    dmin=atof(Fields[Candidates::DMIN]);
    //Preserve variables
    tmin0=tmin;
    dmin0=dmin;

    print1(VSTREAM,"\tEstimated nominal time and distance: dmin=%e, tmin=%e\n",
	    dmin,tmin);

    //If tmin is zero the star is too close
    if(!(fabs(tmin0)>1/*year*/)){
      print0(stdout,"\t***This star is too close***\n");
      print1(VSTREAM,"\t***This star is too close***\n");
      Nstar_close++;
      continue;
    }

    //Range of search
    mint=1.2*tmin;
    maxt=0.8*tmin;
    print2(VSTREAM,"\t\tTesting range = [%.6e,%.6e]\n",mint,maxt);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //STELLAR PROPERTIES
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    mobs[0]=ra=atof(Fields[FC[RA]]);
    dra=atof(Fields[FC[RA_ERROR]])*MAS;

    mobs[1]=dec=atof(Fields[FC[DEC]]);
    ddec=atof(Fields[FC[DEC_ERROR]])*MAS;

    mobs[2]=par=atof(Fields[FC[PARALLAX]]);
    dpar=atof(Fields[FC[PARALLAX_ERROR]]);

    mobs[3]=mura=atof(Fields[FC[PMRA]]);
    dmura=atof(Fields[FC[PMRA_ERROR]]);

    mobs[4]=mudec=atof(Fields[FC[PMDEC]]);
    dmudec=atof(Fields[FC[PMDEC_ERROR]]);

    mobs[5]=vr=atof(Fields[Candidates::RV]);
    dvr=atof(Fields[Candidates::E_RV]);

    //COVARIANCE MATRIX
    /*RA*/cov[0][0]=dra*dra;
    cov[0][1]=atof(Fields[FC[RA_DEC_CORR]])*dra*ddec;
    cov[0][2]=atof(Fields[FC[RA_PARALLAX_CORR]])*dra*dpar;
    cov[0][3]=atof(Fields[FC[RA_PMRA_CORR]])*dra*dmura;
    cov[0][4]=atof(Fields[FC[RA_PMDEC_CORR]])*dra*dmudec;
    cov[0][5]=0.0;
    /*DEC*/cov[1][1]=ddec*ddec;
    cov[1][0]=cov[0][1];
    cov[1][2]=atof(Fields[FC[DEC_PARALLAX_CORR]])*ddec*dpar;
    cov[1][3]=atof(Fields[FC[DEC_PMRA_CORR]])*ddec*dmura;
    cov[1][4]=atof(Fields[FC[DEC_PMDEC_CORR]])*ddec*dmudec;
    cov[1][5]=0.0;
    /*PAR*/cov[2][2]=dpar*dpar;
    cov[2][0]=cov[0][2];
    cov[2][1]=cov[1][2];
    cov[2][3]=atof(Fields[FC[PARALLAX_PMRA_CORR]])*dpar*dmura;
    cov[2][4]=atof(Fields[FC[PARALLAX_PMDEC_CORR]])*dpar*dmudec;
    cov[2][5]=0.0;
    //Distance
    d=starDistance(par);//In parsecs
    /*MURA*/cov[3][3]=dmura*dmura;
    cov[3][0]=cov[0][3];
    cov[3][1]=cov[1][3];
    cov[3][2]=cov[2][3];
    cov[3][4]=atof(Fields[FC[PMRA_PMDEC_CORR]])*dmura*dmudec;
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
    
    //Show properties
    print1(VSTREAM,"\tStellar distance: %e pc\n",d);
    print1(VSTREAM,"\tStellar properties (ra,dec,par,mura,mudec,vr): %s\n",
	    vec2strn(mobs,6,"%.5e "));
    print1(VSTREAM,"\tErrors:");
    for(int i=0;i<6;i++) print1(VSTREAM,"%.6e ",sqrt(cov[i][i]));
    print1(VSTREAM,"\n");
    print1(VSTREAM,"\tStar Covariance Matrix:\n");
    for(int i=0;i<6;i++)
      print1(VSTREAM,"\t\t|%s|\n",vec2strn(cov[i],6,"%-+15.3e"));
    Ncand++;

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //GENERATE SURROGATE STAR PROPERTIES
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    try{
      generateMultivariate(cov,mobs,obs,6,Nsur);
    }catch(int e){
      print0(stdout,"\t***Star %d has a problem in its properties**\n",n);
      print1(VSTREAM,"\t***Star %d has a problem in its properties**\n",n);
      Nstar_prop++;
      continue;
    }
    //Nominal star properties
    obs[0][0]=ra;
    obs[0][1]=dec;
    obs[0][2]=par;
    obs[0][3]=mura;
    obs[0][4]=mudec;
    obs[0][5]=vr;

    //Surrogate properties
    print2(VSTREAM,"\tSurrogate random properties:\n");
    for(int i=Nsur;i-->0;){
      print2(VSTREAM,"\t\tObservation %d: %s\n",i,vec2strn(obs[i],6,"%.10e "));
    }

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INTEGRATE ALL STARS UNTIL TING
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    j=nsysp;
    for(i=0;i<Nsur;i++){
      elapsedTime();
      ip=6*i;

      //Get key properties
      ra=obs[i][0];
      dec=obs[i][1];
      par=obs[i][2];
      mura=obs[i][3];
      mudec=obs[i][4];
      vr=obs[i][5];

      //Convert sky position to cartesian position
      radrec_c(d,ra*DEG,dec*DEG,xg);

      //From sky to galactic
      mxv_c(M_J2000_Galactic,xg,x);

      //Galactic coordinates
      recrad_c(x,&tmp,&l,&b);

      //Velocity
      calcUVW(ra,dec,par,dpar,mura,dmura,mudec,dmudec,vr,dvr,x+3,dx+3);

      //Convert from parsecs to km
      vscl_c(PARSEC/1e3,x,x);

      //Convert to galactic
      LSR2GC(x,xg);

      //Convert from km to system units
      vscl_c(1e3/UL,xg,xg);
      vscl_c(1e3/UV,xg+3,xg+3);

      //Convert to polar 
      cart2polar(xg,xFullp0+nsysp+ip);

      //Integrate until ting
      params[0]=6;
      hstep=fabs(ting)/10;
      try{
	integrateEoM(0,xFullp0+nsysp+ip,hstep,2,ting,6,EoMGalactic,params,ts,xInt);
	copyVec(xFullp0+nsysp+ip,xInt[1],6);
      }catch(int e){
	print0(stdout,"\t\t***No suitable integration for surrogate %d of star %d***\n",i,n);
	print1(VSTREAM,"\t\t***No suitable integration for surrogate %d of star %d***\n",i,n);
	copyVec(xFullpe+nsysp+ip,XFOO,6);
	Nstar_nosur++;
	continue;
      }
      //Actual initial conditions for nominal star
    }

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //MINIMUM DISTANCE AND TIME BETWEEN NOM.
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //Initial condition for test particles
    copyVec(xnom0,xnoms0,6);
    copyVec(xInt0,xFullp0+nsysp,6);
    copyVec(xTraj0,xInt0,6);

    //Minimum distance nominal to nominal
    try{
      minDistance(xInt0,xnom0,tmin0,&dmin,&tmin,params);
      nomdmin=dmin;
      nomtmin=tmin;
    }catch(int e){
      print0(stdout,"\t***No minimum for star %d***\n",n);
      print1(VSTREAM,"\t***No minimum for star %d***\n",n);
      Nstar_nomin++;
      continue;
    }

    //Relative position and velocity at minimum
    polar2cart(xInt0,x);
    polar2cart(xnom0,xg);
    vsubg_c(x,xg,6,dx);
    nomvrel=vnorm_c(dx+3)*UV/1e3;

    print1(VSTREAM,"\tMinimum distance from nominal to nominal (LMA t=%.6e, d=%.6e): t = %.6e, d = %.6e, dv = %.6e\n",tmin0,dmin0,nomtmin,nomdmin,nomvrel);

    //Star has gone too far
    if(fabs(tmin)>tRet){
      print0(stdout,"\t***This star has gone too far (tmin = %e). Excluding it**\n",tmin);
      print1(VSTREAM,"\t***This star has gone too far (tmin = %e). Excluding it**\n",tmin);
      Nstar_far++;
      continue;
    }

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //ATTEMPT AT INTEGRATE SURROGATE STARS
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    print1(VSTREAM,"\tAttempting to integrate %d surrogates\n",Nsur);
    for(i=0;i<Nsur;i++){
      ip=6*i;
      xdum=xFullp0+nsysp+ip;
      if(xdum[0]==99.99){
	copyVec(xFullpe+nsysp+ip,XFOO,6);
	continue;
      }
      print2(VSTREAM,"\t\tAttempting to integrate surrogate %d...\n",i);
      params[0]=6;
      hstep=fabs(nomtmin)/100;
      try{
	integrateEoM(0,xdum,hstep,2,nomtmin,6,EoMGalactic,params,ts,xIntp);
	copyVec(xFullpe+nsysp+ip,xIntp[1],6);
	Nsur_suc++;
      }catch(int e){
	print0(stdout,"\t\t***No integration for surrogate %d of star %d***\n",i,n);
	print1(VSTREAM,"\t\t***No integration for surrogate %d of star %d***\n",i,n);
	copyVec(xFullpe+nsysp+ip,XFOO,6);
	Nstar_noint++;
	continue;
      }
    }
    print1(VSTREAM,"\t\tSuccessful integration for %d surrogates\n",Nsur_suc);
    
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INTEGRATE TEST PARTICLES UNTIL TMIN
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    print1(VSTREAM,"\tAttempting to integrate %d particles\n",Ntest);
    params[0]=nsysp;
    hstep=fabs(nomtmin)/100;
    try{
      integrateEoM(0,xFullp0,hstep,2,nomtmin,nsysp,EoMGalactic,params,ts,xIntp);
      copyVec(xFullpe,xIntp[1],nsysp);
    }catch(int e){
      print0(stdout,"\t\t***No integration for all particles and star %d***\n",n);
      print1(VSTREAM,"\t\t***No integration for all particles and star %d***\n",n);
      Nstar_nointall++;
      continue;
    }
    if(qsingle){
      fprintf(fso,"%s\n",vec2strn(xFullpe,nfullp,"%e "));
      fflush(fso);
    }
    print1(VSTREAM,"\t\tSuccessful integration of %d test particles\n",Ntest);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //CLOUD PARTICLES PROPERTIES
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    cloudProperties2(xFullpe,Ntest,&r90,&v90);
    print1(VSTREAM,"\tSize of the cloud of objects at = %e yr\n",tmin);
    print1(VSTREAM,"\t\tr90 = %e pc\n",r90);
    print1(VSTREAM,"\t\tv90 = %e km/s\n",v90*UV/1e3);

    cloudProperties2(xFullpe+nsysp,Nsur,&r90,&v90);
    print1(VSTREAM,"\tSize of the cloud of stars at = %e yr\n",tmin);
    print1(VSTREAM,"\t\tr90 = %e pc\n",r90);
    print1(VSTREAM,"\t\tv90 = %e km/s\n",v90*UV/1e3);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //PROBABILITY CALCULATION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //Loop on test particle
    nt=0;
    print1(VSTREAM,"\tAttempting probability calculation\n");
    for(j=0;j<Ntest;j++){
      jp=6*j;
      print2(VSTREAM,"\tTest particle %d:\n",j);

      //Position particle
      xdum=xFullpe+jp;

      //Cartesian position
      polar2cart(xdum,xg);

      //Loop on surrogate star
      nsum=0;
      for(i=0;i<Nsur;i++){
	ip=6*i;
	print2(VSTREAM,"\t\tSurrogate star %d:\n",i);

	//Position star
	xdum=xFullpe+nsysp+ip;
	if(xdum[0]==99.99){
	  print2(VSTREAM,"\t\t\tSkipping\n",i);
	  continue;
	}
	//Cartesian position
	polar2cart(xdum,xc);
	
	//Calculate the distance from object to surrogate star
	vsubg_c(xg,xc,6,dx);
	D=vnorm_c(dx);
	Ds[nt]=D;
	print2(VSTREAM,"\t\t\tDistance = %e\n",D);

	//Number density
	ni=wFunction2(D,&hprob);
	print2(VSTREAM,"\t\t\tNumber density contribution = %e\n",ni);
	nsum+=ni;

	nt++;
      }
      print2(VSTREAM,"\t\tSum ni = %.17e\n",nsum);
      if(nsum==0){
	print0(stdout,"\t\t***This star has zero probability***\n");
	print1(VSTREAM,"\t\t***This star has zero probability***\n");
	Nstar_zero++;
	break;
      }

      //Relative velocity with respect to nominal star
      xdum=xFullpe+nsysp;
      polar2cart(xdum,xc);
      vsubg_c(xc,xg,6,dx);
      copyVec(vrelvec[j],dx+3,3);

      //Convert velocity to km/s
      vscl_c(UV/1e3,vrelvec[j],vrelvec[j]);
      print2(VSTREAM,"\t\tRelative velocity = %s\n",vec2strn(vrelvec[j],3,"%e "));

      //Relative speed
      vrel=vnorm_c(vrelvec[j]);
      print2(VSTREAM,"\t\tSpeed (uv): %e\n",vrel*1e3/uv);
      
      //Speed probability
      Pspeed=vinfPosterior(vrel*1e3/uv,&vpar);
      print2(VSTREAM,"\t\tPspeed (uv^-1): %e\n",Pspeed);
      Pspeed*=(1e3/uv);
      print2(VSTREAM,"\t\tPspeed (1/(m/s)): %e\n",Pspeed);

      //Probability contribution
      Pp=nsum;
      Pv=vrel*vrel*Pspeed;
      P=Pp*Pv;
      print2(VSTREAM,"\t\tPosition contribution: %e\n",Pp);
      print2(VSTREAM,"\t\tVelocity contribution: %e\n",Pv);
      print2(VSTREAM,"\t\tPosition-velocity contribution: %e\n",P);

      Psump+=Pp;
      Psumv+=Pv;
      Psum+=P;
    }
    if(nsum==0) continue;
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //ENCOUNTER STATISTICS
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    getPercentile(Ds,nt,0.95,
		  &dminmin,&dminmax,
		  &dminl,&dminmed,&dminu);

    print1(VSTREAM,"\tDistance:\n");
    print1(VSTREAM,"\t\tNominal: %e\n",nomdmin);
    print1(VSTREAM,"\t\tMinimum: %e\n",dminmin);
    print1(VSTREAM,"\t\tMaximum: %e\n",dminmax);
    print1(VSTREAM,"\t\tMedian: %e\n",dminmed);
    print1(VSTREAM,"\t\t95%% CL: %e,%e\n",dminl,dminu);

    //Solid angle subtended by the relative velocity of the objects cloud
    deltaOmegav=solidAngle(vrelvec,Ntest,&vrel,&dv);
    print1(VSTREAM,"\tSolid angle relative velocity:%e\n",deltaOmegav);
    print1(VSTREAM,"\tAverage relative speed (km/s): %e\n",vrel);
    print1(VSTREAM,"\tRelative speed dispesion (km/s): %e\n",dv);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //POSITION-VELOCITY PROBABILITY
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    print1(VSTREAM,"\tProbabilities:\n");

    //Position contribution
    Ppos=log10(sigma*deltaV*Psump/(Nsur_suc*Ntest));
    print1(VSTREAM,"\t\tProbability position: %f\n",Ppos);
    
    //Velocity contribution
    Pvel=log10(deltaOmegav*Psumv/(Nsur_suc*Ntest));
    print1(VSTREAM,"\t\tProbability velocity: %f\n",Pvel);

    //Combined
    Pposvel=log10(sigma*deltaV*dv*deltaOmegav/(Nsur_suc*Ntest)*Psum);
    print1(VSTREAM,"\t\tProbability position-velocity: %f\n",Pposvel);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //DISTANCE PROBABILITY
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //Distance probability
    Pdist=log10(deltaOmegaing*(RTRUNC*AU/PARSEC/d)*(RTRUNC*AU/PARSEC/d)/(4*M_PI));
    print1(VSTREAM,"\t\tProbability by distance: %f\n",Pdist);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //TOTAL
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    IOP=Pposvel+Pdist;
    sprintf(Filename,"IOP: %f",IOP);
    printHeader(stdout,Filename,'&',1,30);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //SAVE RESULTS
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    fprintf(fp,"%e,%e,%e,%e,%e,%e,",nomtmin,nomdmin,nomvrel,dminl,dminmed,dminu);
    fprintf(fp,"%e,%e,%e,%e,%e,",Ppos,Pvel,Pposvel,Pdist,IOP);
    fprintf(fp,"%s",Values);
    fflush(fp);

    Telaps+=elapsedTime();
    Nelaps++;
    //if(Nelaps>10) break;
  }
  Telaps/=Nelaps;
  print0(stdout,"Average time per candidate: %f s\n",Telaps);

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
  
  print0(stdout,"Counts:\n");
  print0(stdout,"\tNumber of candidates: %d\n",Ncand);
  print0(stdout,"\tNumber of potential progenitors: %d\n",Nelaps);
  print0(stdout,"\tNumber of stars rejected by zero probability: %d\n",Nstar_zero);
  print0(stdout,"\tNumber of stars rejected by being close: %d\n",Nstar_close);
  print0(stdout,"\tNumber of stars rejected by being far: %d\n",Nstar_far);
  print0(stdout,"\tNumber of stars with problems in properties: %d\n",Nstar_prop);
  print0(stdout,"\tNumber of surrogates not integrables: %d\n",Nstar_nosur);
  print0(stdout,"\tNumber of stars with no minimum: %d\n",Nstar_nomin);
  print0(stdout,"\tNumber of stars with no objects integration: %d\n",Nstar_nointall);

  Telaps=elapsedTime(0);
  print0(stdout,"Total elapsed time = %.5f (%.5f min)\n",Telaps,Telaps/60.0);
  print0(stdout,"DONE.\n");

  fclose(fc);
  fclose(fp);
  if(qsingle) fclose(fso);
  return 0;
}
