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
# Wanderer dynamics				
#########################################################
*/
#include <iwander.cpp>
using namespace std;

#define VERBOSE 0

int main(int argc,char* argv[])
{
  /*
    Example: ./wanderer.exe
    
    Function: 

    This program perform three different tasks:

    1) Calculate the time t_asymp when the single conic approximation is
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
	  0:NUmber of the object (0 for nominal trajectory) 
	  1-6:Initial random elements, q,e,i,W,w,Mo,to,mu
	  7-12:Asymptotic elements, q,e,i,W,w,Mo,to,mu
	  13:Time of ingress to Solar System
	  14-19:Cartesian position at ingress wrt. Ecliptic J2000
	  20-25:Cartesian position at ingress wrt. J2000
	  26-31:Cartesian position at ingress wrt. Galactic
	  32:Radiant at ingress RA(h) 
	  33:Radiant at ingress DEC (deg)
	  34:Radiant at ingress l(deg)
	  35:Radiant at ingress b(deg)

     * ingress.dat: a summary of the ingress orbit properties
       including the epoch of asymptotic elements and their covariance
       matrix, the time of ingress, the radiant and velocity at
       ingress.

  */

  ////////////////////////////////////////////////////
  //CONFIGURATION
  ////////////////////////////////////////////////////
  #include <iwander.conf>
  #include <wanderer.conf>

  ////////////////////////////////////////////////////
  //INITIALIZE iWANDER
  ////////////////////////////////////////////////////
  initWander();
  
  ////////////////////////////////////////////////////
  //VARIABLES
  ////////////////////////////////////////////////////
  double dt;
  double tp;
  //ELEMENTS
  double q,e,inc,W,w,Mo;
  //POSITION
  SpiceDouble position[6],sun[6],ltmp;
  SpiceDouble M_Eclip_J2000[3][3],M_Eclip_Galactic[3][3];
  SpiceDouble posJ2000[6],RA,DEC,d;
  SpiceDouble posGalactic[6],l,b;
  SpiceDouble posFutureJ2000[6],posFuture[6],posFutureGalactic[6];
  double RAfut,DECfut,lfut,bfut,dfut;
  double h,t_start,t_step,tend,t_stop,t;
  SpiceChar dend[100];
  double h_used,h_next,h_adjust,deltat;
  int i,status;
  //INITIAL CONDITIONS
  double X0[6],X[6],Xu[6],as,tdyn;
  //ELEMENTS
  double elts[8];
  int Ndisp=50,Ncov=500;
  double ds,dp,Xref[6],Xdif[6],dpasymp,elemasymp[8];
  double tdur,dtdur,X0s[6],Xend[6],tfut,durasymp,dasymp,tasymp,dtasymp;
  double during,durold,ting,ding,dold,vasymp;
  char filename[100];

  ////////////////////////////////////////////////////
  //VARIABLES BASED ON CONFIGURATION
  ////////////////////////////////////////////////////
  //GRAVITATIONAL PARAMETER
  double n=ini_n*DEG/DAY;
  double a=ini_a*(-AU/1E3);
  double mu=n*n*a*a*a;
  double munom=mu;
  VPRINT(stdout,"MU Nominal=%.17e\n",munom);

  //EPHEMERIS TIME
  double to=unitim_c(ini_to_jed,"JDTDB","TDB");

  //INTEGRATION
  int npoints=2;
  double tini=to;
  double direction=duration/fabs(duration);
  double params[]={6};
  duration/=UT;
  SpiceDouble elements[8];
  gsl_vector* ielements=gsl_vector_alloc(6);
  gsl_vector_set(ielements,0,ini_e);
  gsl_vector_set(ielements,1,ini_q);
  gsl_vector_set(ielements,2,ini_tp);
  gsl_vector_set(ielements,3,ini_W);
  gsl_vector_set(ielements,4,ini_w);
  gsl_vector_set(ielements,5,ini_i);
  gsl_vector* relements=gsl_vector_alloc(6);
  gsl_matrix* Lo=gsl_matrix_alloc(6,6);
  gsl_matrix* L=gsl_matrix_alloc(6,6);
  for(int i=0;i<6;i++) for(int j=0;j<6;j++) gsl_matrix_set(Lo,i,j,ini_cov[i][j]);
  //DIAGONAL COVARIANCE
  gsl_matrix* D=gsl_matrix_alloc(6,6);
  gsl_matrix_set_zero(D);
  gsl_matrix_set(D,0,0,ini_de*ini_de);
  gsl_matrix_set(D,1,1,ini_dq*ini_dq);
  gsl_matrix_set(D,2,2,ini_cov[2][2]);
  gsl_matrix_set(D,3,3,ini_dW*ini_dW);
  gsl_matrix_set(D,4,4,ini_dw*ini_dw);
  gsl_matrix_set(D,5,5,ini_di*ini_di);
  double ts[2],**Xout=matrixAllocate(2,6);

  ////////////////////////////////////////////////////
  //TRANSFORM MATRICES
  ////////////////////////////////////////////////////
  pxform_c("ECLIPJ2000","GALACTIC",t,M_Eclip_Galactic);
  pxform_c("ECLIPJ2000","J2000",0,M_Eclip_J2000);
  sprintf(Filename,"scratch/ingress-%s.dat",WANDERER);
  FILE* fi=fopen(Filename,"w");
  fprintf(fi,"Epoch=%.1f\n",ini_to_jed);
  fprintf(fi,"Epoch_Date='%s'\n",ini_to_date);
  fprintf(fi,"Reference='%s'\n",reference_solution);
  fprintf(fi,"truncation=%.1e\n",truncation/AU);
  fprintf(fi,"q_nom=%.17e\n",ini_q);
  fprintf(fi,"e_nom=%.17e\n",ini_e);
  fprintf(fi,"i_nom=%.17e\n",ini_i);
  fprintf(fi,"W_nom=%.17e\n",ini_W);
  fprintf(fi,"w_nom=%.17e\n",ini_w);
  fprintf(fi,"M_nom=%.17e\n",ini_M);
  fprintf(fi,"mu_nom=%.17e\n",mu);
  
  ////////////////////////////////////////////////////
  //DISPERSION OF SURROGATE OBJECTS IN THE FUTURE
  ////////////////////////////////////////////////////
  printHeader(stdout,"COMPUTING ERROR IN ELEMENTS");
  for(int j=0;j<=Ndisp+1;j++){
    if(j>0){
      if(qdiagonal) gsl_matrix_memcpy(L,D);
      else gsl_matrix_memcpy(L,Lo);
      gsl_linalg_cholesky_decomp1(L);
      gsl_ran_multivariate_gaussian(RAND,ielements,L,relements);
      n=ini_n+gsl_ran_gaussian(RAND,ini_dn);
      tp=gsl_vector_get(relements,2);
      Mo=n*(ini_to_jed-tp);
      /*q=*/elements[0]=gsl_vector_get(relements,1)*AU/1e3;
      /*e=*/elements[1]=gsl_vector_get(relements,0);
      /*i=*/elements[2]=gsl_vector_get(relements,5)*DEG;
      /*W=*/elements[3]=gsl_vector_get(relements,3)*DEG;
      /*w=*/elements[4]=gsl_vector_get(relements,4)*DEG;
      /*M=*/elements[5]=Mo*DEG;
      /*to=*/elements[6]=to;
      /*mu=*/elements[7]=mu;
    }else{
      /*q=*/elements[0]=ini_q*AU/1e3;
      /*e=*/elements[1]=ini_e;
      /*i=*/elements[2]=ini_i*DEG;
      /*W=*/elements[3]=ini_W*DEG;
      /*w=*/elements[4]=ini_w*DEG;
      /*M=*/elements[5]=ini_M*DEG;
      /*to=*/elements[6]=to;
      /*mu=*/elements[7]=mu;
    }
    conics_c(elements,to,position);
    spkezr_c("SUN",to,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
    vaddg_c(position,sun,6,position);
    vpack_c(position[0]*1E3/UL,position[1]*1E3/UL,position[2]*1E3/UL,X0);
    vpack_c(position[3]*1E3/UV,position[4]*1E3/UV,position[5]*1E3/UV,X0+3);
    as=vnorm_c(X0);
    tdyn=2*M_PI*sqrt(as*as*as/(GGLOBAL*MSUN/UM));
    h=tdyn/1000.0;
    integrateEoM(tini/UT,X0,h,npoints,duration,6,EoM,params,ts,Xout);
    vscl_c(UL/1E3,Xout[1],Xu);vscl_c(UV/1E3,Xout[1]+3,Xu+3);
    if(j==0){
      copyVec(Xref,Xu,6);
    }else{
      vsub_c(Xref,Xu,Xdif);
      ds+=vnorm_c(Xdif);
    }
  }
  ds/=Ndisp;
  ds*=1e3/AU;//AVERAGE DISTANCE IN AU
  fprintf(stdout,"Typical size of the cloud at (t = %e, d = %e): %e AU\n",duration*UT/YEAR,vnorm_c(Xu)*1e3/AU,2*ds);
  fprintf(fi,"t_test=%.17e\n",duration*UT/YEAR);
  fprintf(fi,"d_test=%.17e\n",vnorm_c(Xu)*1e3/AU);
  fprintf(fi,"dr_test=%.17e\n",2*ds);

  //DATE
  tend=tini+duration*UT;
  deltet_c(tend,"ET",&deltat);
  tend+=deltat;
  et2utc_c(tend,"C",2,100,dend);
  VPRINT(stdout,"Position of nominal object at t = %e yr (%.17e, %s): %s\n",duration*UT/YEAR,tend,dend,vec2strn(Xref,6,"%.17e "));
  
  ////////////////////////////////////////////////////
  //CALCULATE TIME OF ASYMPTOTIC ELEMENTS
  ////////////////////////////////////////////////////
  /*q=*/elements[0]=ini_q*AU/1e3;
  /*e=*/elements[1]=ini_e;
  /*i=*/elements[2]=ini_i*DEG;
  /*W=*/elements[3]=ini_W*DEG;
  /*w=*/elements[4]=ini_w*DEG;
  /*M=*/elements[5]=ini_M*DEG;
  /*to=*/elements[6]=to;
  /*mu=*/elements[7]=mu;
  conics_c(elements,to,position);
  spkezr_c("SUN",to,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
  vaddg_c(position,sun,6,position);
  vpack_c(position[0]*1E3/UL,position[1]*1E3/UL,position[2]*1E3/UL,X0);
  vpack_c(position[3]*1E3/UV,position[4]*1E3/UV,position[5]*1E3/UV,X0+3);
  as=vnorm_c(X0);
  tdyn=2*M_PI*sqrt(as*as*as/(GGLOBAL*MSUN/UM));
  h=tdyn/1000.0;

  tfut=to+duration*UT;
  durasymp=duration;
  dtdur=0.05*YEAR;
  copyVec(X0s,X0,6);

  tdur=direction*0.01*YEAR;
  do{
    VPRINT(stdout,"Dur.=%e\n",tdur/YEAR);

    //CALCULATE POSITION ASSUMING DURATION = TDUR
    integrateEoM(tini/UT,X0,h,npoints,tdur/UT,6,EoM,params,ts,Xout);
    vscl_c(UL/1E3,Xout[1],Xu);vscl_c(UV/1E3,Xout[1]+3,Xu+3);
    VPRINT(stdout,"Integration=%s\n",vec2strn(Xref,6,"%e "));
    oscelt_c(Xu,ts[1]*UT,munom,elts);

    //DETERMINE PREDICTED POSITION USING THE CONIC APPROXIMATION AND ITS ERROR
    conics_c(elts,tfut,Xend);
    VPRINT(stdout,"Conic approximation=%s\n",vec2strn(Xend,6,"%e "));
    vsub_c(Xend,Xref,Xdif);
    dp=vnorm_c(Xdif)*1e3/UL;
    VPRINT(stdout,"Difference in position at asymptotic time=%e\n",dp);
      
    //COMPARE DIFFERENCE BETWEEN CONIC AND RIGOROUS SOLUTION
    if(dp<(ds/10)){
      VPRINT(stdout,"Dur.Con. Found!\n");
      durasymp=tdur;
      dasymp=vnorm_c(Xu)*1e3/AU;
      dpasymp=dp;
      break;
    }
    tdur+=direction*dtdur;
    dtdur*=2;
  }while(direction*(tdur-1*duration*UT)<0);
  fprintf(stdout,"Conic approximation can be calculated from t = %e years, when dist. = %e (dp = %e)\n",durasymp/YEAR,dasymp,dpasymp);

  //ASYMPTOTIC ELEMENTS
  copyVec(elemasymp,elts,8);
  tend=tini+tdur;
  deltet_c(tend,"ET",&deltat);
  tend+=deltat;
  et2utc_c(tend,"C",2,100,dend);
  fprintf(stdout,"Asymptotic elements at t = %s: (q = %.10lf AU, e = %.10lf, i = %.10lf deg, Omega = %.10lf deg, omega = %.10lf deg, M = %.10lf deg, to = %.3lf, mu = %.10e)\n",
	  dend,elemasymp[0]*1e3/AU,elemasymp[1],elemasymp[2]*RAD,
	  elemasymp[3]*RAD,elemasymp[4]*RAD,elemasymp[5]*RAD,
	  unitim_c(elemasymp[6],"TDB","JDTDB"),elemasymp[7]);

  fprintf(fi,"dt_asy=%.17e\n",durasymp/YEAR);
  fprintf(fi,"q_asy=%.17e\n",elemasymp[0]*1e3/AU);
  fprintf(fi,"e_asy=%.17e\n",elemasymp[1]);
  fprintf(fi,"i_asy=%.17e\n",elemasymp[2]*RAD);
  fprintf(fi,"W_asy=%.17e\n",elemasymp[3]*RAD);
  fprintf(fi,"w_asy=%.17e\n",elemasymp[4]*RAD);
  fprintf(fi,"M_asy=%.17e\n",elemasymp[5]*RAD);
  fprintf(fi,"mu_asy=%.17e\n",elemasymp[7]);
  fprintf(fi,"t_asy=%.17e\n",unitim_c(elemasymp[6],"TDB","JDTDB"));
  fprintf(fi,"date_asy='%s'\n",dend);

  ////////////////////////////////////////////////////
  //COMPUTE ERROR OF ASYMPTOTIC ELEMENTS
  ////////////////////////////////////////////////////
  fprintf(stdout,"Calculating elements at t = %f\n",durasymp/YEAR);
  double qs[Ncov],es[Ncov],tps[Ncov],Ws[Ncov],ws[Ncov],is[Ncov];
  for(int j=0;j<=Ncov+1;j++){

    if(qdiagonal) gsl_matrix_memcpy(L,D);
    else gsl_matrix_memcpy(L,Lo);
    gsl_linalg_cholesky_decomp1(L);
    gsl_ran_multivariate_gaussian(RAND,ielements,L,relements);
    n=ini_n+gsl_ran_gaussian(RAND,ini_dn);
    tp=gsl_vector_get(relements,2);
    Mo=n*(ini_to_jed-tp);
    /*q=*/elements[0]=gsl_vector_get(relements,1)*AU/1e3;
    /*e=*/elements[1]=gsl_vector_get(relements,0);
    /*i=*/elements[2]=gsl_vector_get(relements,5)*DEG;
    /*W=*/elements[3]=gsl_vector_get(relements,3)*DEG;
    /*w=*/elements[4]=gsl_vector_get(relements,4)*DEG;
    /*M=*/elements[5]=Mo*DEG;
    /*to=*/elements[6]=to;
    /*mu=*/elements[7]=mu;

    /*
    //Uncomment if you want to check covariance computation
    qs[j]=elements[0]*1e3/AU;
    es[j]=elements[1];
    tps[j]=tp;
    is[j]=elements[2]*RAD;
    Ws[j]=elements[3]*RAD;
    ws[j]=elements[4]*RAD;
    */

    conics_c(elements,to,position);
    spkezr_c("SUN",to,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
    vaddg_c(position,sun,6,position);
    vpack_c(position[0]*1E3/UL,position[1]*1E3/UL,position[2]*1E3/UL,X0);
    vpack_c(position[3]*1E3/UV,position[4]*1E3/UV,position[5]*1E3/UV,X0+3);
    as=vnorm_c(X0);
    tdyn=2*M_PI*sqrt(as*as*as/(GGLOBAL*MSUN/UM));
    h=tdyn/1000.0;
    integrateEoM(tini/UT,X0,h,npoints,durasymp/UT,6,EoM,params,ts,Xout);
    vscl_c(UL/1E3,Xout[1],Xu);vscl_c(UV/1E3,Xout[1]+3,Xu+3);

    //CALCULATE ELEMENTS
    oscelt_c(Xu,to*DAY,munom,elts);
    qs[j]=elts[0]*1e3/AU;
    es[j]=elts[1];
    is[j]=elts[2]*RAD;
    Ws[j]=elts[3]*RAD;
    ws[j]=elts[4]*RAD;

    //Compute tp
    a=qs[j]*AU/1e3/(es[j]-1);//km
    n=sqrt(munom/(a*a*a))*RAD*DAY;//deg/day
    Mo=elts[5]*RAD;
    tp=ini_to_jed+durasymp/DAY-Mo/n;
    tps[j]=tp;
    
    if(j==0){
      copyVec(Xref,Xu,6);
    }else{
      vsub_c(Xref,Xu,Xdif);
      ds+=vnorm_c(Xdif);
    }
  }
  ds/=Ncov;
  ds*=1e3/AU;//AVERAGE DISTANCE IN AU
  fprintf(stdout,"Typical size of the cloud at asymptotic time: %e\n",2*ds);

  //Compute covariances
  fprintf(fi,"#actual covariance and original covariance\n");
  fprintf(fi,"cov_ee=%.17e,%.17e\n",gsl_stats_covariance(es,1,es,1,Ncov),ini_cov[0][0]);
  fprintf(fi,"cov_eq=%.17e,%.17e\n",gsl_stats_covariance(es,1,qs,1,Ncov),ini_cov[0][1]);
  fprintf(fi,"cov_et=%.17e,%.17e\n",gsl_stats_covariance(es,1,tps,1,Ncov),ini_cov[0][2]);
  fprintf(fi,"cov_eW=%.17e,%.17e\n",gsl_stats_covariance(es,1,Ws,1,Ncov),ini_cov[0][3]);
  fprintf(fi,"cov_ew=%.17e,%.17e\n",gsl_stats_covariance(es,1,ws,1,Ncov),ini_cov[0][4]);
  fprintf(fi,"cov_ei=%.17e,%.17e\n",gsl_stats_covariance(es,1,is,1,Ncov),ini_cov[0][5]);
  fprintf(fi,"cov_qq=%.17e,%.17e\n",gsl_stats_covariance(qs,1,qs,1,Ncov),ini_cov[1][1]);
  fprintf(fi,"cov_qt=%.17e,%.17e\n",gsl_stats_covariance(qs,1,tps,1,Ncov),ini_cov[1][2]);
  fprintf(fi,"cov_qW=%.17e,%.17e\n",gsl_stats_covariance(qs,1,Ws,1,Ncov),ini_cov[1][3]);
  fprintf(fi,"cov_qw=%.17e,%.17e\n",gsl_stats_covariance(qs,1,ws,1,Ncov),ini_cov[1][4]);
  fprintf(fi,"cov_qi=%.17e,%.17e\n",gsl_stats_covariance(qs,1,is,1,Ncov),ini_cov[1][5]);
  fprintf(fi,"cov_tt=%.17e,%.17e\n",gsl_stats_covariance(tps,1,tps,1,Ncov),ini_cov[2][2]);
  fprintf(fi,"cov_tW=%.17e,%.17e\n",gsl_stats_covariance(tps,1,Ws,1,Ncov),ini_cov[2][3]);
  fprintf(fi,"cov_tw=%.17e,%.17e\n",gsl_stats_covariance(tps,1,ws,1,Ncov),ini_cov[2][4]);
  fprintf(fi,"cov_ti=%.17e,%.17e\n",gsl_stats_covariance(tps,1,is,1,Ncov),ini_cov[2][5]);
  fprintf(fi,"cov_WW=%.17e,%.17e\n",gsl_stats_covariance(Ws,1,Ws,1,Ncov),ini_cov[3][3]);
  fprintf(fi,"cov_Ww=%.17e,%.17e\n",gsl_stats_covariance(Ws,1,ws,1,Ncov),ini_cov[3][4]);
  fprintf(fi,"cov_Wi=%.17e,%.17e\n",gsl_stats_covariance(Ws,1,is,1,Ncov),ini_cov[3][5]);
  fprintf(fi,"cov_ww=%.17e,%.17e\n",gsl_stats_covariance(ws,1,ws,1,Ncov),ini_cov[4][4]);
  fprintf(fi,"cov_wi=%.17e,%.17e\n",gsl_stats_covariance(ws,1,is,1,Ncov),ini_cov[4][5]);
  fprintf(fi,"cov_ii=%.17e,%.17e\n",gsl_stats_covariance(is,1,is,1,Ncov),ini_cov[5][5]);

  ////////////////////////////////////////////////////
  //COMPUTE INGRESS TIME OF NOMINAL OBJECT
  ////////////////////////////////////////////////////
  printHeader(stdout,"COMPUTING INGRESS TIME");
  conics_c(elemasymp,elemasymp[6],Xend);
  vasymp=vnorm_c(Xend+3);
  fprintf(stdout,"Asymptotic velocity: %e\n",vasymp);
  during=truncation/(vasymp*1e3);
  fprintf(stdout,"Estimated time of ingress: %e\n",during/YEAR);
  conics_c(elemasymp,elemasymp[6]+direction*during,Xend);
  ding=vnorm_c(Xend);
  fprintf(stdout,"Estimated distance at ingress: %e\n",ding*1e3/AU);
  while((ding*1e3/AU)<(truncation/AU)){
    durold=during;
    dold=ding;
    during+=1000*YEAR;
    conics_c(elemasymp,elemasymp[6]-during,Xend);
    ding=vnorm_c(Xend);
  }
  during=durold+(during-durold)/(ding-dold)*(truncation/1e3-dold);
  fprintf(stdout,"Time of ingress: %e\n",during/YEAR);
  fprintf(fi,"t_ing=%.17e\n",during/YEAR);

  //COMPUTING THE HELIOCENTRIC ACCELERATION
  conics_c(elemasymp,elemasymp[6]-during,Xend);
  double accel[6],accel_sun[6];
  spkezr_c("SUN",to,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
  vaddg_c(Xend,sun,6,position);
  vpack_c(position[0]*1E3/UL,position[1]*1E3/UL,position[2]*1E3/UL,X0);
  vpack_c(position[3]*1E3/UV,position[4]*1E3/UV,position[5]*1E3/UV,X0+3);
  EoM(0,X0,accel,params);
  vscl_c(UL/(UT*UT),accel+3,accel+3);
  fprintf(stdout,"Heliocentric acceleration at ingress (m/s^2): %s\n",
	  vec2strn(accel+3,3,"%.5e "));

  double asun=vnorm_c(accel+3);
  fprintf(stdout,"Magnitude of the heliocentric acceleration of the object (m/s^2): %e\n",
	  asun);

  //COMPUTING GALACTIC ACCELERATION
  double ul=UL,um=UM,ut=UT,uv=UV,gglobal=GGLOBAL;
  UL=PARSEC;UM=MSUN;UT=YEAR;UV=UL/UT;GCONST/(UL*UL*UL/(UM*UT*UT));

  //OBJECT
  mxv_c(M_Eclip_Galactic,Xend,posGalactic);
  mxv_c(M_Eclip_Galactic,Xend+3,posGalactic+3);
  LSR2GC(posGalactic,Xend);
  vscl_c(1e3/UL,Xend,Xend);//SET UNITS
  vscl_c(1e3/UV,Xend+3,Xend+3);
  cart2polar(Xend,X0,1.0);
  EoMGalactic(0,X0,accel,params);
  vscl_c(UL/(UT*UT),accel+3,accel+3);
  
  fprintf(stdout,"Galactocentric acceleration of the Object (m/s^2): %s\n",
	  vec2strn(accel+3,3,"%.5e "));

  //SUN
  double Xsun[]={0,0,0,0,0,0};
  mxv_c(M_Eclip_Galactic,Xsun,posGalactic);
  mxv_c(M_Eclip_Galactic,Xsun+3,posGalactic+3);
  LSR2GC(posGalactic,Xend);
  vscl_c(1e3/UL,Xend,Xend);//SET UNITS
  vscl_c(1e3/UV,Xend+3,Xend+3);
  cart2polar(Xend,X0,1.0);
  EoMGalactic(0,X0,accel_sun,params);
  vscl_c(UL/(UT*UT),accel_sun+3,accel_sun+3);

  fprintf(stdout,"Galactocentric acceleration of the Sun (m/s^2): %s\n",
	  vec2strn(accel_sun+3,3,"%.5e "));

  //Galactic acceleration of the object respect to the Sun
  vsub_c(accel+3,accel_sun+3,accel+3);
  fprintf(stdout,"Galactic acceleration of the object respect to the Sun (galactic coordinates) (m/s^2): %s\n",
	  vec2strn(accel+3,3,"%.5e "));

  double agal=vnorm_c(accel+3);
  fprintf(stdout,"Magnitude of the galactic acceleration of the object respect to the Sun (m/s^2): %e\n",
	  agal);
  
  fprintf(stdout,"Ratio Sun/Gal (>1 still in Solar System) = %e\n",asun/agal);

  //RECOVERING ORIGINAL UNITS
  UL=ul;UM=um;UT=ut;UV=uv;GGLOBAL=gglobal;

  ////////////////////////////////////////////////////
  //COMPUTE POSITION AT INGRESS
  ////////////////////////////////////////////////////
  printHeader(stdout,"GENERATING SURROGATE OBJECTS POSITION");

  sprintf(Filename,"scratch/wanderer-%s.csv",WANDERER);
  FILE* fc=fopen(Filename,"w");

  int Nfreq=ceil(Npart/10);
  Nfreq=Nfreq==0?1:Nfreq;

  //FILE HEADER
  fprintf(fc,"i,qo,eo,inco,Wo,wo,Mo,to,muo,qasy,easy,incasy,Wasy,wasy,Masy,tasy,mu,ting,xecl,yecl,zecl,vxecl,vyecl,vzecl,xsky,ysky,zsky,vxsky,vysky,vzsky,xgal,ygal,zgal,vxgal,vygal,vzgal,RA,DEC,l,b,dummy\n");
  for(int j=0;j<Npart;j++){
    
    if ((j%Nfreq)==0)
      fprintf(stdout,"Particle %d...\n",j+1);

    if(j>0){
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      //GENERATE INITIAL ELEMENTS (MULTIVAR)
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      //GENERATE ELEMENTS BY GAUSSAN MULTIVARIATE
      if(qdiagonal)
	gsl_matrix_memcpy(L,D);
      else
	gsl_matrix_memcpy(L,Lo);
      
      gsl_linalg_cholesky_decomp1(L);
      gsl_ran_multivariate_gaussian(RAND,ielements,L,relements);
      n=ini_n+gsl_ran_gaussian(RAND,ini_dn);
      tp=gsl_vector_get(relements,2);
      Mo=n*(ini_to_jed-tp);

      //EXTRACT ELEMENTS
      /*q=*/elements[0]=gsl_vector_get(relements,1)*AU/1e3;
      /*e=*/elements[1]=gsl_vector_get(relements,0);
      /*i=*/elements[2]=gsl_vector_get(relements,5)*DEG;
      /*W=*/elements[3]=gsl_vector_get(relements,3)*DEG;
      /*w=*/elements[4]=gsl_vector_get(relements,4)*DEG;
      /*M=*/elements[5]=Mo*DEG;
      /*to=*/elements[6]=to;
      /*mu=*/elements[7]=mu;
    }else{
      /*q=*/elements[0]=ini_q*AU/1e3;
      /*e=*/elements[1]=ini_e;
      /*i=*/elements[2]=ini_i*DEG;
      /*W=*/elements[3]=ini_W*DEG;
      /*w=*/elements[4]=ini_w*DEG;
      /*M=*/elements[5]=ini_M*DEG;
      /*to=*/elements[6]=to;
      /*mu=*/elements[7]=mu;
    }      

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INITIAL POSITION @ SSB
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    conics_c(elements,to,position);
    VPRINT(stdout,"\tInitial position @ Sun (EJ2000) : %s\n",
			vec2strn(position,3,"%.17e "));
    VPRINT(stdout,"\tInitial velocity @ Sun (EJ2000) : %s\n",
			vec2strn(position+3,3,"%.17e "));
    
    spkezr_c("SUN",to,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
    vaddg_c(position,sun,6,position);
    VPRINT(stdout,"\tInitial position @ SSB (EJ2000) : %s\n",
	    vec2strn(position,3,"%.17e "));
    VPRINT(stdout,"\tInitial velocity @ SSB (EJ2000) : %s\n",
	    vec2strn(position+3,3,"%.17e "));

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INITIAL CONDITIONS
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    vpack_c(position[0]*1E3/UL,position[1]*1E3/UL,position[2]*1E3/UL,X0);
    vpack_c(position[3]*1E3/UV,position[4]*1E3/UV,position[5]*1E3/UV,X0+3);
    as=vnorm_c(X0);
    tdyn=2*M_PI*sqrt(as*as*as/(GGLOBAL*MSUN/UM));

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INTEGRATION UNTIL TASYMP
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    h=tdyn/1000.0;
    integrateEoM(tini/UT,X0,h,npoints,durasymp/UT,6,EoM,params,ts,Xout);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //FINAL POSITION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    vscl_c(UL/1E3,Xout[1],Xu);vscl_c(UV/1E3,Xout[1]+3,Xu+3);
    VPRINT(stdout,"\tFinal position : %s\n",
	    vec2strn(Xu,3,"%.17e "));
    double d=vnorm_c(Xu)*1e3/AU;
    VPRINT(stdout,"\tDistance : %e\n",d);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //ASYMPTOTIC ELEMENTS
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    oscelt_c(Xu,ts[1],munom,elts);
    VPRINT(stdout,"\tAsymptotic elements : %s\n",vec2strn(elts,8,"%e "));

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INGRESS POSITION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    conics_c(elts,elts[6]+direction*during,Xu);
    ding=vnorm_c(Xu);
    VPRINT(stdout,"\tIngress distance : %e\n",ding*1e3/AU);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //POSITION AT INGRESS
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //J2000
    mxv_c(M_Eclip_J2000,Xu,posJ2000);
    mxv_c(M_Eclip_J2000,Xu+3,posJ2000+3);
    recrad_c(posJ2000,&d,&RA,&DEC);
    VPRINT(stdout,"\tRA(+HH:MM:SS) = %s, DEC(DD:MM:SS) = %s, d = %.3lf AU\n",
	    dec2sex(RA*RAD/15.0),dec2sex(DEC*RAD),d*1e3/AU);
    
    //GALACTIC
    mxv_c(M_Eclip_Galactic,Xu,posGalactic);
    mxv_c(M_Eclip_Galactic,Xu+3,posGalactic+3);
    recrad_c(posGalactic,&d,&l,&b);
    VPRINT(stdout,"\tl(+DD:MM:SS) = %s, b(DD:MM:SS) = %s\n",
	    dec2sex(l*RAD),dec2sex(b*RAD));
    
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //SAVE INFORMATION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //ELEMENTOS
    fprintf(fc,"%d,",j);
    //INITIAL ELEMENTS
    fprintf(fc,"%s",vec2strn(elements,8,"%.17e,"));
    //ASYMPTOTIC ELEMENTS
    fprintf(fc,"%s",vec2strn(elements,8,"%.17e,"));
    //TIME OF INGRESS
    fprintf(fc,"%.17e,",elts[6]+direction*during);
    //POSITION ECLIPJ2000
    fprintf(fc,"%s",vec2strn(Xu,6,"%.17e,"));
    //POSITION J2000
    fprintf(fc,"%s",vec2strn(posJ2000,6,"%.17e,"));
    //POSITION GALACTIC
    fprintf(fc,"%s",vec2strn(posGalactic,6,"%.17e,"));
    //J2000
    fprintf(fc,"%.17e,%.17e,",RA*RAD,DEC*RAD);
    //GALACTIC
    fprintf(fc,"%.17e,%.17e,",l*RAD,b*RAD);
    fprintf(fc,"\n");
  }
  fclose(fc);
  fclose(fi);

  printHeader(stdout,"DONE.",'!');
}
