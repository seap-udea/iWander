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

#define VERBOSE 2 //Verbosity level
#define OSTREAM stdout //Stream where the output is redirected
#define VSTREAM stderr //Stream where the error output is redirected

int main(int argc,char* argv[])
{
  /*
    Example: ./wanderer.exe
    
    Function: 

    This program perform three different tasks:

    1) Calculate the time t_asy when the single conic approximation is
       good enough to predict the future position of the interstellar
       object.

    2) Calculate the time t_ing at which the object was at a distance
       equivalent to the truncation tidal radius of the Solar System.

    3) Predict the position and velocity of the surrogate objects at
       t_ing.

    Input: None

    Output: 

    * wanderer.csv: properties of all the surrogate objects.

    * ingress.dat: a summary of the ingress orbit properties including
      the epoch of asymptotic elements and their covariance matrix,
      the time of ingress, the radiant and velocity at ingress.

  */

  ////////////////////////////////////////////////////
  //CONFIGURATION
  ////////////////////////////////////////////////////
  #include <iwander.conf>
  #include <wanderer.conf>
  printHeader(OSTREAM,"WANDERER SOLAR SYSTEM DYNAMICS",'*');

  ////////////////////////////////////////////////////
  //INITIALIZE iWANDER
  ////////////////////////////////////////////////////
  initWander();
  
  ////////////////////////////////////////////////////
  //VARIABLES DECLARATION
  ////////////////////////////////////////////////////
  //COUNTERS
  int i,j;

  //INTEGRATOR VARIABLES
  double t,to,tp,dt,h,ts[2];
  double tini,t_start,t_step,tend,t_stop;
  double h_used,h_next,h_adjust,deltat;
  double during,durold,ting,ding,dold,vasymp,direction;
  double tdur,dtdur,tfut,durasymp,dasymp,tasymp,dtasymp;
  double params[]={6};
  int npoints=2;

  //VECTOR AND MATRICES

  //ORBITAL COMPUTATION
  double q,a,e,n,mu,munom,inc,W,w,Mo;
  SpiceDouble position[6],sun[6];
  SpiceDouble posJ2000[6],posGalactic[6];
  SpiceDouble M_Eclip_J2000[3][3],M_Eclip_Galactic[3][3];
  double elts[8],elemasymp[8];
  SpiceDouble elements[8];
  double as,tdyn;
  SpiceChar dend[100];
  double ds,dp,dpasymp;
  double qs[Ncov],es[Ncov],tps[Ncov],Ws[Ncov],ws[Ncov],is[Ncov];
  double truncation=RTRUNC*AU;
  pxform_c("ECLIPJ2000","GALACTIC",t,M_Eclip_Galactic);
  pxform_c("ECLIPJ2000","J2000",0,M_Eclip_J2000);
  
  //POSITION PARAMETERS
  double X[6],X0[6],Xu[6],Xref[6],Xdif[6],X0s[6],Xend[6];
  
  //COORDINATES
  double RA,DEC,d,l,b;

  //TEMPORAL
  double ltmp;

  //BEHAVIOR
  int status;

  //ALLOCATION
  gsl_vector* relements=gsl_vector_alloc(6);
  gsl_matrix* Lo=gsl_matrix_alloc(6,6);
  gsl_matrix* L=gsl_matrix_alloc(6,6);
  gsl_vector* ielements=gsl_vector_alloc(6);
  gsl_matrix* D=gsl_matrix_alloc(6,6);
  double **Xout=matrixAllocate(2,6);

  //ACCELERATION
  double accel[6],accel_sun[6],agal;
  double ul,um,ut,uv,gglobal;
  double Xsun[]={0,0,0,0,0,0};

  ////////////////////////////////////////////////////
  //UNITS
  ////////////////////////////////////////////////////
  UL=Ul;
  UM=Um;
  UT=sqrt(UL*UL*UL/(GCONST*UM));
  UV=UL/UT;
  GGLOBAL=Gglobal;
  
  ////////////////////////////////////////////////////
  //INPUT ELEMENTS
  ////////////////////////////////////////////////////
  printHeader(stdout,"READING INPUT ELEMENTS",'-');
  print1(VSTREAM,"READING INPUT ELEMENTS",'-');

  //GRAVITATIONAL PARAMETER
  n=ini_n*DEG/DAY;
  a=ini_a*(-AU/1E3);
  mu=n*n*a*a*a;
  munom=mu;

  //EPHEMERIS TIME
  to=unitim_c(ini_to_jed,"JDTDB","TDB");
  print1(VSTREAM,"\tTime of elements=%.17e s\n",to);

  //INTEGRATION PARAMETERS
  tini=to;
  direction=Duration/fabs(Duration);
  Duration/=UT;
  gsl_vector_set(ielements,0,ini_e);
  gsl_vector_set(ielements,1,ini_q);
  gsl_vector_set(ielements,2,ini_tp);
  gsl_vector_set(ielements,3,ini_W);
  gsl_vector_set(ielements,4,ini_w);
  gsl_vector_set(ielements,5,ini_i);
  for(i=0;i<6;i++) for(j=0;j<6;j++) gsl_matrix_set(Lo,i,j,ini_cov[i][j]);

  //DIAGONAL COVARIANCE
  gsl_matrix_set_zero(D);
  gsl_matrix_set(D,0,0,ini_de*ini_de);
  gsl_matrix_set(D,1,1,ini_dq*ini_dq);
  gsl_matrix_set(D,2,2,ini_cov[2][2]);
  gsl_matrix_set(D,3,3,ini_dW*ini_dW);
  gsl_matrix_set(D,4,4,ini_dw*ini_dw);
  gsl_matrix_set(D,5,5,ini_di*ini_di);

  sprintf(FILENAME,"scratch/ingress-%s.dat",Wanderer);
  FILE* fi=fopen(FILENAME,"w");
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

  print1(VSTREAM,"\tEpoch=%.1f\n",ini_to_jed);
  print1(VSTREAM,"\tEpoch_Date='%s'\n",ini_to_date);
  print1(VSTREAM,"\tReference='%s'\n",reference_solution);
  print1(VSTREAM,"\ttruncation=%.1e\n",truncation/AU);
  print1(VSTREAM,"\tq_nom=%.17e+/-%.17e\n",ini_q,ini_dq);
  print1(VSTREAM,"\te_nom=%.17e+/-%.17e\n",ini_e,ini_de);
  print1(VSTREAM,"\ti_nom=%.17e+/-%.17e\n",ini_i,ini_di);
  print1(VSTREAM,"\tW_nom=%.17e+/-%.17e\n",ini_W,ini_dW);
  print1(VSTREAM,"\tw_nom=%.17e+/-%.17e\n",ini_w,ini_dw);
  print1(VSTREAM,"\tM_nom=%.17e+/-%.17e\n",ini_M,ini_dM);
  
  ////////////////////////////////////////////////////
  //COMPUTING CLOUD SIZE IN THE PAST
  ////////////////////////////////////////////////////
  printHeader(stdout,"COMPUTING CLOUD SIZE",'-');
  print1(VSTREAM,"COMPUTING CLOUD SIZE\n");

  print0(OSTREAM,"\tBack time, t = %.1f years\n",Duration*UT/YEAR);

  for(j=0;j<=Ndisp+1;j++){

    double q,a,n;
    if(j>0){

      //Qdiagonal?
      if(Qdiagonal) gsl_matrix_memcpy(L,D);
      else gsl_matrix_memcpy(L,Lo);
      
      //Generate random elements
      gsl_linalg_cholesky_decomp1(L);
      gsl_ran_multivariate_gaussian(RAND,ielements,L,relements);
      double e=gsl_vector_get(relements,0);
      if(e<1){
	j--;
	continue;
      }
      
      //Generate mean orbital motion
      //double ng=ini_n+gsl_ran_gaussian(RAND,ini_dn);
      q=gsl_vector_get(relements,1)*AU;//m
      a=q/(1-e);//m
      n=sqrt((MUTOT*1e9)/(fabs(a)*fabs(a)*fabs(a)))*RAD*DAY;//deg/day
      tp=gsl_vector_get(relements,2);
      Mo=n*(ini_to_jed-tp);
      //printf("\ttp = %e, to = %e\n",tp,ini_to_jed);
      //printf("\tM = %e\n",Mo);
      
      /*q=*/elements[0]=gsl_vector_get(relements,1)*AU/1e3;
      /*e=*/elements[1]=gsl_vector_get(relements,0);
      /*i=*/elements[2]=gsl_vector_get(relements,5)*DEG;
      /*W=*/elements[3]=gsl_vector_get(relements,3)*DEG;
      /*w=*/elements[4]=gsl_vector_get(relements,4)*DEG;
      /*M=*/elements[5]=Mo*DEG;
      /*to=*/elements[6]=to;
      /*mu=*/elements[7]=mu;
      print2(OSTREAM,"\tSurrogate %d: %s\n",j,vec2strn(elements,8,"%.2e "));

    }else{
      //The nominal object use the nominal elements
      /*q=*/elements[0]=ini_q*AU/1e3;
      /*e=*/elements[1]=ini_e;
      /*i=*/elements[2]=ini_i*DEG;
      /*W=*/elements[3]=ini_W*DEG;
      /*w=*/elements[4]=ini_w*DEG;
      /*M=*/elements[5]=ini_M*DEG;
      /*to=*/elements[6]=to;
      /*mu=*/elements[7]=mu;
      print2(OSTREAM,"\tNominal: %s\n",vec2strn(elements,8,"%.2e "));
    }

    //Convert from conic elements to rectangular position
    conics_c(elements,to,position);

    //Calculate the position of the Sun
    spkezr_c("SUN",to,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
    vaddg_c(position,sun,6,position);
    vpack_c(position[0]*1E3/UL,position[1]*1E3/UL,position[2]*1E3/UL,X0);
    vpack_c(position[3]*1E3/UV,position[4]*1E3/UV,position[5]*1E3/UV,X0+3);

    //Calculate the semimajor axis
    as=vnorm_c(X0);

    //Estimate the dynamic time for this orbit
    tdyn=2*M_PI*sqrt(as*as*as/(GGLOBAL*MSUN/UM));

    //Initial time-step
    h=MIN(tdyn/1000.0,-Duration/100);

    //Integrate surrogate objects
    integrateEoM(tini/UT,X0,h,npoints,Duration,6,EoM,params,ts,Xout);
    
    //Convert to km, km/s
    vscl_c(UL/1E3,Xout[1],Xu);vscl_c(UV/1E3,Xout[1]+3,Xu+3);

    //Save the reference position
    if(j==0){
      copyVec(Xref,Xu,6);
    }else{
      //Distance to reference objects
      vsub_c(Xref,Xu,Xdif);
      ds+=vnorm_c(Xdif);
    }
  }
  ds/=Ndisp;
  ds*=1e3/AU;//AVERAGE DISTANCE IN AU
  print0(OSTREAM,"\tSize of the cloud at t = %e, d = %e AU: %e AU\n",
	 Duration*UT/YEAR,vnorm_c(Xu)*1e3/AU,2*ds);
  
  //Save report in ingress file
  fprintf(fi,"t_test=%.17e\n",Duration*UT/YEAR);
  fprintf(fi,"d_test=%.17e\n",vnorm_c(Xu)*1e3/AU);
  fprintf(fi,"dr_test=%.17e\n",2*ds);

  //Calculate date of end
  tend=tini+Duration*UT;
  deltet_c(tend,"ET",&deltat);
  tend+=deltat;
  et2utc_c(tend,"C",2,100,dend);
  print1(VSTREAM,"\tPosition of nominal object at t = %e yr (%.17e, %s):\n",
	 Duration*UT/YEAR,tend,dend);

  vscl_c(1e3/AU,Xref,Xend);
  print1(VSTREAM,"\t\tr (AU) = %s\n",vec2strn(Xend,3,"%.17e "));
  print1(VSTREAM,"\t\tv (km/s) = %s\n",vec2strn(Xend+3,3,"%.17e "));

  ////////////////////////////////////////////////////
  //CALCULATE TIME OF ASYMPTOTIC ELEMENTS
  ////////////////////////////////////////////////////
  printHeader(stdout,"COMPUTING ASYMPTOTIC ELEMENTS",'-');
  print1(VSTREAM,"COMPUTING ASYMPTOTIC ELEMENTS\n");

  //Nominal elements
  /*q=*/elements[0]=ini_q*AU/1e3;
  /*e=*/elements[1]=ini_e;
  /*i=*/elements[2]=ini_i*DEG;
  /*W=*/elements[3]=ini_W*DEG;
  /*w=*/elements[4]=ini_w*DEG;
  /*M=*/elements[5]=ini_M*DEG;
  /*to=*/elements[6]=to;
  /*mu=*/elements[7]=mu;

  //Cartesian position
  conics_c(elements,to,position);

  //Solar position
  spkezr_c("SUN",to,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
  vaddg_c(position,sun,6,position);
  vpack_c(position[0]*1E3/UL,position[1]*1E3/UL,position[2]*1E3/UL,X0);
  vpack_c(position[3]*1E3/UV,position[4]*1E3/UV,position[5]*1E3/UV,X0+3);
  as=vnorm_c(X0);
  tdyn=2*M_PI*sqrt(as*as*as/(GGLOBAL*MSUN/UM));
  h=MIN(tdyn/1000.0,-Duration/100);

  //Estimate time of asymptotic elements
  tfut=to+Duration*UT;
  durasymp=Duration;
  dtdur=0.05*YEAR;
  copyVec(X0s,X0,6);

  tdur=direction*0.01*YEAR;

  print0(VSTREAM,"\tSearching for asymptotic elements\n");
  do{
    print1(VSTREAM,"\t\tDur.=%e\n",tdur/YEAR);

    //CALCULATE POSITION ASSUMING DURATION = TDUR
    integrateEoM(tini/UT,X0,h,npoints,tdur/UT,6,EoM,params,ts,Xout);
    vscl_c(UL/1E3,Xout[1],Xu);vscl_c(UV/1E3,Xout[1]+3,Xu+3);
    print1(VSTREAM,"\t\t\tIntegration=%s\n",vec2strn(Xu,6,"%e "));
    oscelt_c(Xu,ts[1]*UT,munom,elts);

    //DETERMINE PREDICTED POSITION USING THE CONIC APPROXIMATION AND ITS ERROR
    conics_c(elts,tfut,Xend);
    print1(VSTREAM,"\t\t\tConic approximation=%s\n",vec2strn(Xend,6,"%e "));
    vsub_c(Xend,Xref,Xdif);
    dp=vnorm_c(Xdif)*1e3/UL;
    print1(VSTREAM,"\t\t\tDifference in position at asymptotic time=%e\n",dp);
      
    //COMPARE DIFFERENCE BETWEEN CONIC AND RIGOROUS SOLUTION
    if(dp<(ds/10)){
      print1(VSTREAM,"\t\tAsymptotic elements found!\n");
      durasymp=tdur;
      dasymp=vnorm_c(Xu)*1e3/AU;
      dpasymp=dp;
      break;
    }
    tdur+=direction*dtdur;
    dtdur*=2;
  }while(direction*(tdur-1*Duration*UT)<0);
  print0(OSTREAM,"\tConic approximation can be calculated from t = %e years, when dist. = %e AU (Cloud size = %e AU)\n",durasymp/YEAR,dasymp,dpasymp);
  
  //ASYMPTOTIC ELEMENTS
  copyVec(elemasymp,elts,8);
  tend=tini+tdur;
  deltet_c(tend,"ET",&deltat);
  tend+=deltat;
  et2utc_c(tend,"C",2,100,dend);

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

  print1(VSTREAM,"\tAsymptotic elements at t = %s:\n",dend);
  print1(VSTREAM,"\t\tdt_asy=%.17e\n",durasymp/YEAR);
  print1(VSTREAM,"\t\tq_asy=%.17e\n",elemasymp[0]*1e3/AU);
  print1(VSTREAM,"\t\te_asy=%.17e\n",elemasymp[1]);
  print1(VSTREAM,"\t\ti_asy=%.17e\n",elemasymp[2]*RAD);
  print1(VSTREAM,"\t\tW_asy=%.17e\n",elemasymp[3]*RAD);
  print1(VSTREAM,"\t\tw_asy=%.17e\n",elemasymp[4]*RAD);
  print1(VSTREAM,"\t\tM_asy=%.17e\n",elemasymp[5]*RAD);
  print1(VSTREAM,"\t\tmu_asy=%.17e\n",elemasymp[7]);
  print1(VSTREAM,"\t\tt_asy=%.17e\n",unitim_c(elemasymp[6],"TDB","JDTDB"));
  print1(VSTREAM,"\t\tdate_asy='%s'\n",dend);

  if(fabs(durasymp/YEAR)<0.1){
    //Compute covariances
    fprintf(fi,"#actual covariance and original covariance\n");
    fprintf(fi,"cov_ee=%.17e,%.17e\n",ini_cov[0][0],ini_cov[0][0]);
    fprintf(fi,"cov_eq=%.17e,%.17e\n",ini_cov[0][1],ini_cov[0][1]);
    fprintf(fi,"cov_et=%.17e,%.17e\n",ini_cov[0][2],ini_cov[0][2]);
    fprintf(fi,"cov_eW=%.17e,%.17e\n",ini_cov[0][3],ini_cov[0][3]);
    fprintf(fi,"cov_ew=%.17e,%.17e\n",ini_cov[0][4],ini_cov[0][4]);
    fprintf(fi,"cov_ei=%.17e,%.17e\n",ini_cov[0][5],ini_cov[0][5]);
    fprintf(fi,"cov_qq=%.17e,%.17e\n",ini_cov[1][1],ini_cov[1][1]);
    fprintf(fi,"cov_qt=%.17e,%.17e\n",ini_cov[1][2],ini_cov[1][2]);
    fprintf(fi,"cov_qW=%.17e,%.17e\n",ini_cov[1][3],ini_cov[1][3]);
    fprintf(fi,"cov_qw=%.17e,%.17e\n",ini_cov[1][4],ini_cov[1][4]);
    fprintf(fi,"cov_qi=%.17e,%.17e\n",ini_cov[1][5],ini_cov[1][5]);
    fprintf(fi,"cov_tt=%.17e,%.17e\n",ini_cov[2][2],ini_cov[2][2]);
    fprintf(fi,"cov_tW=%.17e,%.17e\n",ini_cov[2][3],ini_cov[2][3]);
    fprintf(fi,"cov_tw=%.17e,%.17e\n",ini_cov[2][4],ini_cov[2][4]);
    fprintf(fi,"cov_ti=%.17e,%.17e\n",ini_cov[2][5],ini_cov[2][5]);
    fprintf(fi,"cov_WW=%.17e,%.17e\n",ini_cov[3][3],ini_cov[3][3]);
    fprintf(fi,"cov_Ww=%.17e,%.17e\n",ini_cov[3][4],ini_cov[3][4]);
    fprintf(fi,"cov_Wi=%.17e,%.17e\n",ini_cov[3][5],ini_cov[3][5]);
    fprintf(fi,"cov_ww=%.17e,%.17e\n",ini_cov[4][4],ini_cov[4][4]);
    fprintf(fi,"cov_wi=%.17e,%.17e\n",ini_cov[4][5],ini_cov[4][5]);
    fprintf(fi,"cov_ii=%.17e,%.17e\n",ini_cov[5][5],ini_cov[5][5]);
    fflush(fi);
  }else{
    ////////////////////////////////////////////////////
    //COMPUTE ERROR OF ASYMPTOTIC ELEMENTS
    ////////////////////////////////////////////////////
    printHeader(stdout,"COMPUTING ERRORS IN ASYMPTOTIC ELEMENTS",'-');
    print1(VSTREAM,"COMPUTING ERRORS IN ASYMPTOTIC ELEMENTS\n");
  
    print0(OSTREAM,"\tUsing %d objects\n",Ncov);
    for(j=0;j<=Ncov+1;j++){

      //Generate elements
      if(Qdiagonal) gsl_matrix_memcpy(L,D);
      else gsl_matrix_memcpy(L,Lo);
      gsl_linalg_cholesky_decomp1(L);
      gsl_ran_multivariate_gaussian(RAND,ielements,L,relements);
      n=ini_n+gsl_ran_gaussian(RAND,ini_dn);
      tp=gsl_vector_get(relements,2);
      Mo=n*(ini_to_jed-tp);

      //Save elements
      /*q=*/elements[0]=gsl_vector_get(relements,1)*AU/1e3;
      /*e=*/elements[1]=gsl_vector_get(relements,0);
      /*i=*/elements[2]=gsl_vector_get(relements,5)*DEG;
      /*W=*/elements[3]=gsl_vector_get(relements,3)*DEG;
      /*w=*/elements[4]=gsl_vector_get(relements,4)*DEG;
      /*M=*/elements[5]=Mo*DEG;
      /*to=*/elements[6]=to;
      /*mu=*/elements[7]=mu;

      //Compute position
      conics_c(elements,to,position);
      spkezr_c("SUN",to,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
      vaddg_c(position,sun,6,position);
      vpack_c(position[0]*1E3/UL,position[1]*1E3/UL,position[2]*1E3/UL,X0);
      vpack_c(position[3]*1E3/UV,position[4]*1E3/UV,position[5]*1E3/UV,X0+3);

      //Integrate 
      as=vnorm_c(X0);
      tdyn=2*M_PI*sqrt(as*as*as/(GGLOBAL*MSUN/UM));
      h=MIN(tdyn/1000.0,-durasymp/1000.0);
      integrateEoM(tini/UT,X0,h,npoints,durasymp/UT,6,EoM,params,ts,Xout);
      vscl_c(UL/1E3,Xout[1],Xu);vscl_c(UV/1E3,Xout[1]+3,Xu+3);
    
      //Calculate asymptotic elements
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
    print0(OSTREAM,"\tTypical size of the cloud at asymptotic time: %e\n",2*ds);

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
    fflush(fi);
  }

  ////////////////////////////////////////////////////
  //COMPUTE INGRESS TIME OF NOMINAL OBJECT
  ////////////////////////////////////////////////////
  printHeader(stdout,"COMPUTING INGRESS TIME",'-');
  print1(VSTREAM,"COMPUTING INGRESS TIME\n");

  //Compute asymptotic position
  conics_c(elemasymp,elemasymp[6],Xend);
  vasymp=vnorm_c(Xend+3);
  print0(OSTREAM,"\tAsymptotic velocity: %e\n",vasymp);
  during=truncation/(vasymp*1e3);
  print0(OSTREAM,"\tEstimated time of ingress: %e\n",during/YEAR);
  conics_c(elemasymp,elemasymp[6]+direction*during,Xend);
  ding=vnorm_c(Xend);
  print0(OSTREAM,"\tEstimated distance at ingress: %e\n",ding*1e3/AU);
  while((ding*1e3/AU)<(truncation/AU)){
    durold=during;
    dold=ding;
    during+=1000*YEAR;
    conics_c(elemasymp,elemasymp[6]-during,Xend);
    ding=vnorm_c(Xend);
  }
  during=durold+(during-durold)/(ding-dold)*(truncation/1e3-dold);
  print0(OSTREAM,"\tTime of ingress: %e\n",during/YEAR);
  fprintf(fi,"t_ing=%.17e\n",during/YEAR);
  
  ////////////////////////////////////////////////////
  //COMPUTE POSITION AT INGRESS
  ////////////////////////////////////////////////////
  printHeader(stdout,"COMPUTING ACCELERATION AT INGRESS TIME",'-');
  print1(VSTREAM,"COMPUTING ACCELERATION AT INGRESS TIME\n");

  conics_c(elemasymp,elemasymp[6]-during,Xend);
  spkezr_c("SUN",to,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
  vaddg_c(Xend,sun,6,position);
  vpack_c(position[0]*1E3/UL,position[1]*1E3/UL,position[2]*1E3/UL,X0);
  vpack_c(position[3]*1E3/UV,position[4]*1E3/UV,position[5]*1E3/UV,X0+3);
  EoM(0,X0,accel,params);
  vscl_c(UL/(UT*UT),accel+3,accel+3);
  print0(OSTREAM,"\tHeliocentric acceleration at ingress (m/s^2): %s\n",
	  vec2strn(accel+3,3,"%.5e "));

  double asun=vnorm_c(accel+3);
  print0(OSTREAM,"\tMagnitude of the heliocentric acceleration of the object (m/s^2): %e\n",
	  asun);

  //COMPUTING GALACTIC ACCELERATION
  ul=UL;um=UM;ut=UT;uv=UV;gglobal=GGLOBAL;
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
  
  print0(OSTREAM,"\tGalactocentric acceleration of the Object (m/s^2): %s\n",
	  vec2strn(accel+3,3,"%.5e "));

  //SUN
  mxv_c(M_Eclip_Galactic,Xsun,posGalactic);
  mxv_c(M_Eclip_Galactic,Xsun+3,posGalactic+3);
  LSR2GC(posGalactic,Xend);
  vscl_c(1e3/UL,Xend,Xend);//SET UNITS
  vscl_c(1e3/UV,Xend+3,Xend+3);
  cart2polar(Xend,X0,1.0);
  EoMGalactic(0,X0,accel_sun,params);
  vscl_c(UL/(UT*UT),accel_sun+3,accel_sun+3);

  print0(OSTREAM,"\tGalactocentric acceleration of the Sun (m/s^2): %s\n",
	  vec2strn(accel_sun+3,3,"%.5e "));

  //Galactic acceleration of the object respect to the Sun
  vsub_c(accel+3,accel_sun+3,accel+3);
  print0(OSTREAM,"\tGalactic acceleration of the object respect to the Sun (galactic coordinates) (m/s^2): %s\n",
	  vec2strn(accel+3,3,"%.5e "));

  agal=vnorm_c(accel+3);
  print0(OSTREAM,"\tMagnitude of the galactic acceleration of the object respect to the Sun (m/s^2): %e\n",
	  agal);
  
  print0(OSTREAM,"\tRatio Sun/Gal (>1 still in Solar System) = %e\n",asun/agal);

  //RECOVERING ORIGINAL UNITS
  UL=ul;UM=um;UT=ut;UV=uv;GGLOBAL=gglobal;

  ////////////////////////////////////////////////////
  //COMPUTE POSITION AT INGRESS
  ////////////////////////////////////////////////////
  printHeader(stdout,"GENERATING SURROGATE OBJECTS POSITION");
  print1(VSTREAM,"GENERATING SURROGATE OBJECTS POSITION\n");
  
  sprintf(FILENAME,"scratch/wanderer-%s.csv",Wanderer);
  FILE* fc=fopen(FILENAME,"w");

  int Nfreq=ceil(Npart/10);
  Nfreq=Nfreq==0?1:Nfreq;

  //FILE HEADER
  fprintf(fc,"i,qo,eo,inco,Wo,wo,Mo,to,muo,qasy,easy,incasy,Wasy,wasy,Masy,tasy,mu,ting,xecl,yecl,zecl,vxecl,vyecl,vzecl,xsky,ysky,zsky,vxsky,vysky,vzsky,xgal,ygal,zgal,vxgal,vygal,vzgal,RA,DEC,l,b,dummy\n");
  for(j=0;j<Npart;j++){
    
    print2(VSTREAM,"\tParticle %d...\n",j+1);
    if ((j%Nfreq)==0)
      print0(OSTREAM,"\tParticle %d...\n",j+1);

    if(j>0){
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      //GENERATE INITIAL ELEMENTS (MULTIVAR)
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      //GENERATE ELEMENTS BY GAUSSAN MULTIVARIATE
      if(Qdiagonal)
	gsl_matrix_memcpy(L,D);
      else
	gsl_matrix_memcpy(L,Lo);
      
      gsl_linalg_cholesky_decomp1(L);
      gsl_ran_multivariate_gaussian(RAND,ielements,L,relements);
      n=ini_n+gsl_ran_gaussian(RAND,ini_dn);
      tp=gsl_vector_get(relements,2);
      Mo=n*(ini_to_jed-tp);

      double e=gsl_vector_get(relements,0);
      if(e<1){j--;continue;}
      
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
    print2(VSTREAM,"\t\tInitial position @ Sun (EJ2000) : %s\n",
			vec2strn(position,3,"%.17e "));
    print2(VSTREAM,"\t\tInitial velocity @ Sun (EJ2000) : %s\n",
			vec2strn(position+3,3,"%.17e "));
    
    spkezr_c("SUN",to,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
    vaddg_c(position,sun,6,position);
    print2(VSTREAM,"\t\tInitial position @ SSB (EJ2000) : %s\n",
	    vec2strn(position,3,"%.17e "));
    print2(VSTREAM,"\t\tInitial velocity @ SSB (EJ2000) : %s\n",
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
    try{
      integrateEoM(tini/UT,X0,h,npoints,durasymp/UT,6,EoM,params,ts,Xout);
    }catch(int e){
      print2(VSTREAM,"\t\tIntegration failed.\n");
      continue;
    }

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //FINAL POSITION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    vscl_c(UL/1E3,Xout[1],Xu);vscl_c(UV/1E3,Xout[1]+3,Xu+3);
    print2(VSTREAM,"\t\tFinal position : %s\n",
	    vec2strn(Xu,3,"%.17e "));
    double d=vnorm_c(Xu)*1e3/AU;
    print2(VSTREAM,"\t\tDistance : %e\n",d);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //ASYMPTOTIC ELEMENTS
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    oscelt_c(Xu,ts[1],munom,elts);
    print2(VSTREAM,"\t\tAsymptotic elements : %s\n",vec2strn(elts,8,"%e "));

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INGRESS POSITION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    conics_c(elts,elts[6]+direction*during,Xu);
    ding=vnorm_c(Xu);
    print2(VSTREAM,"\t\tIngress distance : %e\n",ding*1e3/AU);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //POSITION AT INGRESS
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //J2000
    mxv_c(M_Eclip_J2000,Xu,posJ2000);
    mxv_c(M_Eclip_J2000,Xu+3,posJ2000+3);
    recrad_c(posJ2000,&d,&RA,&DEC);
    print2(VSTREAM,"\t\tRA(+HH:MM:SS) = %s, DEC(DD:MM:SS) = %s, d = %.3lf AU\n",
	    dec2sex(RA*RAD/15.0),dec2sex(DEC*RAD),d*1e3/AU);
    
    //GALACTIC
    mxv_c(M_Eclip_Galactic,Xu,posGalactic);
    mxv_c(M_Eclip_Galactic,Xu+3,posGalactic+3);
    recrad_c(posGalactic,&d,&l,&b);
    print2(VSTREAM,"\t\tl(+DD:MM:SS) = %s, b(DD:MM:SS) = %s\n",
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

  TELAPS=elapsedTime(0);
  print0(OSTREAM,"Total elapsed time = %.5f (%.5f min)\n",TELAPS,TELAPS/60.0);
  print0(OSTREAM,"DONE.\n");
  return 0;
}
