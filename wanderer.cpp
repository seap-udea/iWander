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
  double h_used,h_next,h_adjust,deltat;
  int i,status;
  //INITIAL CONDITIONS
  double X0[6],X[6],Xu[6],as,tdyn;
  //ELEMENTS
  double elts[8];
  double ds,dp,Xref[6],Xdif[6],Ndisp=50,dpasymp,elemasymp[8];
  double tdur,dtdur,X0s[6],Xend[6],tfut,durasymp,dasymp,tasymp,dtasymp;
  double during,durold,ting,ding,dold,vasymp;

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

  ////////////////////////////////////////////////////
  //DISPERSION OF SURROGATE OBJECTS IN THE FUTURE
  ////////////////////////////////////////////////////
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
  fprintf(stdout,"Average distance between surrograte at (t = %e, d = %e): %e AU\n",duration*UT/YEAR,vnorm_c(Xu)*1e3/AU,ds);

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
  for(tdur=-0.01*YEAR;tdur>=duration*UT;tdur-=dtdur){

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
    VPRINT(stdout,"Dur.aymp.=%e\n",dp);
      
    //COMPARE DIFFERENCE BETWEEN CONIC AND RIGOROUS SOLUTION
    if(dp<(ds/10)){
      VPRINT(stdout,"Dur.Con. Found!\n");
      durasymp=tdur;
      dasymp=vnorm_c(Xu)*1e3/AU;
      dpasymp=dp;
      break;
    }
  }
  fprintf(stdout,"Conic approximation can be calculated from t = %e years, when dist. = %e (dp = %e)\n",durasymp/YEAR,dasymp,dpasymp);

  //ASYMPTOTIC ELEMENTS
  copyVec(elemasymp,elts,8);
  fprintf(stdout,"Asymptotic elements: (q = %.10lf AU, e = %.10lf, i = %.10lf deg, Omega = %.10lf deg, omega = %.10lf deg, M = %.10lf deg, to = %.3lf, mu = %.10e\n",
	  elemasymp[0]*1e3/AU,elemasymp[1],elemasymp[2]*RAD,
	  elemasymp[3]*RAD,elemasymp[4]*RAD,elemasymp[5]*RAD,
	  unitim_c(elemasymp[6],"TDB","JDTDB"),elemasymp[7]);

  ////////////////////////////////////////////////////
  //COMPUTE INGRESS TIME OF NOMINAL OBJECT
  ////////////////////////////////////////////////////
  conics_c(elemasymp,elemasymp[6],Xend);
  vasymp=vnorm_c(Xend+3);
  fprintf(stdout,"Asymptotic velocity: %e\n",vasymp);
  during=truncation/(vasymp*1e3);
  fprintf(stdout,"Estimated time of ingress: %e\n",during/YEAR);
  conics_c(elemasymp,elemasymp[6]-during,Xend);
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
  
  ////////////////////////////////////////////////////
  //COMPUTE POSITION AT INGRESS
  ////////////////////////////////////////////////////
  FILE* fc=fopen("wanderer.csv","w");
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
    conics_c(elts,elts[6]-during,Xu);
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
    fprintf(fc,"%s",vec2strn(elements,8,"%.17e,"));
    //ASYMPTOTIC ELEMENTS
    fprintf(fc,"%s",vec2strn(elements,8,"%.17e,"));
    //TIME OF INGRESS
    fprintf(fc,"%.17e,",elts[6]-during);
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
}
