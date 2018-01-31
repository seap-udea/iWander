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

    Function: This program compute the position of the interstellar
    objects and their surrogates at the time when the object entered
    the Solar System.

    Input: None

    Output: 

    * wanderer.data
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
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initWander();
  
  ////////////////////////////////////////////////////
  //GLOBAL DECLARATIONS
  ////////////////////////////////////////////////////
  //GRAVITATIONAL PARAMETER
  double n=ini_n*DEG/DAY;
  double a=ini_a*(-AU/1E3);
  double mu=n*n*a*a*a;
  //EPHEMERIS TIME
  double dt;
  double to=unitim_c(ini_to_jed,"JDTDB","TDB");
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
  //INTEGRATION
  int npoints=2;
  double tini=to;
  double direction=duration/fabs(duration);
  double params[]={6};
  duration/=UT;
  double h,t_start,t_step,tend,t_stop,t;
  double h_used,h_next,h_adjust,deltat;
  int i,status;
  //INITIAL CONDITIONS
  double X0[6],X[6],Xu[6],as,tdyn;
  //ELEMENTS
  double elts[8];

  ////////////////////////////////////////////////////
  //PREPARE ELEMENTS GENERATION
  ////////////////////////////////////////////////////
  SpiceDouble elements[6];
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
  //LOOP OVER PARTICLES
  ////////////////////////////////////////////////////
  pxform_c("ECLIPJ2000","GALACTIC",t,M_Eclip_Galactic);
  pxform_c("ECLIPJ2000","J2000",0,M_Eclip_J2000);

  FILE* fc=fopen("wanderer.csv","w");
  int Nfreq=ceil(Npart/10);
  Nfreq=Nfreq==0?1:Nfreq;

  //FILE HEADER
  fprintf(fc,"i,qo,eo,inco,Wo,wo,Moo,too,muo,t,xecl,yecl,zecl,vxecl,vyecl,vzecl,xsky,ysky,zsky,vxsky,vysky,vzsky,xgal,ygal,zgal,vxgal,vygal,vzgal,ra,dec,l,b,d,q,e,inc,W,w,Mo,to,mu,dummy\n");
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
    if(VERBOSE) fprintf(stdout,"\tInitial position @ Sun (EJ2000) : %s\n",
			vec2strn(position,3,"%.17e "));
    if(VERBOSE) fprintf(stdout,"\tInitial velocity @ Sun (EJ2000) : %s\n",
			vec2strn(position+3,3,"%.17e "));
    
    spkezr_c("SUN",to,"ECLIPJ2000","NONE",SSB,sun,&ltmp);
    vaddg_c(position,sun,6,position);
    if(VERBOSE) fprintf(stdout,"\tInitial position @ SSB (EJ2000) : %s\n",
	    vec2strn(position,3,"%.17e "));
    if(VERBOSE) fprintf(stdout,"\tInitial velocity @ SSB (EJ2000) : %s\n",
	    vec2strn(position+3,3,"%.17e "));

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INITIAL CONDITIONS
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    vpack_c(position[0]*1E3/UL,position[1]*1E3/UL,position[2]*1E3/UL,X0);
    vpack_c(position[3]*1E3/UV,position[4]*1E3/UV,position[5]*1E3/UV,X0+3);
    as=vnorm_c(X0);
    tdyn=2*M_PI*sqrt(as*as*as/(GGLOBAL*MSUN/UM));
    h=direction*tdyn/1000.0;
    t_start=tini/UT;
    t_step=duration/(npoints-1);
    tend=t_start+duration;
    t_stop=tend;
    t=t_start;

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INTEGRATION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    h=tdyn/1000.0;
    integrateEoM(tini/UT,X0,h,npoints,duration,6,EoM,params,ts,Xout);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //FINAL POSITION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //vscl_c(UL/1E3,X0,Xu);vscl_c(UV/1E3,X0+3,Xu+3);
    vscl_c(UL/1E3,Xout[1],Xu);vscl_c(UV/1E3,Xout[1]+3,Xu+3);
    if(VERBOSE) fprintf(stdout,"\tFinal position : %s\n",
	    vec2strn(Xu,3,"%.17e "));
    double d=vnorm_c(Xu)*1e3/AU;
    //fprintf(stdout,"\tDistance : %e\n",d);

    //ASYMPTOTIC ELEMENTS
    oscelt_c(Xu,t,MUTOT,elts);
    
    //J2000
    mxv_c(M_Eclip_J2000,Xu,posJ2000);
    mxv_c(M_Eclip_J2000,Xu+3,posJ2000+3);
    recrad_c(posJ2000,&d,&RA,&DEC);
    if(VERBOSE) fprintf(stdout,"\tRA(+HH:MM:SS) = %s, DEC(DD:MM:SS) = %s, d = %.3lf AU\n",
	    dec2sex(RA*RAD/15.0),dec2sex(DEC*RAD),d*1e3/AU);

    //GALACTIC
    mxv_c(M_Eclip_Galactic,Xu,posGalactic);
    mxv_c(M_Eclip_Galactic,Xu+3,posGalactic+3);
    recrad_c(posGalactic,&d,&l,&b);
    if(VERBOSE) fprintf(stdout,"\tl(+DD:MM:SS) = %s, b(DD:MM:SS) = %s\n",
	    dec2sex(l*RAD),dec2sex(b*RAD));

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //SAVE POSITION
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //ELEMENTOS
    fprintf(fc,"%d,%s",j,vec2strn(elements,8,"%.17e,"));
    //PARTICLE TEND
    fprintf(fc,"%.17e,",(tend*UT-tini));
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
    //DISTANCE
    fprintf(fc,"%.17e,",d*1E3/AU);
    //ASYMPTOTIC ELEMENTS
    fprintf(fc,"%s",vec2strn(elts,8,"%.17e,"));
    fprintf(fc,"\n");
  }
  fclose(fc);
}
