#include <gravray.cpp>
using namespace std;

#define VERBOSE 0

double uniform(double a=0,double b=1)
{
  double u=gsl_rng_uniform(RAND);
  return a+(b-a)*u;
}

double cart2sph(double x,double y,double z,double *lat,double *lon)
{
  double rho=sqrt(x*x+y*y);
  double r=sqrt(rho*rho+z*z);
  double psi=PI/2-acos(z/r);
  double phi;
  double phir=atan(y/x);
  if(x>0) phi=phir;
  else if(y>0) phi=phir+PI;
  else phi=phir-PI;
  *lat=psi*180/PI;
  *lon=phi*180/PI;
  return 0;
}

/*
int velocityDistribution(double Ms=1.0,double ap=1.0,
			 double Mp=1e-3,double Rk=7e7,
			 double *vmean,double *vstd)
*/
int velocityDistribution(double Ms,double ap,
			 double Mp,double Rk,
			 double *vm,double *vs)
{
  //////////////////////////////////////////////////////////////
  //UNITS AND CONSTANTS
  //////////////////////////////////////////////////////////////
  double G=1.0;
  double UL=1*AU;
  double UM=1*MSUN;
  double UT=sqrt(G*UL*UL*UL/(GCONST*UM));
  double UV=UL/UT;

  double TOLMEAN=1e-5;
  double TOLSTD=1e-3;
  double TOLRATIO=1e-3;

  //////////////////////////////////////////////////////////////
  //SYSTEM PROPERTIES
  //////////////////////////////////////////////////////////////
  double Rp=Rk/UL;
  double mu=G*Ms;
  double fh=0.1;
  double RH=fh*ap*pow(Mp/(3*Ms),1./3);
  double vpm=sqrt(mu/ap);
  double vesc=sqrt(2*mu/ap);
  double mup=G*Mp;

  double vpx=-vpm;
  double vpy=0;
  double vpz=0;
  double vp[]={vpx,vpy,vpz};

  //#############################################################
  //#GENERATE RANDOM ELEMENTS
  //#############################################################
  int Npart=5000;
  int n=0;
  double vinfs[Npart];

  double vmean,vstd;
  double vmean_old=0;
  double vstd_old=0;
  double ratio_old=0;

  FILE *fv=fopen("vinfs.data","w");
  while(n<Npart){
    //#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //#GENERATE RANDOM ASTROCENTRIC VELOCITIES FOR TEST PARTICLES
    //#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    double coswb=2;
    double ab,eb,pb,ib,hop,Ob,wb,wpf;
    while(fabs(coswb)>1){
      ab=uniform(0.5*ap,1.5*ap);
      eb=uniform(0.0,1.0);
      pb=ab*(1-eb*eb);
      ib=uniform(0,90);
      hop=sqrt(mu/pb);
      Ob=0.0;
      if(uniform()>0.5) Ob=180.0;
      if(Ob==0){
	coswb=(pb-ap)/(ap*eb);
	wpf=0.0;
      }else{
	coswb=(ap-pb)/(ap*eb);
	wpf=180.0;
      }
    }
    wb=acos(coswb)*RAD;
    double vi=sqrt(2*mu/ap-mu/ab);
    double xdot=-hop*(cos(Ob*DEG)*(eb*sin(wb*DEG)+sin(wpf*DEG))+
		      sin(Ob*DEG)*cos(ib*DEG)*(eb*cos(wb*DEG)+cos(wpf*DEG)));
    double ydot=-hop*(sin(Ob*DEG)*(eb*sin(wb*DEG)+cos(wpf*DEG))-
		      cos(Ob*DEG)*cos(ib*DEG)*(eb*cos(wb*DEG)+cos(wpf*DEG)));
    double zdot=+hop*sin(ib*DEG)*(eb*cos(wb*DEG)+
				  cos(wpf*DEG));
  
    double xdotrel=(-ydot)-vpx;
    double ydotrel=(+xdot)-vpy;
    double zdotrel=(+zdot)-vpz;
    double Vi[]={xdotrel,ydotrel,zdotrel};
    double Vin=vnorm_c(Vi);
    double vr=vnorm_c(Vi);
  
    //#PLANETOCENTRIC DIRECTION
    double psi,phi;
    cart2sph(xdotrel,ydotrel,zdotrel,&psi,&phi);
  
    //#GENERATE RANDOM INCOMING IMPACT PARAMETER
    double xtp=uniform(-RH,RH);
    double ytp=uniform(-RH,RH);
    double beta=asin(xtp/RH)*RAD,csi=asin(ytp/RH)*RAD;

    //#INCOMING POSITION
    double Ri[]={-RH*sin((phi+beta)*DEG)*cos((psi+csi)*DEG),
		 +RH*cos((phi+beta)*DEG)*cos((psi+csi)*DEG),
		 +RH*sin((psi+csi)*DEG)};

    //#COMPUTE U
    double U[3],u[3],Umag;
    vcrss_c(Ri,Vi,U);
    unorm_c(U,u,&Umag);

    //#IS THE PLANET APPROACHING?
    double qap=vdot_c(Ri,Vi);

    if(qap>0) continue;

    //#IMPACT PARAMETER
    double B=Umag/vnorm_c(Vi);

    if(B<1.1*Rp) continue;


    //#GAMMA
    double gamma=2*atan(mup/(B*Vin*Vin))*RAD;

    //#PLANETOCENTRIC ECCENTRICITY
    double e=1/sin(gamma*DEG/2);

    //#PERICENTER
    double q=mup*(e-1)/(Vin*Vin);

    if(q<1.1*Rp) continue;

    //#ROTATION OF INCOMING VELOCITY VECTOR
    double c=cos(gamma*DEG),s=sin(gamma*DEG);
    double ux=u[0],uy=u[1],uz=u[2];
    double M[][3]={{c+(1-c)*ux*ux,(1-c)*uy*ux-s*uz,(1-c)*uz*ux+s*uy},
		   {(1-c)*ux*uy+s*uz,c+(1-c)*uy*uy,(1-c)*ux*uy-s*ux},
		   {(1-c)*ux*uz-s*uy,(1-c)*uy*uz+s*ux,c+(1-c)*uz*uz}};
    double Vf[3];
    mxv_c(M,Vi,Vf);

    //#ASTROCENTRIC OUTBOUND VELOCITY
    double vf[3];
    vadd_c(Vf,vp,vf);
    double vfn=vnorm_c(vf);

    //#CHECK IF OBJECT IS BOUND
    if(vfn<vesc) continue;

    //#INFINITE VELOCITY
    double vinf=sqrt(vfn*vfn-vesc*vesc);

    VPRINT(stdout,"n= %d, gamma = %e, e = %e, q = %e, Rp = %e, vfn = %e, vesc = %e, vinf = %e\n",n,gamma,e,q,Rp,vfn,vesc,vinf*UV/1e3);

    vinfs[n]=vinf;//*UV/1e3;
    n++;

    double vmean=gsl_stats_mean(vinfs,1,n);
    double vstd=sqrt(gsl_stats_variance(vinfs,1,n));
    double ratio=vstd/vmean;

    //VPRINT(stdout,"%s\n",vec2strn(vinfs,n,"%.3e "));
    VPRINT(stdout,"vmean = %e, vstd = %e, ratio = %e\n",vmean,vstd,ratio);

    if(fabs(vmean-vmean_old)/vmean<TOLMEAN && 
       fabs(vstd-vstd_old)/vstd<TOLSTD && 
       fabs(ratio-ratio_old)/ratio<TOLRATIO)
      break;

    fprintf(fv,"%e\n",vinf);//*UV/1e3);

    vmean_old=vmean;
    vstd_old=vstd;
    ratio_old=ratio;
  }
  *vm=vmean_old;
  *vs=vstd_old;
  fclose(fv);
  return n;
}

int main(int argc,char* argv[])
{
  ////////////////////////////////////////////////////
  //INITIALIZE SPICE
  ////////////////////////////////////////////////////
  initSpice();

  ////////////////////////////////////////////////////
  //DECLARE
  ////////////////////////////////////////////////////
  int navg;
  double vmean,vstd;

  ////////////////////////////////////////////////////
  //PARAMETERS
  ////////////////////////////////////////////////////
  double Ms=1.0;
  double ap=5.0;
  double Mp=1.0e-3,logMp;
  double Rp=7e5;//km

  /*
  navg=velocityDistribution(Ms,ap,Mp,Rp,&vmean,&vstd);
  fprintf(stdout,"Converge after %d bodies\n",navg);
  fprintf(stdout,"vmean = %e, vstd = %e, ratio = %e\n",vmean,vstd,vstd/vmean);
  exit(0);
  //*/

  double ap1=0.4,ap2=5.0,dap=0.2;
  double logMp1=log10(1e-4),logMp2=log10(5e-3),dlogMp=(logMp2-logMp1)/20.0;

  FILE* fe=fopen("ejection.data","w");
  int n=1;
  for(ap=0.5;ap<=5.0+dap;ap+=dap){
    for(logMp=logMp1;logMp<=logMp2+dlogMp;logMp+=dlogMp){
      Mp=pow(10.0,logMp);
      navg=velocityDistribution(Ms,ap,Mp,Rp,&vmean,&vstd);
      fprintf(stdout,"n = %d: ap = %lf, Mp = %e\n",n,ap,Mp);
      fprintf(fe,"%e %e %e %e\n",ap,Mp,vmean,vstd);
      n++;
    }
  }
  return 0;
}
