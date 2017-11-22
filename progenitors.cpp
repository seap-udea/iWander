#include <iwander.cpp>
using namespace std;

#define VERBOSE 1

int main(int argc,char* argv[])
{
  /*
    Example: ./progenitors.exe

    Function: select the potential progenitors from a list of
    candidates.

    Input:
    * candidates.csv

    Output: 

    * encounters.csv
      Cols:
          0: n
          1-6: postar (present)
	  7-9: pos.body periastro
	  10-12: pos.star periastro
	  13: dmin (pc)
	  14: tmin1 (yr)
	  14: tmin2 (yr)
	  15-17: relative velocity (km/s)
	  18: relative speed (km/s)

    * candidates.csv
      Cols:
          0: n
          1-6: postar (present)
	  7-9: pos.body periastro
	  10-12: pos.star periastro
	  13: dmin (pc)
	  14: tmin1 (yr)
	  14: tmin2 (yr)
	  15-17: relative velocity (km/s)
	  18: relative speed (km/s)
	  19-??:: All the fields in the AstroRV catalogue
  */

  ////////////////////////////////////////////////////
  //CONFIGURATION
  ////////////////////////////////////////////////////
  #include <progenitors.conf>

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initWander();

  ////////////////////////////////////////////////////
  //UNITS
  ////////////////////////////////////////////////////
  VPRINT(stdout,"Units:\n\tUM = %.5e kg=%.5e Msun\n\tUL = %.17e m = %.17e pc\n\tUT = %.5e s = %.5e yr\n\tUV = %.5e m/s = %.5e km/s\n\tG = %.5e\n",
	 UM,UM/MSUN,UL,UL/PARSEC,UT,UT/YEAR,UV,UV/1e3,GGLOBAL);
  VPRINT(stdout,"\n");

  ////////////////////////////////////////////////////
  //GENERAL VARIABLES
  ////////////////////////////////////////////////////
  double tmp;
  char ctmp[100],line[MAXLINE];
  int Ntimes,Ntimesp,nsys,nsysp; 
  int ip,n;
  int nfields;
  double params[10];
  double hstep,duration;
  //MATRICES WITH INTEGRATIONS
  double *xIntp0,**xIntp,**xIntc;
  double *xInt0,**xInt;
  double *tsp,*ts;
  //INITIAL CONDITIONS
  double *dxIntdt,*x,*xg,*xp1,*xp2,*xpmin,*dx,*x0;
  double dmin,tmin,ftmin,dyn_tmin,dyn_dmin,dyn_vrel;
  double G;
  double t;
  int it;

  ////////////////////////////////////////////////////
  //GLOBAL ALLOCATION
  ////////////////////////////////////////////////////
  nsysp=6*Npart;
  nsys=6;

  Ntimesp=10000;
  Ntimes=100;

  x=(double*)malloc(6*sizeof(double));//LSR STATE VECTOR
  xg=(double*)malloc(6*sizeof(double));//GGLOBALC STATE VECTOR
  dx=(double*)malloc(6*sizeof(double));//LSR STATE VECTOR
  xpmin=(double*)malloc(6*sizeof(double));//GGLOBALC STATE VECTOR

  xIntp0=(double*)malloc(nsysp*sizeof(double));
  xInt0=(double*)malloc(nsys*sizeof(double));

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

  ////////////////////////////////////////////////////
  //GGLOBALLOBAL PROPERTIES FOR INTEGGLOBALRATION
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
  
  
  return 0;
}
