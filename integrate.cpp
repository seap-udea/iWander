#include <iwander.cpp>
using namespace std;
#define VERBOSE 1

int main(int argc,char* argv[])
{
  /*
    Example: ./integrate.exe

    Input: 
    * wanderer.csv: test particles compatible with A/2017U1 orbit.

    Output:
    * cloud.csv
  */

  ////////////////////////////////////////////////////
  //CONFIGURATION
  ////////////////////////////////////////////////////
  #include <integrate.conf>

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initWander();
  
  ////////////////////////////////////////////////////
  //GENERAL VARIABLES
  ////////////////////////////////////////////////////
  double tmp;
  char ctmp[100],line[MAXLINE];
  int nsys,nsysp; 
  int ip,n;
  int nfields;
  double params[10];
  //MATRICES WITH INTEGRATIONS
  double *xIntp0,**xIntp,**xIntc;
  double *xInt0,**xInt;
  double *tsp,*ts;
  //INITIAL CONDITIONS
  double *dxIntdt,*x,*xg,*x0;
  double G;

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

  x=(double*)malloc(6*sizeof(double));//LSR STATE VECTOR
  xg=(double*)malloc(6*sizeof(double));//GC STATE VECTOR

  xIntp0=(double*)malloc(nsysp*sizeof(double));
  xInt0=(double*)malloc(nsys*sizeof(double));

  xIntp=(double**)malloc(Ntimesp*sizeof(double*));
  for(int j=0;j<Ntimesp;j++) xIntp[j]=(double*)malloc(nsysp*sizeof(double));

  xIntc=(double**)malloc(Ntimesp*sizeof(double*));
  for(int j=0;j<Ntimesp;j++) xIntc[j]=(double*)malloc(nsysp*sizeof(double));
  
  char **fields=(char**)malloc(MAXCOLS*sizeof(char*));
  for(int i=0;i<MAXCOLS;i++) fields[i]=(char*)malloc(MAXTEXT*sizeof(char));

  tsp=(double*)malloc(Ntimesp*sizeof(double));
  ts=(double*)malloc(Ntimesp*sizeof(double));

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
  //READING TEST PARTICLES
  ////////////////////////////////////////////////////
  VPRINT(stdout,"Reading %d test particles\n",Npart);
  VPRINT(stdout,"\n");
  FILE *fc=fopen("wanderer.csv","r");
  fgets(line,MAXLINE,fc);//HEADER

  int i=0;
  while(fgets(line,MAXLINE,fc)!=NULL){
    parseLine(line,fields,&nfields);
    ip=6*i;
    n=Wanderer::XGAL;
    for(int k=0;k<6;k++) x[k]=atof(fields[n++]);
    //TRANSFORM COORDINATES FROM LSR->GC
    LSR2GC(x,xg);
    vscl_c(1e3/UL,xg,xg);//SET UNITS
    vscl_c(1e3/UV,xg+3,xg+3);
    //CONVERT TO CYLINDRICAL GALACTIC COORDINATES
    cart2polar(xg,xIntp0+ip,1.0);
    i++;
    if(i==Npart) break;
  }
  fclose(fc);

  ////////////////////////////////////////////////////
  //INTEGRATING TEST PARTICLES FOR A LONG PERIOD 
  ////////////////////////////////////////////////////
  //INTEGRATE
  VPRINT(stdout,"Integrating %d test particles (nsys = %d)\n",Npart,nsysp);
  VPRINT(stdout,"\n");

  params[0]=nsysp;
  integrateEoM(0,xIntp0,hstep,Ntimesp,duration,
	       nsysp,EoMGalactic,params,
	       tsp,xIntp);
  
  ////////////////////////////////////////////////////
  //SAVING POSITIONS
  ////////////////////////////////////////////////////
  VPRINT(stdout,"Saving %d test particles\n",Npart);
  VPRINT(stdout,"\n");

  fc=fopen("cloud.csv","w");
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //HEADER
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf(fc,"t,");
  for(int j=0;j<Npart;j++){
    fprintf(fc,"part%d-R,part%d-phi,part%d-Z,part%d-vR,part%d-dphi,part%d-vZ,",j,j,j,j,j,j);
    fprintf(fc,"part%d-x,part%d-y,part%d-z,part%d-vx,part%d-xy,part%d-vz,",j,j,j,j,j,j);
  }
  fprintf(fc,"dummy\n");

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //SAVING PARTICLE POSITIONS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for(int i=0;i<Ntimesp;i++){
    fprintf(fc,"%.5e,",tsp[i]);
    for(int j=0;j<Npart;j++){
      ip=6*j;
      fprintf(fc,"%s",vec2strn(xIntp[i]+ip,6,"%.17e,"));
      polar2cart(xIntp[i]+ip,xIntc[i]+ip,1.0);
      fprintf(fc,"%s",vec2strn(xIntc[i]+ip,6,"%.17e,"));
    }
    fprintf(fc,"\n");
  }
  fclose(fc);
  
  VPRINT(stdout,"Done.\n");
  return 0;
}
