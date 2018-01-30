#include <iwander.cpp>
using namespace std;

#define VERBOSE 0

int main(int argc,char* argv[])
{
  /*
    Example: ./reconstruct.exe
    
    Function: 

      This program try to reconstruct the trajectory of the wanderer
      among the nearby stars and the progenitor candidates.
      
    Input:
    * wanderer.csv
    * encounters.csv
    * progenitors.csv

    Output: 

    * simstars.csv
    * simulation.csv
  */

  ////////////////////////////////////////////////////
  //CONFIGURATION
  ////////////////////////////////////////////////////
  #include <iwander.conf>
  #include <reconstruct.conf>

  ////////////////////////////////////////////////////
  //INITIALIZE iWANDER
  ////////////////////////////////////////////////////
  initWander();

  ////////////////////////////////////////////////////
  //VARIABLES
  ////////////////////////////////////////////////////
  int i,j,k,n,ip;
  int Nprog;
  int nfields;
  FILE *fc,*fs,*fe;
  char ctmp[100],line[MAXLINE],values[MAXLINE];
  double x[6],xnom[6],xmax[6],dx[6],xgc[6],xpgc[6],xg[6],xnull[]={0,0,0,0,0,0};
  double xsun[6],xrel[6];
  double dmax,dnom,ting,dsun,dstar,dmin;
  double params[10];
  double tmin,tminmax,hstep;

  ////////////////////////////////////////////////////
  //ALLOCATION AND INITIALIZATION
  ////////////////////////////////////////////////////
  char **fields=charMatrixAllocate(MAXCOLS,MAXTEXT);
  double *xInt0=(double*)malloc(6*Nobjs*sizeof(double));
  double *ts=(double*)malloc(Ntimes*sizeof(double));
  double **xInt=matrixAllocate(Ntimes,6*Nobjs);
  k=0;
  sprintf(Filename,"simulation-%s.csv",WANDERER);
  fs=fopen(Filename,"w");
  fprintf(fs,"t,");

  sprintf(Filename,"simstars-%s.csv",WANDERER);
  fe=fopen(Filename,"w");

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
  //SOLAR PROPERTIES
  ////////////////////////////////////////////////////
  copyVec(x,xnull,6);
  LSR2GC(x,xgc);
  VPRINT(stdout,"\tSolar position (cartesian): %s\n",vec2strn(xgc,6,"%e "));
  vscl_c(1e3/UL,xgc,xgc);//SET UNITS
  vscl_c(1e3/UV,xgc+3,xgc+3);
  cart2polar(xgc,xpgc);
  VPRINT(stdout,"\tSolar position (polar): %s\n",vec2strn(xpgc,6,"%e "));
  copyVec(xInt0+6*k,xpgc,6);
  for(int i=0;i<6;i++) fprintf(fs,"x%d_%d,",k,i);
  k++;

  ////////////////////////////////////////////////////
  //READ SURROGATE OBJECTS
  ////////////////////////////////////////////////////
  printHeader(stdout,"READING SURROGATE OBJECTS");
  fprintf(stdout,"Object: %s\n",WANDERER_NAME);
  sprintf(Filename,"wanderer-%s.csv",WANDERER);
  fc=fopen(Filename,"r");
  fgets(line,MAXLINE,fc);//HEADER
  i=0;
  while(fgets(line,MAXLINE,fc)!=NULL){
    parseLine(line,fields,&nfields);
    ting=atof(fields[Wanderer::TING])/UT;
    n=Wanderer::XGAL;
    for(int k=0;k<6;k++) x[k]=atof(fields[n++]);
    LSR2GC(x,xgc);
    vscl_c(1e3/UL,xgc,xgc);//SET UNITS
    vscl_c(1e3/UV,xgc+3,xgc+3);
    cart2polar(xgc,xpgc);
    if(i==0){
      VPRINT(stdout,"Nominal (cartesian): %s\n",vec2strn(x,6,"%e "));
      copyVec(xnom,xpgc,6);
      copyVec(xInt0+6*k,xnom,6);
      for(int j=0;j<6;j++) fprintf(fs,"x%d_%d,",k,j);
      k++;
    }else{
      polar2cart(xnom,x);
      polar2cart(xpgc,xg);
      vsub_c(x,xg,dx);
      dnom=vnorm_c(dx);
      if(dnom>dmax){
	copyVec(xmax,xpgc,6);
	dmax=dnom;
      }
    }
    i++;
  }
  copyVec(xInt0+6*k,xmax,6);
  k++;
  fprintf(stdout,"Ingress time: %e year\n",ting*UT/YEAR);
  fprintf(stdout,"Nominal position (pc,rad,pc,3*pc/year): %s\n",vec2strn(xnom,6,"%e "));
  fprintf(stdout,"Farthest object (pc,rad,pc,3*pc/year): %s\n",vec2strn(xmax,6,"%e "));
  fprintf(stdout,"Maximum distance: %e pc\n",dmax);

  ////////////////////////////////////////////////////
  //READ PROGENITOR CANDIDATES
  ////////////////////////////////////////////////////
  printHeader(stdout,"READING PROGENITOR CANDIDATES");
  sprintf(Filename,"progenitors-%s.csv",WANDERER);
  fc=fopen(Filename,"r");
  fgets(line,MAXLINE,fc);//HEADER
  fprintf(fe,"%s",line);

  i=0;
  tminmax=0;
  while(fgets(line,MAXLINE,fc)!=NULL){
    parseLine(line,fields,&nfields);

    VPRINT(stdout,"\tProgenitor %d (%s,%s,%s):\n",i,
	   fields[Progenitors::HIP],fields[Progenitors::TYCHO2_ID],fields[Progenitors::NAME_SIMBAD]);
    tmin=atof(fields[Progenitors::NOMTMIN]);
    dmin=atof(fields[Progenitors::NOMDMIN]);

    //if(strcmp(fields[Progenitors::NAME_SIMBAD],"HD_200325")!=0) continue;

    VPRINT(stdout,"\t\tStar : %s,%s,%s\n",
	   fields[Progenitors::HIP],
	   fields[Progenitors::TYCHO2_ID],
	   fields[Progenitors::NAME_SIMBAD]);

    for(j=0;j<Progenitors::CAT;j++)
      fprintf(fe,"%s,",fields[j]);
    fprintf(fe,"%s",fields[j]);

    if(fabs(tmin)>fabs(tminmax)){tminmax=tmin;}
    VPRINT(stdout,"\t\ttmin = %e, dmin = %e\n",tmin,dmin);

    n=Progenitors::POSTARX;
    for(j=0;j<6;j++) x[j]=atof(fields[n++]);
    vscl_c(UL/1e3,x,x);//SET UNITS
    VPRINT(stdout,"\t\tProgenitor %d (cartesian): %s\n",i,vec2strn(x,6,"%e "));

    LSR2GC(x,xgc);
    vscl_c(1e3/UL,xgc,xgc);//SET UNITS
    vscl_c(1e3/UV,xgc+3,xgc+3);
    cart2polar(xgc,xpgc);
    VPRINT(stdout,"\t\tProgenitor %d (polar): %s\n",i,vec2strn(xpgc,6,"%e "));

    //INTEGRATE BACK TO TIME OF INGRESS
    hstep=fabs(ting)/10;
    VPRINT(stdout,"\t\tIntegrating ting = %e, hstep = %e...\n",ting,hstep);
    params[0]=6;
    integrateEoM(0,xpgc,hstep,2,ting,6,EoMGalactic,params,ts,xInt);
    VPRINT(stdout,"\t\tTimes: %s\n",vec2strn(ts,2,"%e "));
    copyVec(xpgc,xInt[1],6);
    VPRINT(stdout,"\t\tProgenitor %d (polar): %s\n",i,vec2strn(xpgc,6,"%e "));
    //exit(0);

    //COPY TO INITIAL CONDITIONS FILE
    copyVec(xInt0+6*k,xpgc,6);
    for(j=0;j<6;j++) fprintf(fs,"x%d_%d,",k,j);

    k++;
    i++;
  }
  fclose(fe);
  Nprog=i;
  fprintf(stdout,"Number of progenitors: %d\n",Nprog);
  fprintf(stdout,"Maximum tmin: %e years\n",tminmax);

  Nobjs=k;
  fprintf(stdout,"Number of objects: %d\n",Nobjs);
  fprintf(fs,"dummy\n");

  /*
    fprintf(stdout,"Initial conditions: %s\n",
    vec2strn(xInt0,6*Nobjs,"%e "));
  */

  ////////////////////////////////////////////////////
  //INTEGRATE
  ////////////////////////////////////////////////////
  printHeader(stdout,"INTEGRATING OBJECTS");
  hstep=fabs(tminmax)/1000.0;
  VPRINT(stdout,"Initial conditions: %s...\n",vec2strn(xInt0,6*4,"%e "));
    
  fprintf(stdout,"hstep = %e\n",hstep);
  params[0]=6*Nobjs;
  integrateEoM(0,xInt0,hstep,Ntimes,1.1*tminmax,6*Nobjs,EoMGalactic,params,ts,xInt);

  VPRINT(stdout,"Integration times (years): %s...\n",vec2strn(ts,10,"%e "));
  VPRINT(stdout,"Results (years): %s...\n",vec2strn(xInt[Ntimes-1],6*3,"%e "));

  ////////////////////////////////////////////////////
  //ANALYSE AND SAVE POSITIONS
  ////////////////////////////////////////////////////
  sprintf(Filename,"simulation-%s.dat",WANDERER);
  fc=fopen(Filename,"w");
  printHeader(stdout,"ANALYSING AND SAVING POSITIONS");
  for(int i=0;i<Ntimes;i++){
    VPRINT(stdout,"t = %e yr:\n",ts[i]);
    fprintf(fs,"%e,",ts[i]);

    k=0;
    //SOLAR POSITION
    copyVec(xpgc,xInt[i]+6*k,6);k++;
    polar2cart(xpgc,xsun);
    VPRINT(stdout,"Solar position: %s\n",vec2strn(xsun,6,"%e "));
    vscl_c(1,xsun,x);
    vscl_c(UV/1e3,xsun+3,x+3);
    fprintf(fs,"%s",vec2strn(x,6,"%e,"));

    //OBJECTS POSITION
    copyVec(xpgc,xInt[i]+6*k,6);k++;
    copyVec(xgc,xInt[i]+6*k,6);k++;

    polar2cart(xpgc,xnom);

    vsubg_c(xnom,xsun,6,xrel);
    vscl_c(UV/1e3,xrel+3,xrel+3);
    fprintf(fs,"%s",vec2strn(xrel,6,"%e,"));
    VPRINT(stdout,"\tRelative position of object: %s\n",vec2strn(xrel,6,"%e "));

    polar2cart(xgc,x);

    //DISTANCE BETWEEN OBJECTS
    vsub_c(x,xnom,dx);
    dnom=vnorm_c(dx);
    VPRINT(stdout,"\tDistance between test particles = %e pc\n",ts[i],dnom);

    //DISTANCE TO THE SUN
    vsub_c(xsun,xnom,dx);
    dsun=vnorm_c(dx);
    VPRINT(stdout,"\tDistance to the Sun = %e pc\n",dsun);

    VPRINT(stdout,"\tDistance to the Stars (Nstar = %d):\n",Nobjs-3);
    for(int j=0;j<Nobjs-3;j++){
      //STARS POSITIONS
      VPRINT(stdout,"\t\tStar %d:\n",j);
      copyVec(xpgc,xInt[i]+6*k,6);k++;
      polar2cart(xpgc,x);

      vsubg_c(x,xsun,6,xrel);
      vscl_c(UV/1e3,xrel+3,xrel+3);
      VPRINT(stdout,"\t\tRelative position: %s\n",vec2strn(xrel,6,"%e "));

      fprintf(fs,"%s",vec2strn(xrel,6,"%e,"));

      vsub_c(xnom,x,dx);
      dstar=vnorm_c(dx);
      VPRINT(stdout,"\t\t\tDistance object-star %d = %e pc\n",j,dstar);
      vsub_c(xsun,x,dx);
      dstar=vnorm_c(dx);
      VPRINT(stdout,"\t\t\tDistance sun-star %d = %e pc\n",j,dstar);
    }
    //getchar();
    fprintf(fs,"\n");
    //break;
  }
  exit(0);
  ////////////////////////////////////////////////////
  //FINALIZE
  ////////////////////////////////////////////////////
  fclose(fs);

  return 0;
}
