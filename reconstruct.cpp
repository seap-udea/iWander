#include <iwander.cpp>
using namespace std;

#define VERBOSE 2 //Verbosity level
#define OSTREAM stdout //Stream where the output is redirected
#define VSTREAM stderr //Stream where the error output is redirected

int main(int argc,char* argv[])
{
  /*
    Example: ./reconstruct.exe
    
    Function: 

      This program try to reconstruct the trajectory of the wanderer
      among the nearby stars and the progenitor candidates.
      
    Input:
    * wanderer-<object>.csv
    * encounters-<object>.csv
    * progenitors-<object>.csv

    Output: 

    * simstars-<object>.csv
    * simulation-<object>.csv
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
  int NFIELDS;
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
  char **FIELDS=charMatrixAllocate(MAXCOLS,MAXTEXT);
  double *xInt0=(double*)malloc(6*Nobjs*sizeof(double));
  double *ts=(double*)malloc(Ntimes*sizeof(double));
  double **xInt=matrixAllocate(Ntimes,6*Nobjs);
  k=0;
  sprintf(FILENAME,"simulation-%s.csv",Wanderer);
  fs=fopen(FILENAME,"w");
  fprintf(fs,"t,");

  sprintf(FILENAME,"simstars-%s.csv",Wanderer);
  fe=fopen(FILENAME,"w");

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
  print0(OSTREAM,"Params = %s\n",vec2strn(params,10,"%e "));

  ////////////////////////////////////////////////////
  //SOLAR PROPERTIES
  ////////////////////////////////////////////////////
  copyVec(x,xnull,6);
  LSR2GC(x,xgc);
  print1(VSTREAM,"\tSolar position (cartesian): %s\n",vec2strn(xgc,6,"%e "));
  vscl_c(1e3/UL,xgc,xgc);//SET UNITS
  vscl_c(1e3/UV,xgc+3,xgc+3);
  cart2polar(xgc,xpgc);
  print1(VSTREAM,"\tSolar position (polar): %s\n",vec2strn(xpgc,6,"%e "));
  copyVec(xInt0+6*k,xpgc,6);
  for(int i=0;i<6;i++) fprintf(fs,"x%d_%d,",k,i);
  k++;

  ////////////////////////////////////////////////////
  //READ SURROGATE OBJECTS
  ////////////////////////////////////////////////////
  printHeader(stdout,"READING SURROGATE OBJECTS");
  print0(OSTREAM,"Object: %s\n",Wanderer_Name);
  sprintf(FILENAME,"wanderer-%s.csv",Wanderer);
  fc=fopen(FILENAME,"r");
  fgets(LINE,MAXLINE,fc);//HEADER
  i=0;
  while(fgets(LINE,MAXLINE,fc)!=NULL){
    parseLine(LINE,FIELDS,&NFIELDS);
    ting=atof(FIELDS[Wanderer::TING])/UT;
    n=Wanderer::XGAL;
    for(int k=0;k<6;k++) x[k]=atof(FIELDS[n++]);
    LSR2GC(x,xgc);
    vscl_c(1e3/UL,xgc,xgc);//SET UNITS
    vscl_c(1e3/UV,xgc+3,xgc+3);
    cart2polar(xgc,xpgc);
    if(i==0){
      print1(VSTREAM,"Nominal (cartesian): %s\n",vec2strn(x,6,"%e "));
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
  print0(OSTREAM,"Ingress time: %e year\n",ting*UT/YEAR);
  print0(OSTREAM,"Nominal position (pc,rad,pc,3*pc/year): %s\n",vec2strn(xnom,6,"%e "));
  print0(OSTREAM,"Farthest object (pc,rad,pc,3*pc/year): %s\n",vec2strn(xmax,6,"%e "));
  print0(OSTREAM,"Maximum distance: %e pc\n",dmax);

  ////////////////////////////////////////////////////
  //READ PROGENITOR CANDIDATES
  ////////////////////////////////////////////////////
  printHeader(stdout,"READING PROGENITOR CANDIDATES");
  sprintf(FILENAME,"progenitors-%s.csv",Wanderer);
  fc=fopen(FILENAME,"r");
  fgets(LINE,MAXLINE,fc);//HEADER
  fprintf(fe,"%s",line);

  i=0;
  tminmax=0;
  while(fgets(LINE,MAXLINE,fc)!=NULL){
    parseLine(LINE,FIELDS,&NFIELDS);

    print1(VSTREAM,"\tProgenitor %d (%s,%s,%s):\n",i,
	   FIELDS[Progenitors::HIP],FIELDS[Progenitors::TYCHO2_ID],FIELDS[Progenitors::NAME_SIMBAD]);
    tmin=atof(FIELDS[Progenitors::NOMTMIN]);
    dmin=atof(FIELDS[Progenitors::NOMDMIN]);

    //if(strcmp(FIELDS[Progenitors::NAME_SIMBAD],"HD_200325")!=0) continue;

    print1(VSTREAM,"\t\tStar : %s,%s,%s\n",
	   FIELDS[Progenitors::HIP],
	   FIELDS[Progenitors::TYCHO2_ID],
	   FIELDS[Progenitors::NAME_SIMBAD]);

    for(j=0;j<Progenitors::CAT;j++)
      fprintf(fe,"%s,",FIELDS[j]);
    fprintf(fe,"%s",FIELDS[j]);

    if(fabs(tmin)>fabs(tminmax)){tminmax=tmin;}
    print1(VSTREAM,"\t\ttmin = %e, dmin = %e\n",tmin,dmin);

    n=Progenitors::POSTARX;
    for(j=0;j<6;j++) x[j]=atof(FIELDS[n++]);
    vscl_c(UL/1e3,x,x);//SET UNITS
    print1(VSTREAM,"\t\tProgenitor %d (cartesian): %s\n",i,vec2strn(x,6,"%e "));

    LSR2GC(x,xgc);
    vscl_c(1e3/UL,xgc,xgc);//SET UNITS
    vscl_c(1e3/UV,xgc+3,xgc+3);
    cart2polar(xgc,xpgc);
    print1(VSTREAM,"\t\tProgenitor %d (polar): %s\n",i,vec2strn(xpgc,6,"%e "));

    //INTEGRATE BACK TO TIME OF INGRESS
    hstep=fabs(ting)/10;
    print1(VSTREAM,"\t\tIntegrating ting = %e, hstep = %e...\n",ting,hstep);
    params[0]=6;
    integrateEoM(0,xpgc,hstep,2,ting,6,EoMGalactic,params,ts,xInt);
    print1(VSTREAM,"\t\tTimes: %s\n",vec2strn(ts,2,"%e "));
    copyVec(xpgc,xInt[1],6);
    print1(VSTREAM,"\t\tProgenitor %d (polar): %s\n",i,vec2strn(xpgc,6,"%e "));
    //exit(0);

    //COPY TO INITIAL CONDITIONS FILE
    copyVec(xInt0+6*k,xpgc,6);
    for(j=0;j<6;j++) fprintf(fs,"x%d_%d,",k,j);

    k++;
    i++;
  }
  fclose(fe);
  Nprog=i;
  print0(OSTREAM,"Number of progenitors: %d\n",Nprog);
  print0(OSTREAM,"Maximum tmin: %e years\n",tminmax);

  Nobjs=k;
  print0(OSTREAM,"Number of objects: %d\n",Nobjs);
  fprintf(fs,"dummy\n");

  ////////////////////////////////////////////////////
  //INTEGRATE
  ////////////////////////////////////////////////////
  printHeader(stdout,"INTEGRATING OBJECTS");
  hstep=fabs(tminmax)/1000.0;
  print1(VSTREAM,"Initial conditions: %s...\n",vec2strn(xInt0,6*4,"%e "));
    
  print0(OSTREAM,"hstep = %e\n",hstep);
  params[0]=6*Nobjs;
  integrateEoM(0,xInt0,hstep,Ntimes,1.1*tminmax,6*Nobjs,EoMGalactic,params,ts,xInt);

  print1(VSTREAM,"Integration times (years): %s...\n",vec2strn(ts,10,"%e "));
  print1(VSTREAM,"Results (years): %s...\n",vec2strn(xInt[Ntimes-1],6*3,"%e "));

  ////////////////////////////////////////////////////
  //ANALYSE AND SAVE POSITIONS
  ////////////////////////////////////////////////////
  sprintf(FILENAME,"simulation-%s.dat",Wanderer);
  fc=fopen(FILENAME,"w");
  printHeader(stdout,"ANALYSING AND SAVING POSITIONS");
  for(int i=0;i<Ntimes;i++){
    print1(VSTREAM,"t = %e yr:\n",ts[i]);
    fprintf(fs,"%e,",ts[i]);

    k=0;
    //SOLAR POSITION
    copyVec(xpgc,xInt[i]+6*k,6);k++;
    polar2cart(xpgc,xsun);
    print1(VSTREAM,"Solar position: %s\n",vec2strn(xsun,6,"%e "));
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
    print1(VSTREAM,"\tRelative position of object: %s\n",vec2strn(xrel,6,"%e "));

    polar2cart(xgc,x);

    //DISTANCE BETWEEN OBJECTS
    vsub_c(x,xnom,dx);
    dnom=vnorm_c(dx);
    print1(VSTREAM,"\tDistance between test particles = %e pc\n",ts[i],dnom);

    //DISTANCE TO THE SUN
    vsub_c(xsun,xnom,dx);
    dsun=vnorm_c(dx);
    print1(VSTREAM,"\tDistance to the Sun = %e pc\n",dsun);

    print1(VSTREAM,"\tDistance to the Stars (Nstar = %d):\n",Nobjs-3);
    for(int j=0;j<Nobjs-3;j++){
      //STARS POSITIONS
      print1(VSTREAM,"\t\tStar %d:\n",j);
      copyVec(xpgc,xInt[i]+6*k,6);k++;
      polar2cart(xpgc,x);

      vsubg_c(x,xsun,6,xrel);
      vscl_c(UV/1e3,xrel+3,xrel+3);
      print1(VSTREAM,"\t\tRelative position: %s\n",vec2strn(xrel,6,"%e "));

      fprintf(fs,"%s",vec2strn(xrel,6,"%e,"));

      vsub_c(xnom,x,dx);
      dstar=vnorm_c(dx);
      print1(VSTREAM,"\t\t\tDistance object-star %d = %e pc\n",j,dstar);
      vsub_c(xsun,x,dx);
      dstar=vnorm_c(dx);
      print1(VSTREAM,"\t\t\tDistance sun-star %d = %e pc\n",j,dstar);
    }
    fprintf(fs,"\n");
  }

  ////////////////////////////////////////////////////
  //FINALIZE
  ////////////////////////////////////////////////////
  fclose(fs);

  return 0;
}
