#include <iwander.cpp>
using namespace std;

#define VERBOSE 0

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
  int ip,n,np;
  int nfields;
  double params[10];
  double hstep,duration;
  //MATRICES WITH INTEGRATIONS
  double *xIntp0,**xIntp,**xIntc;
  double *xInt0,**xInt;
  double *tsp,*ts;
  double t1,t2;
  //INITIAL CONDITIONS
  double *dxIntdt,*x,*xg,*xp1,*xp2,*xpmin,*dx,*x0;
  double dmin,tmin,ftmin,dyn_tmin,dyn_dmin,dyn_vrel;
  double G;
  double t;
  int it;
  char values[MAXLINE];

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
  VPRINT(stdout,"Reading %d test particles\n",Npart);
  
  FILE *fc;
  fc=fopen("cloud.csv","r");
  fgets(line,MAXLINE,fc);//HEADER
  int i=0;
  while(fgets(line,MAXLINE,fc)!=NULL){
    parseLine(line,fields,&nfields);
    tsp[i]=atof(fields[0]);
    n=1;
    for(int j=0;j<Npart;j++){
      ip=6*j;
      x=xIntp[i]+ip;
      for(int k=0;k<6;k++) x[k]=atof(fields[n++]);
      x=xIntc[i]+ip;
      for(int k=0;k<6;k++) x[k]=atof(fields[n++]);
    }
    i++;
  }
  fclose(fc);
  
  ////////////////////////////////////////////////////
  //READING SELECTED OBJECTS
  ////////////////////////////////////////////////////
  params[0]=nsys;

  fc=fopen("candidates.csv","r");
  fgets(line,MAXLINE,fc);//HEADER

  FILE *fp=fopen("potential.csv","w");
  fprintf(fp,"dyntmin,dyndmin,dynvrel,%s",line);

  n=0;
  np=0;
  while(fgets(line,MAXLINE,fc)!=NULL){
    
    //PARSE FIELDS
    strcpy(values,line);
    parseLine(line,fields,&nfields);
    n++;

    fprintf(stdout,"Star %d,%s,%s,%s:\n",n,
	    fields[Candidates::HIP],
	    fields[Candidates::TYCHO2_ID],
	    fields[Candidates::NAME_SIMBAD]
	    );
    
    //LMA ENCOUNTER TIME
    duration=atof(fields[Candidates::TMIN])*YEAR/UT;
    fprintf(stdout,"\tLMA time:%e (compared to %e)\n",duration,tsp[Ntimesp-1]);
    dmin=atof(fields[Candidates::DMIN]);
    fprintf(stdout,"\tLMA distance:%e\n",dmin);
    hstep=fabs(duration)/(10*Ntimes);

    if(1.1*duration<tsp[Ntimesp-1]){
      fprintf(stderr,"\t\tEncounters happens out of the integration range (%e)...\n",
	      tsp[Ntimesp-1]);
      continue;
    }
    
    //SPATIAL COORDINATES AND LSR VELOCITY
    for(int k=0;k<6;k++) x[k]=atof(fields[Candidates::POSTARX+k]);
    vscl_c(PARSEC/1e3,x,x);

    //CONVERT TO GC
    LSR2GC(x,xg);
    vscl_c(1e3/UL,xg,xg);//SET UNITS
    vscl_c(1e3/UV,xg+3,xg+3);
    cart2polar(xg,xInt0,1.0);//CONVERSION TO CYLINDRICAL
    fprintf(stdout,"\tInitial position: %s\n",vec2strn(xInt0,6,"%.5e "));

    //INTEGRATE 
    duration*=1.1;
    fprintf(stdout,"\tIntegrating until t = %e, time step = %e\n",duration,hstep);

    integrateEoM(0,xInt0,hstep,Ntimes,duration,
		 nsys,EoMGalactic,params,
		 ts,xInt);

    //FIND TIME OF MINIMUM DISTANCE
    dyn_dmin=1.0e+100;
    for(int i=1;i<Ntimes;i++){

      //STAR POSITION 
      polar2cart(xInt[i],x,1.0);
      VPRINT(stdout,"\tPosition %d @ t=%e:%s\n",i,ts[i],vec2strn(x,6,"%.5e,"));

      //FIND TIME IN PARTICLE INTEGRATION
      t=ts[i];
      it=findTime(t,tsp,Ntimesp);
      VPRINT(stdout,"\t\tCorrespond to interval %d: [%e,%e]\n",it,tsp[it],tsp[it+1]);

      //FIND POSITION OF NOMINAL PARTICLE AT INITIAL POINT OF INTERVAL
      xp1=xIntc[it];
      xp2=xIntc[it+1];
      VPRINT(stdout,"\t\tNominal particle position 1 (%e):%s\n",tsp[it],vec2strn(xp1,6,"%.5e "));
      VPRINT(stdout,"\t\tNominal particle position 2 (%e):%s\n",tsp[it+1],vec2strn(xp2,6,"%.5e "));

      //DIFFERENCE
      vsubg_c(x,xp1,6,dx);
      VPRINT(stdout,"\t\tDifference 1: [%s] (pos:%.5e,vel:%.5e]\n",
	      vec2strn(dx,6,"%.5e,"),vnorm_c(dx),vnorm_c(dx+3)*UV);
      vsubg_c(x,xp2,6,dx);
      VPRINT(stdout,"\t\tDifference 2: [%s] (pos:%.5e,vel:%.5e]\n",
	      vec2strn(dx,6,"%.5e,"),vnorm_c(dx),vnorm_c(dx+3)*UV);

      //CALCULATE DISTANCE FROM PARTICLE TO INTERVAL
      dmin=distancePointLine(x,xp1,xp2,&ftmin);
      VPRINT(stdout,"\t\tPosition factor: %.5e\n",ftmin);
      VPRINT(stdout,"\t\tDistance at minimum: %.5e\n",dmin);
      
      //CORRECTED TIME
      t1=tsp[it]<tsp[it+1]?tsp[it]:tsp[it+1];
      t2=tsp[it]>tsp[it+1]?tsp[it]:tsp[it+1];
      tmin=t1+ftmin*(t2-t1);
      VPRINT(stdout,"\t\tCorrected time at minimum:%e\n",tmin);

      //INTERPOLATED POSITION OF NOMINAL PARTICLE
      vsubg_c(xp2,xp1,6,dx);
      vscl_c((tmin-tsp[it])/(tsp[it+1]-tsp[it]),dx,dx);
      vscl_c((tmin-tsp[it])/(tsp[it+1]-tsp[it]),dx+3,dx+3);
      vaddg_c(xp1,dx,6,xpmin);

      //DISTANCE TO STAR
      vsubg_c(xpmin,x,6,dx);
      VPRINT(stdout,"\t\tDifference at minimum: [%s] (pos:%.5e,vel:%.5e]\n",
	      vec2strn(dx,6,"%.5e,"),vnorm_c(dx),vnorm_c(dx+3)*UV/1e3);

      //COMPUTE MINIMUM DISTANCE
      if(dmin<=dyn_dmin){
	dyn_dmin=dmin;
	dyn_tmin=tmin;
	dyn_vrel=vnorm_c(dx+3)*UV/1e3;
      }
    }

    fprintf(stdout,"\tDynamical minimum distance: %e\n",dyn_dmin);
    fprintf(stdout,"\tDynamical minimum time: %e\n",dyn_tmin);
    
    //STORE THE BEST CANDIDATES
    if(dyn_dmin<dmax2){
      np++;
      fprintf(fp,"%.5lf,%.5lf,%.5lf,%s",dyn_tmin,dyn_dmin,dyn_vrel,values);
    }

    if(VERBOSE) break;
  }
  fclose(fc);
  fclose(fp);

  fprintf(stdout,"Number of potential progenitors identified: %d\n",np);

  return 0;
}
