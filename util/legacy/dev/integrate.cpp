#include <gravray.cpp>
using namespace std;

#define VERBOSE 0

//FIELDS OF THE CANDIDATES.CSV FILE
enum {N,POSTARX,POSTARY,POSTARZ,VELSTARX,VELSTARY,VELSTARZ,POSBODYPERIX,POSBODYPERIY,POSBODYPERIZ,POSTARPERIX,POSTARPERIY,POSTARPERIZ,DMIN,TMIN,VRELX,VRELY,VRELZ,VREL,HIP,TYCHO2_ID,SOLUTION_ID,SOURCE_ID,RANDOM_INDEX,REF_EPOCH,RA,RA_ERROR,DEC,DEC_ERROR,PARALLAX,PARALLAX_ERROR,PMRA,PMRA_ERROR,PMDEC,PMDEC_ERROR,RA_DEC_CORR,RA_PARALLAX_CORR,RA_PMRA_CORR,RA_PMDEC_CORR,DEC_PARALLAX_CORR,DEC_PMRA_CORR,DEC_PMDEC_CORR,PARALLAX_PMRA_CORR,PARALLAX_PMDEC_CORR,PMRA_PMDEC_CORR,ASTROMETRIC_N_OBS_AL,ASTROMETRIC_N_OBS_AC,ASTROMETRIC_N_GOOD_OBS_AL,ASTROMETRIC_N_GOOD_OBS_AC,ASTROMETRIC_N_BAD_OBS_AL,ASTROMETRIC_N_BAD_OBS_AC,ASTROMETRIC_DELTA_Q,ASTROMETRIC_EXCESS_NOISE,ASTROMETRIC_EXCESS_NOISE_SIG,ASTROMETRIC_PRIMARY_FLAG,ASTROMETRIC_RELEGATION_FACTOR,ASTROMETRIC_WEIGHT_AL,ASTROMETRIC_WEIGHT_AC,ASTROMETRIC_PRIORS_USED,MATCHED_OBSERVATIONS,DUPLICATED_SOURCE,SCAN_DIRECTION_STRENGTH_K1,SCAN_DIRECTION_STRENGTH_K2,SCAN_DIRECTION_STRENGTH_K3,SCAN_DIRECTION_STRENGTH_K4,SCAN_DIRECTION_MEAN_K1,SCAN_DIRECTION_MEAN_K2,SCAN_DIRECTION_MEAN_K3,SCAN_DIRECTION_MEAN_K4,PHOT_G_N_OBS,PHOT_G_MEAN_FLUX,PHOT_G_MEAN_FLUX_ERROR,PHOT_G_MEAN_MAG,PHOT_VARIABLE_FLAG,L,B,ECL_LON,ECL_LAT,RAJ2000,DEJ2000,RV,ERV,CAT};

int main(int argc,char* argv[])
{
  /*
    Example: ./integrate.exe 10

    Where: 
    10: Number of particles

    Input: 
    * cloud.data: test particles compatible with A/2017U1 orbit.
    * candidates.csv: candidates to close encounters.

    Output:
    *
  */

  ////////////////////////////////////////////////////
  //INPUT PARAMETERS
  ////////////////////////////////////////////////////
  int Npart=1;
  if(argc>1){
    Npart=atoi(argv[1]);
  }
  VPRINT(stdout,"Analysing %d test particles against candidate stars\n",Npart);
  VPRINT(stdout,"\n");

  ////////////////////////////////////////////////////
  //INITIALIZE CSPICE
  ////////////////////////////////////////////////////
  initSpice();
  
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
  //UNITS
  ////////////////////////////////////////////////////
  UM=MSUN;
  UL=PARSEC;
  UT=YEAR;
  UV=UL/UT;
  G=GCONST/(UL*UL*UL/(UM*UT*UT));
  
  VPRINT(stdout,"Units:\n\tUM = %.5e kg=%.5e Msun\n\tUL = %.17e m = %.17e pc\n\tUT = %.5e s = %.5e yr\n\tUV = %.5e m/s = %.5e km/s\n\tG = %.5e\n",
	 UM,UM/MSUN,UL,UL/PARSEC,UT,UT/YEAR,UV,UV/1e3,G);
  VPRINT(stdout,"\n");

  ////////////////////////////////////////////////////
  //GLOBAL ALLOCATION
  ////////////////////////////////////////////////////
  nsysp=6*Npart;
  nsys=6;

  Ntimesp=10000;
  Ntimes=100;

  x=(double*)malloc(6*sizeof(double));//LSR STATE VECTOR
  xg=(double*)malloc(6*sizeof(double));//GC STATE VECTOR
  dx=(double*)malloc(6*sizeof(double));//LSR STATE VECTOR
  xpmin=(double*)malloc(6*sizeof(double));//GC STATE VECTOR

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
  //GLOBAL PROPERTIES FOR INTEGRATION
  ////////////////////////////////////////////////////
  ip=1;
  params[ip++]=G*MDISK*MSUN/UM;
  params[ip++]=ADISK*PARSEC/UL;
  params[ip++]=BDISK*PARSEC/UL;
  params[ip++]=G*MBULGE*MSUN/UM;
  params[ip++]=0.0;
  params[ip++]=BBULGE*PARSEC/UL;
  params[ip++]=G*MHALO*MSUN/UM;
  params[ip++]=0.0;
  params[ip++]=BHALO*PARSEC/UL;

  ////////////////////////////////////////////////////
  //READ PARTICLES POSITION
  ////////////////////////////////////////////////////
  FILE *fc;
  if((fc=fopen("cloud-int.csv","r"))==NULL){
    fprintf(stdout,"Houston we've got a problem\n");
    exit(0);
  }
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

  n=0;
  fc=fopen("candidates.csv","r");
  FILE *fp=fopen("potential.csv","w");
  fgets(line,MAXLINE,fc);//HEADER
  fprintf(fp,"dyntmin,dyndmin,dynvrel,%s",line);
  char values[MAXLINE];
  while(fgets(line,MAXLINE,fc)!=NULL){
    
    //PARSE FIELDS
    strcpy(values,line);
    parseLine(line,fields,&nfields);
    n++;

    //if(strcmp(fields[HIP],"40170.0")!=0) continue;
    //if(strcmp(fields[TYCHO2_ID],"6436-395-1")!=0) continue;
    //if(strcmp(fields[HIP],"95319.0")!=0) continue;
    fprintf(stdout,"Star %d,%s,%s:\n",n,fields[HIP],fields[TYCHO2_ID]);

    //LMA ENCOUNTER TIME
    duration=atof(fields[TMIN])*YEAR/UT;
    fprintf(stdout,"\tLMA time:%e (compared to %e)\n",duration,tsp[Ntimesp-1]);
    dmin=atof(fields[DMIN]);
    fprintf(stdout,"\tLMA distance:%e\n",dmin);
    hstep=fabs(duration)/(10*Ntimes);
    if(1.1*duration<tsp[Ntimesp-1]){
      fprintf(stdout,"\t\tEncounters happens out of the integration range (%e)...\n",tsp[Ntimesp-1]);
      continue;
    }

    //SPATIAL COORDINATES AND LSR VELOCITY
    for(int k=0;k<6;k++) x[k]=atof(fields[POSTARX+k]);
    vscl_c(PARSEC/1e3,x,x);

    //CONVERT TO GC
    LSR2GC(x,xg);
    vscl_c(1e3/UL,xg,xg);//SET UNITS
    vscl_c(1e3/UV,xg+3,xg+3);
    cart2polar(xg,xInt0,1.0);//CONVERSION TO CYLINDRICAL
    fprintf(stdout,"\tInitial position: %s\n",vec2strn(xInt0,6,"%.5e "));
    //exit(0);

    //INTEGRATE 
    duration+=duration/10;
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
      VPRINT(stdout,"\t\tNominal particle position 1:%s\n",vec2strn(xp1,6,"%.5e "));
      VPRINT(stdout,"\t\tNominal particle position 2:%s\n",vec2strn(xp2,6,"%.5e "));

      //DIFFERENCE
      vsubg_c(x,xp1,6,dx);
      VPRINT(stdout,"\t\tDifference 1: [%s] (pos:%.5e,vel:%.5e]\n",
	      vec2strn(dx,6,"%.5e,"),vnorm_c(dx),vnorm_c(dx+3)*UV);
      vsubg_c(x,xp2,6,dx);
      VPRINT(stdout,"\t\tDifference 2: [%s] (pos:%.5e,vel:%.5e]\n",
	      vec2strn(dx,6,"%.5e,"),vnorm_c(dx),vnorm_c(dx+3)*UV);

      //CALCULATE DISTANCE FROM PARTICLE TO INTERVAL
      dmin=distancePointLine(x,xp1,xp2,&ftmin);
      VPRINT(stdout,"\t\tDistance: %.5e\n",dmin);
      
      //CORRECTED TIME
      tmin=tsp[it]+ftmin*(tsp[it+1]-tsp[it]);
      VPRINT(stdout,"\t\tCorrected time:%e\n",tmin);

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
    if(dyn_dmin<2.0){
      fprintf(fp,"%.5lf,%.5lf,%.5lf,%s",dyn_tmin,dyn_dmin,dyn_vrel,values);
    }
    
    //break;
  }
  fclose(fc);
  fclose(fp);
    
  return 0;
}
