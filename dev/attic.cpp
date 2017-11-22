      if(j==0){
	//CONVERT POSITION TO J2000 AND SAVE POSITION
	vscl_c(UL/1E3,X,Xu);vscl_c(UV/1E3,X+3,Xu+3);

	mxv_c(M_Eclip_J2000,Xu,posJ2000);

	mxv_c(M_Eclip_J2000,Xu+3,posJ2000+3);

	recrad_c(posJ2000,&d,&RA,&DEC);

	mxv_c(M_Eclip_Galactic,Xu,posGalactic);

	mxv_c(M_Eclip_Galactic,Xu+3,posGalactic+3);

	recrad_c(posGalactic,&d,&l,&b);
	
	/*
	fprintf(stdout,"%s %e\n",vec2strn(posGalactic,3,"%e "),d);
	exit(0);
	*/

	fprintf(ftraj,"%e %e %e %e %e %e\n",
		(t*UT-tini)/YEAR,d,l*RAD,b*RAD,RA*RAD,DEC*RAD);
      }

    if(j==0) fclose(ftraj);

    FILE* ftraj;
    if(j==0){
      npoints=100;
      ftraj=fopen("nominal-trajectory.dat","w");
    }else{
      npoints=2;
    }

    /*
    t_step=duration/(npoints-1);
    h_used=h;
    for(i=0;i<npoints;i++) {
      deltat=(t-tini/UT)*UT/YEAR;
      if(direction*((t_start+t_step)-tend)>0) t_step=(tend-t_start);
      t_stop = t_start + t_step;
      h_used = h;
      do {
	while(1){
	  status=Gragg_Bulirsch_Stoer(EoM,X0,X,t,h_used,&h_next,1.0,
				      TOLERANCE,EXTMET,params);
	  if(status) h_used/=4.0;
	  else break;
	}
	t+=h_used;
	copyVec(X0,X,6);
	if(direction*(t+h_next-t_stop)>0) h_used=t_stop-t;
	else h_used=h_next;
      }while(direction*(t-(t_stop-direction*1.e-10))<0);

      if(direction*(t-t_stop)>0){
	h_adjust=(t_stop-t);
	status=Gragg_Bulirsch_Stoer(EoM,X0,X,t,h_adjust,&h_next,1.0,
				    TOLERANCE,EXTMET,params);
	copyVec(X0,X,6);
	t=t_stop;
      }

      t_start = t;
      if(direction*(t_start-tend)>0) break;
    }
    fprintf(stdout,"%e %s\n",t,vec2strn(X,6,"%e "));
    exit(0);
    //*/
    //*


  ////////////////////////////////////////////////////
  //DETERMINE TIME OF ASYMPTOTIC ELEMENTS
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
  integrateEoM(tini/UT,X0,h,npoints,duration,6,EoM,params,ts,Xout);
  vscl_c(UL/1E3,Xout[1],Xu);vscl_c(UV/1E3,Xout[1]+3,Xu+3);
  oscelt_c(Xu,ts[1],MUTOT,elts);
  conics_c(elements,tini+,position);

    
    //DETERMINE PREDICTED POSITION USING LMAAA
    Xend[0]=Xu[0]+(tfut-ts[1]*UT)*Xu[3];
    Xend[1]=Xu[1]+(tfut-ts[1]*UT)*Xu[4];
    Xend[2]=Xu[2]+(tfut-ts[1]*UT)*Xu[5];
    vsub_c(Xend,Xref,Xdif);
    dp=vnorm_c(Xdif)*1e4/UL;
    VPRINT(stdout,"Dur.LMA.=%e\n",dp);

    if(dp<(ds/10) && !qlma){
      VPRINT(stdout,"Dur.LMA.Found!\n");
      durlma=tdur;
      dlma=vnorm_c(Xu)*1e3/AU;
      dplma=dp;
      qlma=1;
      break;
    }

    /*
    //EPOCH
    tstar=atof(fields[REF_EPOCH]);
    VPRINT(stdout,"\tEpoch (year): = %lf\n",tstar);
    sprintf(ctmp,"01/01/%d 00:00:00.000",(int)tstar);
    str2et_c(ctmp,&tstar);
    deltet_c(tstar,"et",&dt);
    tstar-=dt;
    VPRINT(stdout,"\tEpoch: tdb = %lf\n",tstar);
    */


    /*
    //HD39587,HIP27913
    mura=-184.0;dmura=1.4;
    mudec=-87.0;dmura=1.1;
    vr=-13.3;dvr=0.3;
    par=104.1;dpar=5.8;
    /*
    //HR 4867, HIP62512
    mura=+107.0;dmura=3.5;
    mudec=-5.0;dmudec=2.5;
    vr=-12.0;dvr=3.6;
    par=41.4;dpar=12.1;
    //*/
    /*
    //HD115043, HIP64532
    mura=+110.0;dmura=3.8;
    mudec=-35.0;dmudec=3.3;
    vr=-8.8;dvr=2.0;
    par=49.5;dpar=15.7;
    //*/
    /*
    ra=(1+49./60+23.35579/3600)*15;
    dec=-(10+42./60+12.8593/3600.);
    mura=-144;dmura=3.4;
    mudec=-88;dmudec=3.2;
    par=46.4;dpar=6.7;
    vr=-1.5;dvr=1.6;
    //*/
    calcUVW(ra,dec,par,dpar,mura,dmura,mudec,dmudec,vr,dvr,UVW,dUVW);


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//INTEGRATION
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////
  //INTEGRATING TEST PARTICLES FOR A LONG PERIOD 
  ////////////////////////////////////////////////////
  //TIME CONDITIONS
  duration=-1e8*YEAR/UT;

  //INTEGRATE
  params[0]=nsysp;
  integrateEoM(0,xIntp0,hstep,Ntimesp,duration,
	       nsysp,EoMGalactic,params,
	       tsp,xIntp);

  //SAVE POSITIONS IN FILE
  fc=fopen("cloud-int.csv","w");

  //HEADER
  fprintf(fc,"t,");
  for(int j=0;j<Npart;j++){
    fprintf(fc,"part%d-R,part%d-phi,part%d-Z,part%d-vR,part%d-dphi,part%d-vZ,",j,j,j,j,j,j);
    fprintf(fc,"part%d-x,part%d-y,part%d-z,part%d-vx,part%d-xy,part%d-vz,",j,j,j,j,j,j);
  }
  fprintf(fc,"dummy\n");

  //SAVE PARTICLE POSITIONS IN POLAR AND CARTESIAN
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


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//PROGENITORS
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////
  //GLOBAL ALLOCATION
  ////////////////////////////////////////////////////
  nsysp=6*Npart;
  nsys=6;

  Ntimesp=10000;
  Ntimes=100;
  Nobs=10;

  x=(double*)malloc(6*sizeof(double));//LSR STATE VECTOR
  xg=(double*)malloc(6*sizeof(double));//GC STATE VECTOR
  dx=(double*)malloc(6*sizeof(double));//LSR STATE VECTOR
  xpmin=(double*)malloc(6*sizeof(double));//GC STATE VECTOR

  xIntp0=(double*)malloc(nsysp*sizeof(double));
  xIntc0=(double*)malloc(nsysp*sizeof(double));
  xInt0=(double*)malloc(nsys*sizeof(double));
  xnom0=(double*)malloc(nsysp*sizeof(double));
  xnoms0=(double*)malloc(nsysp*sizeof(double));
  
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

  mobs=(double*)malloc(6*sizeof(double));

  obs=(double**)malloc(Nobs*sizeof(double*));
  for(int i=0;i<Nobs;i++) obs[i]=(double*)malloc(6*sizeof(double));

  cov=(double**)malloc(6*sizeof(double*));
  for(int i=0;i<6;i++) cov[i]=(double*)malloc(6*sizeof(double));

  pxform_c("J2000","GALACTIC",0,M_J2000_Galactic);

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
      x=xIntp0+ip;
      for(int k=0;k<6;k++) x[k]=atof(fields[n++]);
      x=xIntc0+ip;
      for(int k=0;k<6;k++) x[k]=atof(fields[n++]);
    }
    //READ ONLY INITIAL CONDITIONS
    break;
  }
  fclose(fc);
  copyVec(xnoms0,xIntp0,6);

  ////////////////////////////////////////////////////
  //READING POTENTIAL OBJECTS
  ////////////////////////////////////////////////////
  params[0]=nsys;
  int Nsur=10;

  //CHOOSE H FOR PROBABILITY CALCULATIONS
  double hprob=1.0*PARSEC/UL;//PC
  double sigma=wNormalization(hprob);

  n=0;
  fc=fopen("potential.csv","r");
  fgets(line,MAXLINE,fc);//HEADER

  while(fgets(line,MAXLINE,fc)!=NULL){
    
    //PARSE FIELDS
    parseLine(line,fields,&nfields);
    n++;
    //if(strcmp(fields[TYCHO2_ID],"7774-308-1")!=0) continue;

    fprintf(stdout,"Star %d,%s,%s:\n",n,fields[HIP],fields[TYCHO2_ID]);

    //ESTIMATED TIME OF ENCOUNTER
    tmin=atof(fields[DYNTMIN]);
    dmin=atof(fields[DYNDMIN]);
    double tmin0=tmin;
    double tmins=tmin0;

    mint=1.2*tmin;
    maxt=0.8*tmin;

    fprintf(stdout,"\tEncounter estimated time, tmin = %.6e\n",tmin);
    fprintf(stdout,"\t\tTesting range = [%.6e,%.6e]\n",mint,maxt);
    fprintf(stdout,"\tEncounter estimated distance, dmin = %.6e\n",dmin);

    if(!(fabs(tmin0)>0)){
      fprintf(stdout,"This star is too close\n");
      getchar();
      continue;
    }
    
    //INFORMATION REQUIRED
    mobs[0]=ra=atof(fields[RA]);
    dra=atof(fields[RA_ERROR])*MAS;

    mobs[1]=dec=atof(fields[DEC]);
    ddec=atof(fields[DEC_ERROR])*MAS;

    mobs[2]=par=atof(fields[PARALLAX]);
    dpar=atof(fields[PARALLAX_ERROR]);

    mobs[3]=mura=atof(fields[PMRA]);
    dmura=atof(fields[PMRA_ERROR]);

    mobs[4]=mudec=atof(fields[PMDEC]);
    dmudec=atof(fields[PMDEC_ERROR]);

    mobs[5]=vr=atof(fields[RV]);
    dvr=atof(fields[ERV]);

    //COVARIANCE MATRIX
    /*RA*/cov[0][0]=dra*dra;
    cov[0][1]=atof(fields[RA_DEC_CORR])*dra*ddec;
    cov[0][2]=atof(fields[RA_PARALLAX_CORR])*dra*dpar;
    cov[0][3]=atof(fields[RA_PMRA_CORR])*dra*dmura;
    cov[0][4]=atof(fields[RA_PMDEC_CORR])*dra*dmudec;
    cov[0][5]=0.0;
    /*DEC*/cov[1][1]=ddec*ddec;
    cov[1][0]=cov[0][1];
    cov[1][2]=atof(fields[DEC_PARALLAX_CORR])*ddec*dpar;
    cov[1][3]=atof(fields[DEC_PMRA_CORR])*ddec*dmura;
    cov[1][4]=atof(fields[DEC_PMDEC_CORR])*ddec*dmudec;
    cov[1][5]=0.0;
    /*PAR*/cov[2][2]=dpar*dpar;
    cov[2][0]=cov[0][2];
    cov[2][1]=cov[1][2];
    cov[2][3]=atof(fields[PARALLAX_PMRA_CORR])*dpar*dmura;
    cov[2][4]=atof(fields[PARALLAX_PMDEC_CORR])*dpar*dmudec;
    cov[2][5]=0.0;
    /*MURA*/cov[3][3]=dmura*dmura;
    cov[3][0]=cov[0][3];
    cov[3][1]=cov[1][3];
    cov[3][2]=cov[2][3];
    cov[3][4]=atof(fields[PMRA_PMDEC_CORR])*dmura*dmudec;
    cov[3][5]=0.0;
    /*MUDEC*/cov[4][4]=dmudec*dmudec;
    cov[4][0]=cov[0][4];
    cov[4][1]=cov[1][4];
    cov[4][2]=cov[2][4];
    cov[4][3]=cov[3][4];
    cov[4][5]=0.0;
    /*RV*/cov[5][5]=dvr*dvr;
    cov[5][0]=cov[0][5];
    cov[5][1]=cov[1][5];
    cov[5][2]=cov[2][5];
    cov[5][3]=cov[3][5];
    cov[5][4]=cov[4][5];
    
    VPRINT(stdout,"\tStellar properties: %s\n",vec2strn(mobs,6,"%.5e "));
    VPRINT(stdout,"\t\tErrors:");
    for(int i=0;i<6;i++) VPRINT(stdout,"%.6e ",sqrt(cov[i][i]));
    VPRINT(stdout,"\n");
    VPRINT(stdout,"\tGalactic coordinates: l = %lf, b = %lf\n",
	   atof(fields[L]),atof(fields[B]));

    VPRINT(stdout,"\tStar Covariance Matrix:\n");
    for(int i=0;i<6;i++)
      fprintf(stdout,"\t\t|%s|\n",vec2strn(cov[i],6,"%-+15.3e"));

    generateMultivariate(cov,mobs,obs,6,Nobs);

    VPRINT(stdout,"\tSurrogate random properties:\n");
    for(int i=Nobs;i-->0;){
      VPRINT(stdout,"\t\tObservation %d: %s\n",i,vec2strn(obs[i],6,"%.10e "));
    }

    //CALCULATE PROBABILITIES
    Pprob=0;
    for(int i=0;i<Nsur;i++){

      VPRINT(stdout,"\tSurrogate %d:\n",i);
      copyVec(xnom0,xnoms0,6);

      //GET KEY PROPERTIES OF SURROGATE
      ra=obs[i][0];
      dec=obs[i][1];
      par=obs[i][2];
      mura=obs[i][3];
      mudec=obs[i][4];
      vr=obs[i][5];

      //INITIAL POSITION RELATIVE TO SUN
      d=AU/tan(par/(60*60*1000.0)*DEG)/PARSEC;
      radrec_c(d,ra*DEG,dec*DEG,xg);
      mxv_c(M_J2000_Galactic,xg,x);
      recrad_c(x,&tmp,&l,&b);
      calcUVW(ra,dec,par,dpar,mura,dmura,mudec,dmudec,vr,dvr,x+3,dx+3);
      //INITIAL POSITION RELATIVE TO GALACTIC CENTER
      vscl_c(PARSEC/1e3,x,x);
      LSR2GC(x,xg);
      vscl_c(1e3/UL,xg,xg);//SET UNITS
      vscl_c(1e3/UV,xg+3,xg+3);
      //INITIAL POLAR COORDINATES OF SURROGATE
      cart2polar(xg,xInt0,1.0);

      VPRINT(stdout,"\t\ttmin0: %.6e\n",tmin0);
      VPRINT(stdout,"\t\tObservations: %s\n",vec2strn(obs[i],6,"%.5e "));
      VPRINT(stdout,"\t\tDistance: %e pc\n",d);
      VPRINT(stdout,"\t\tGalactic coordinates: l = %lf, b = %lf\n",l*RAD,b*RAD);
      VPRINT(stdout,"\t\tInitial position cartesian: %s\n",vec2strn(x,6,"%.5e "));
      VPRINT(stdout,"\t\tInitial position cylindrical: %s\n",vec2strn(xnom0,6,"%.5e "));

      //CALCULATE MINIMUM DISTANCE AND TIME OF NOMINAL SOLUTION TO SURROGATE
      params[0]=6;

      try{
	minDistance2(xInt0,xnom0,tmin0,&dmin,&tmin,params);
      }catch(int e){
	fprintf(stdout,"¡No minimum!\n");
	getchar();
	continue;
      }
      fprintf(stdout,"\t\tMinimum distance at: t = %.6e, d = %.6e\n",tmin,dmin);
      VPRINT(stdout,"\t\tFinal position star: %s\n",vec2strn(xInt0,6,"%.5e "));
      VPRINT(stdout,"\t\tFinal position nominal: %s\n",vec2strn(xnom0,6,"%.5e "));

      //PROPAGATE ALL TEST PARTICLE
      params[0]=6*Npart;
      h=fabs(tmin)/100;
      VPRINT(stdout,"\t\tInitial conditions for all particles: %s\n",vec2strn(xIntp0,6*Npart,"%.5e "));
      integrateEoM(0,xIntp0,h,2,tmin,6*Npart,EoMGalactic,params,ts,xIntp);
      VPRINT(stdout,"\t\tIntegration result for all particles: %s\n",vec2strn(xIntp[1],6*Npart,"%.5e "));

      double D,Dmax=0,*xt1,*xt2,vrel;
      
      //COMPUTE SPH-LIKE PROBABILITY
      double pd,pv;
      Psur=0.0;
      fvel=0.0;
      fprintf(stdout,"\t\tComparing test particle position with star position\n");
      for(int j=0;j<Npart;j++){
	vsubg_c(xInt0,xIntp[1]+6*j,6,dx);
	D=vnorm_c(dx);
	vrel=vnorm_c(dx+3);
	fprintf(stdout,"\t\t\tDistance to test particle %d: v=%.6e,vrel=%.6e\n",j,D,vrel*UV/1e3);
	//CONTRIBUTION TO P FROM DISTANCE
	pd=sigma*wFunction(D,&hprob);

	//CONTRIBUTION TO P FROM VELOCITY
	pv=1.0;
	
	fprintf(stdout,"\t\t\t\tDistance probability: %.6e\n",pd);
	fprintf(stdout,"\t\t\t\tVelocity probability: %.6e\n",pv);

	Psur+=pd*pv;
	//COMPUTE CORRECTION FOR RELATIVE STELLAR VELOCITY
      }

      //COMPUTE CORRECTION FOR STELLAR DISTANCE
      fdist=1/(d*d);
      Psur*=fdist;
      
      //SURROGATE PROBABILITY
      fprintf(stdout,"\t\tSurrogate probability: %.6e\n",Psur);

      //ACCUMULATE
      Pprob+=Psur;
      //getchar();
    }
    fprintf(stdout,"Probability for star: %.6e\n",Pprob);
    getchar();
    //break;
  }
  fclose(fc);


