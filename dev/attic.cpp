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

