    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //INTEGRATE ALL PARTICLES UNTIL TMIN
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    params[0]=6*Ntest;
    hstep=fabs(nomtmin)/100;
    try{
      integrateEoM(0,xIntp0,hstep,2,tmin,6*Ntest,EoMGalactic,params,ts,xIntp);
    }catch(int e){
      fprintf(stdout,"\t\t****No integration for all particles and star %d***\n",n,i);
      Nstar_nointall++;
      continue;
    }
    if(qsingle){
      fprintf(fso,"%s",vec2strn(xInt0,6,"%e "));
      fprintf(fso,"%s\n",vec2strn(xIntp[1],6*Ntest,"%e "));
    }
    fprintf(stdout,"Test particles: %s\n",vec2strn(xIntp[1],6*Ntest,"%e "));

