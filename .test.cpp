#include <iwander.cpp>
using namespace std;

#define VERBOSE 1

int main(int argc,char* argv[])
{
  #include <iwander.conf>
  #include <wanderer.conf>
  initWander();

  //RUN TESTS
  if(argc>1){
    
    ////////////////////////////////////////////////////
    //TEST ROUTINES
    ////////////////////////////////////////////////////
    int Ndir=200;
    double **rvec=matrixAllocate(Ndir,3);
    double vd,qd,fd;
    double qmin=70.0*DEG,qmax=80.0*DEG;
    double fmin=0.0*DEG,fmax=1.0*DEG;
    int idir=0;
    while(idir<Ndir){
      gsl_ran_dir_3d(RAND,&rvec[idir][0],&rvec[idir][1],&rvec[idir][2]);
      recsph_c(rvec[idir],&vd,&qd,&fd);
      if((qd<qmin)||(qd>qmax)||(fd<fmin)||(fd>fmax)) continue;
      idir++;
    }
    fprintf(stdout,"Done\n");
    double dO=solidAngle(rvec,Ndir);
    fprintf(stdout,"Solid angle:%e\n",dO);
    fprintf(stdout,"Solid angle real:%e\n",(fmax-fmin)*(cos(qmin)-cos(qmax)));
    exit(0);
  }
  

  return 0;
}
