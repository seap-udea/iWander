//////////////////////////////////////////
//HEADERS
//////////////////////////////////////////
#include <stdexcept>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <SpiceUsr.h>
#include <eph_manager.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_statistics.h>

//////////////////////////////////////////
//MACROS
//////////////////////////////////////////

//VERBOSE OUTPUT
#define VPRINT if(VERBOSE) fprintf
#define print0 fprintf
#define print1 if(VERBOSE>=1) fprintf
#define print2 if(VERBOSE>=2) fprintf
#define print3 if(VERBOSE>=3) fprintf

//COMMON MACROS
#define D2R(x) (x*M_PI/180)
#define R2D(x) (x*180/M_PI)
#define POWI(x,n) gsl_pow_int(x,n)
#define SGN(x) (x<0?-1:+1)
#define MAX(x,y) (x>y?x:y)
#define MIN(x,y) (x<y?x:y)

//NUMERICAL CONSTANTS
#define DEG ((M_PI/180))
#define RAD ((180/M_PI))
#define PI M_PI

//////////////////////////////////////////
//CSPICE CONSTANTS
//////////////////////////////////////////
#define EARTH_ID "EARTH"
#define ATTEMPTS 12 /*SEE NUMBER_OF_STEPS*/
#define SSB "SOLAR SYSTEM BARYCENTER"

//FOR ABSOLUTE EPHEMERIS
#define ECJ2000 "ECLIPJ2000"
//FOR LOCAL EPHEMERIS
#define J2000 "J2000"

/*
  Range of Ephemeris Times where data to calculate precise ITRF93
  frame directions is availanle
 */

// 01/01/2000 00:00:00.000 UTC
//#define ETINI -4.313581609e+04

// 01/21/1962 00:00:00.000 UTC
#define ETINI -1.197460759e+09
// 07/17/2037 00:00:00.000 UTC
#define ETEND 1.184673668e+09

//////////////////////////////////////////
//CONSTANTS
//////////////////////////////////////////
#define GCONST GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT
#define GKMS (GCONST*1E-9)
#define MSUN 1.9885E30/*kg*/
#define YEAR (365.25*GSL_CONST_MKSA_DAY)
#define DAY GSL_CONST_MKSA_DAY
#define AU GSL_CONST_MKSA_ASTRONOMICAL_UNIT
#define PARSEC GSL_CONST_MKSA_PARSEC
#define VESC_EARTH 11.217 //km/s
#define MAS 1.5915494309189534E-5 //mas -> rad

//////////////////////////////////////////
//BEHAVIOR
//////////////////////////////////////////
#define TOLERANCE 1E-10
#define EXTMET 1
#define VERBOSE 0
#define HTOL 1E-6
#define MAXSTALL 100
#define MAXCOLS 10000
#define MAXTEXT 200
#define MAXLINE 100000

//////////////////////////////////////////
//OBJECTS
//////////////////////////////////////////
#define NUMOBJS 10

double XFOO[]={99.99,99.99,99.99,99.99,99.99,99.99};
double XNULL[]={0.0,0.0,0.0,0.0,0.0,0.0};

//THESE ARE THE LABELS FOR KERNEL de430
static char* LABELS[]={
  "SUN",
  "MERCURY",
  "VENUS",
  "EARTH",
  "MOON",
  "MARS",
  "JUPITER",
  "SATURN",
  "URANUS",
  "NEPTUNE"
};


static char* OBJS[]={
  "10",/*SUN*/
  "1",/*MERCURY*/
  "2",/*VENUS*/
  "399",/*EARTH*/
  "301",/*MOON*/
  "4",/*MARS*/
  "5",/*JUPITER*/
  "6",/*SATURN*/
  "7",/*URANUS*/
  "8"/*NEPTUNE*/
};

//SEE WIKIPEDIA
static double MASSES[]={
  1.9891E30/*SUN*/,
  3.3022E23/*MERCURY*/,
  4.8685E24/*VENUS*/,
  5.9736E24/*EARTH*/,
  7.349E22/*MOON*/,
  6.4185E23/*MARS*/,
  1.8986E27/*JUPITER*/,
  5.6846E26/*SATURN*/,
  8.6810E25/*URANUS*/,
  1.0243E26/*NEPTUNE*/
};

/*
  Source: documentation DE421
  http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421_announcement.pdf
 */
//DE421
static double GMASSES[]={
  132712440040.944000/*SUN*/,
  22032.090000/*MERCURY*/,
  324858.592000/*VENUS*/,
  398600.436233/*EARTH*/,
  4902.800076/*MOON*/,
  42828.375214/*MARS*/,
  126712764.800000/*JUPITER*/,
  37940585.200000/*SATURN*/,
  5794548.600000/*URANUS*/,
  6836535.000000/*NEPTUNE*/
};
double MUTOT;

//DE431
// static double GMASSES[]={
//   132712440040.944000/*SUN*/,
//   22032.090000/*MERCURY*/,
//   324858.592000/*VENUS*/,
//   398600.436233/*EARTH*/,
//   4902.800076/*MOON*/,
//   42828.375214/*MARS*/,
//   126712764.800000/*JUPITER*/,
//   37940585.200000/*SATURN*/,
//   5794548.600000/*URANUS*/,
//   6836535.000000/*NEPTUNE*/
// };

/*
  Gotten from kernes DE430
 */
static double RADII[]={
  696000/*SUN*/,
  2439.70/*MERCURY*/,
  6051.80/*VENUS*/,
  6378.14/*EARTH*/,
  1737.40/*MOON*/,
  3396.19/*MARS*/,
  71492.00/*JUPITER*/,
  37940585.200000/*SATURN*/,
  5794548.600000/*URANUS*/,
  6836535.000000/*NEPTUNE*/
};

//ALL PLANETS 
static int ACTIVE[]={
  1,//SUN
  1,//MERCURY
  1,//VENUS
  1,//EARTH
  1,//MOON
  1,//MARS
  1,//JUPITER
  1,//SATURN
  1,//URANUS
  1,//NEPTUNE
};

struct ObserverStruct{

  SpiceDouble t;
  SpiceDouble lat,lon,alt;

  //Conversion matrix from ITRF93 to ECLIPJ2000 at time t
  SpiceDouble MEJ[3][3];

  //Conversion matrix from ECLIPJ2000 to EPOCHEQUINOX at time t
  SpiceDouble MEE[3][3];

  //hm convert from itrf93 to local and hi is the inverse
  SpiceDouble hm[3][3];
  SpiceDouble hi[3][3];

  //Position with respect to ITRF93 (Earth) J2000
  SpiceDouble posearth[6];

  //Position with respect to ITRF93 (Earth) ECLIPJ2000
  SpiceDouble posj2000[6];

  //Position with respect to SSB in ECLIPJ2000
  SpiceDouble posabs[6];

  //Position with respect to SSB in ECLIPEPOCH
  SpiceDouble posepoch[6];

  //Rotaion velocity of a still observer with respect to ITRF93
  SpiceDouble v[3];

  //Direction of velocity in the direction (A,a) with respect to ECLIPJ2000
  SpiceDouble uv[3];

  //Earth position at observer epoch
  SpiceDouble earth[6];
};

//////////////////////////////////////////
//GLOBAL VARIABLES
//////////////////////////////////////////
enum COMPONENTS {CX,CY,CZ,CVX,CVY,CVZ};
static int NUMBER_OF_STEPS[]={2,4,6,8,12,16,24,32,48,64,96,128};
double REARTH;
double RPEARTH;
double FEARTH;
gsl_rng* RAND;
double GGLOBAL;
double UL,UM,UT,UV;
double INI_TIME,LAST_TIME;

char FILENAME[10000];
double TELAPS=0.0;
int NELAPS=0;
char LINE[MAXLINE],SLINE[MAXLINE];
char VALUES[MAXLINE];
int NFIELDS=0;
char **FIELDS;

//////////////////////////////////////////
//ROUTINES
//////////////////////////////////////////
double *vectorAllocate(int n)
{
  double *v;
  v=(double*)malloc(n*sizeof(double));
  return v;
}

char *charVectorAllocate(int n)
{
  char *v;
  v=(char*)malloc(n*sizeof(char));
  return v;
}

double **matrixAllocate(int n,int m)
{
  double **M;
  
  M=(double**)malloc(n*sizeof(double*));
  for(int i=0;i<n;i++) M[i]=(double*)malloc(m*sizeof(double));
  
  return M;
}

char **charMatrixAllocate(int n,int m)
{
  char **M;
  
  M=(char**)malloc(n*sizeof(char*));
  for(int i=0;i<n;i++) M[i]=(char*)malloc(m*sizeof(char));
  
  return M;
}

int freeMatrix(double **mat,int n,int m)
{
  for(int i=n;i-->0;){
    free(mat[i]);
  }
  free(mat);
  return 0;
}

double elapsedTime(int iprev=1,int iflag=1)
{
  static double times[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  struct timeval  tv;
  double t,dt;
  gettimeofday(&tv, NULL);
  t=(double)(tv.tv_usec/1000000.+tv.tv_sec);
  dt=t-times[iprev];
  times[iflag]=t;
  return dt;
}

void errorGSL(const char * reason,const char * file,int line,int gsl_errno){throw(1);}
int initWander(void)
{
  SpiceInt i,n;
  SpiceDouble radii[3];
  
  //INITIALIZE GSL ERROR HANDLER
  gsl_set_error_handler(&errorGSL); 

  //KERNELS
  furnsh_c("db/kernels/kernels.txt");

  //PLANETARY RADII
  for(i=0;i<10;i++){
    bodvrd_c(LABELS[i],"RADII",3,&n,radii);
    RADII[i]=radii[0];
  }

  //TOTAL MU
  MUTOT=0.0;
  for(i=NUMOBJS;i-->0;){
    if(!ACTIVE[i]) continue;
    MUTOT+=GMASSES[i];
  }

  //RANDOM NUMBERS
  RAND=gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(RAND,time(NULL));
  gsl_rng_set(RAND,3);

  //TIME
  elapsedTime(0,0);
  
  //ALLOCATE GLOBAL MATRICES
  FIELDS=charMatrixAllocate(MAXCOLS,MAXTEXT);

  return 0;
}

int initRadii(void)
{

}

char* dec2sex(double dec)
{
  double d,m,s;
  int id,im,sgn;
  char *str=(char*)calloc(sizeof(char),100); 
  d=fabs(dec);
  sgn=dec/d;
  id=floor(d);
  m=(d-id)*60;
  im=floor(m);
  s=(m-im)*60;
  sprintf(str,"%+d:%02d:%.3f",sgn*id,im,s);
  return str;
}

double sex2dec(double d,double m,double s)
{
  double s2d;
  s2d=d+m/60.0+s/3600.0;
  return s2d;
}

char* vec2str(double vec[],char frm[]="%.8e ")
{
  char format[100];
  char *str=(char*)calloc(sizeof(char),100); 
  sprintf(format,"%s%s%s",frm,frm,frm);
  sprintf(str,format,vec[0],vec[1],vec[2]);
  return str;
}

char* vec2strn(double vec[],int n,char frm[]="%.8e ")
{
  int i;
  char *format=charVectorAllocate(100);
  char *str=(char*)calloc(sizeof(char),100*n);
  sprintf(format,"%ss%s","%",frm);
  for(i=0;i<n;i++) sprintf(str,format,str,vec[i]);

  free(format);
  return str;
}

double greatCircleDistance(double lam1,double lam2,
			   double phi1,double phi2)
{
  double d;

  //HARVESINE FORMULA
  double sf,sl;
  sf=sin((phi2-phi1)/2);
  sl=sin((lam2-lam1)/2);
  d=2*asin(sqrt(sf*sf+cos(phi1)*cos(phi2)*sl*sl));

  return d;
}

int bodyEphemerisApparent(ConstSpiceChar *body,
			  SpiceDouble t,
			  SpiceDouble lon,SpiceDouble lat,SpiceDouble alt,
			  SpiceDouble *range,
			  SpiceDouble *ltime,
			  SpiceDouble *raJ2000,
			  SpiceDouble *decJ2000,
			  SpiceDouble *ra,
			  SpiceDouble *dec
			  )
{
  SpiceDouble earthSSBJ2000[6];
  SpiceDouble bodyJ2000[6],bodySSBJ2000[6],ltbody;
  SpiceDouble bodyTOPOJ2000[3],bodyTOPOEpoch[3];
  SpiceDouble Dbody,RAbody,DECbody,RAbodyJ2000,DECbodyJ2000;
  SpiceDouble observerITRF93[3],observerJ2000[3],observerSSBJ2000[3];
  SpiceDouble M_J2000_Epoch[3][3]={{1,0,0},{0,1,0},{0,0,1}};
  SpiceDouble M_ITRF93_J2000[3][3];
  SpiceDouble d,lt,ltmp,ltold,lttol=1E-2;
  int i,ie=0,ncn=10;
  double cspeed=clight_c();

  //ROTATION MATRIX AT THE TIME OF EPHEMERIS
  pxform_c("J2000","EARTHTRUEEPOCH",t,M_J2000_Epoch);
  pxform_c("ITRF93","J2000",t,M_ITRF93_J2000);

  //OBSERVER POSITION J2000 RELATIVE TO EARTH CENTER
  georec_c(D2R(lon),D2R(lat),alt/1000.0,REARTH,FEARTH,observerITRF93);
  mxv_c(M_ITRF93_J2000,observerITRF93,observerJ2000);

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //ASTROMETRIC POSITION
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  i=0;
  lt=0.0;ltold=1.0;
  spkezr_c(EARTH_ID,t,"J2000","NONE",SSB,
	   earthSSBJ2000,&ltmp);
  vadd_c(earthSSBJ2000,observerJ2000,observerSSBJ2000);
  while((fabs(lt-ltold)/lt)>=lttol && i<ncn){
    ltold=lt;
    spkezr_c(body,t-lt,"J2000","NONE",SSB,bodySSBJ2000,&ltmp);
    vsub_c(bodySSBJ2000,observerSSBJ2000,bodyTOPOJ2000);
    d=vnorm_c(bodyTOPOJ2000);
    lt=d/cspeed;
    i++;
  }
  recrad_c(bodyTOPOJ2000,&d,&RAbodyJ2000,&DECbodyJ2000);
  
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //CORRECTED POSITION
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //OBSERVER POSITION J2000
  spkezr_c(EARTH_ID,t,"J2000","NONE","EARTH BARYCENTER",
	   earthSSBJ2000,&lt);
  vadd_c(earthSSBJ2000,observerJ2000,observerSSBJ2000);

  //CORRECTED OBJECT POSITION
  spkezr_c(body,t,"J2000","LT+S","EARTH BARYCENTER",
	   bodySSBJ2000,&lt);
  vsub_c(bodySSBJ2000,observerSSBJ2000,bodyTOPOJ2000);

  //PRECESS POSITION
  mxv_c(M_J2000_Epoch,bodyTOPOJ2000,bodyTOPOEpoch);

  //RA & DEC PRECESSED
  recrad_c(bodyTOPOEpoch,&d,&RAbody,&DECbody);

  //RETURNED
  *range=d;
  *ltime=lt;
  *ra=RAbody*180/M_PI/15;
  *dec=DECbody*180/M_PI;
  *raJ2000=RAbodyJ2000*180/M_PI/15;
  *decJ2000=DECbodyJ2000*180/M_PI;
  return 0;
}

/*
  Calculates the julian date.  If et=0 it gives the Ephemeris Julian
  Date.  If et=1 it gives the Jul ian Date in the International Atomic
  time reference.
 */
SpiceDouble t2jd(SpiceDouble t,int et=0)
{
  SpiceDouble deltat;
  deltet_c(t,"ET",&deltat);
  double tjd=unitim_c(t,"ET","JED");
  return tjd-et*deltat/86400.0;
}

/*
  Matrix to convert from planetocentric celestial coordinates to
  horizontal coordinates and viceversa.

  See discussion at:
  https://naif.jpl.nasa.gov/pipermail/spice_discussion/2010-July/000307.html

  h2m: converts from geocentric to topocentric
  h2i: converts from topocentric to geocentric
 */
void hormat(SpiceDouble lat,SpiceDouble lon,SpiceDouble t,SpiceDouble h2m[3][3],SpiceDouble h2i[3][3])
{
  SpiceDouble geopos[3],normal[3],normalJ2000[3],normalEpoch[3];
  SpiceDouble ux[]={1,0,0},uy[]={0,1,0},uz[]={0,0,1},uzJ2000[3],uzEpoch[3];
  SpiceDouble M_J2000_Epoch[3][3]={{1,0,0},{0,1,0},{0,0,1}};
  SpiceDouble M_ITRF93_J2000[3][3];

  //TRANSFORM MATRICES
  pxform_c("J2000","EARTHTRUEEPOCH",t,M_J2000_Epoch);
  pxform_c("ITRF93","J2000",t,M_ITRF93_J2000);
  
  //NORMAL ITRF93
  georec_c(D2R(lon),D2R(lat),0.0,REARTH,FEARTH,geopos);
  surfnm_c(REARTH,REARTH,RPEARTH,geopos,normal);

  //NORMAL EPOCH
  mxv_c(M_ITRF93_J2000,normal,normalJ2000);
  mxv_c(M_J2000_Epoch,normalJ2000,normalEpoch);

  //Z EPOCH
  mxv_c(M_ITRF93_J2000,uz,uzJ2000);
  mxv_c(M_J2000_Epoch,uzJ2000,uzEpoch);

  //TRANSFORM MATRICES
  ucrss_c(normalEpoch,uzEpoch,uy);
  ucrss_c(uy,normalEpoch,ux);
  h2m[0][0]=ux[0];h2m[0][1]=ux[1];h2m[0][2]=ux[2];
  h2m[1][0]=uy[0];h2m[1][1]=uy[1];h2m[1][2]=uy[2];
  h2m[2][0]=normalEpoch[0];h2m[2][1]=normalEpoch[1];h2m[2][2]=normalEpoch[2];
  invert_c(h2m,h2i);
}

/*
  Matrix to convert from planetocentric coordinates to horizontal
  coordinates and viceversa.

  See discussion at:
  https://naif.jpl.nasa.gov/pipermail/spice_discussion/2010-July/000307.html

  h2m: converts from geocentric to topocentric
  h2i: converts from topocentric to geocentric
 */
void horgeo(SpiceDouble lat,SpiceDouble lon,SpiceDouble h2m[3][3],SpiceDouble h2i[3][3])
{
  SpiceDouble geopos[3],normal[3];
  SpiceDouble ux[]={1,0,0},uy[]={0,1,0},uz[]={0,0,1};

  //NORMAL ITRF93
  georec_c(D2R(lon),D2R(lat),0.0,REARTH,FEARTH,geopos);
  surfnm_c(REARTH,REARTH,RPEARTH,geopos,normal);

  //TRANSFORM MATRICES
  ucrss_c(normal,uz,uy);
  ucrss_c(uy,normal,ux);
  h2m[0][0]=ux[0];h2m[0][1]=ux[1];h2m[0][2]=ux[2];
  h2m[1][0]=uy[0];h2m[1][1]=uy[1];h2m[1][2]=uy[2];
  h2m[2][0]=normal[0];h2m[2][1]=normal[1];h2m[2][2]=normal[2];
  invert_c(h2m,h2i);
}

/*
  Transform from rectangular position to spherical position using the
  astronomical convention of Azimuth.
 */
int rec2hor(double pos[],double *Az,double *h)
{
  double phi,tmp;
  pos[1]*=-1;
  reclat_c(pos,&tmp,Az,h);
  if(*Az<0) *Az+=2*M_PI;
  *h=R2D(*h);
  *Az=R2D(*Az);
  return 0;
}

int copyVec(double tgt[],double src[],int n)
{
  memcpy(tgt,src,n*sizeof(double));
  return 0;
}

int copyVecInt(int tgt[],int src[],int n)
{
  memcpy(tgt,src,n*sizeof(int));
  return 0;
}

int sumVec(double c[],double ca,double a[],double cb,double b[],int n)
{
  int i;
  for(i=n;i-->0;) c[i]=ca*a[i]+cb*b[i];
  return 0;
}

double maxAbsVec(double a[],int n)
{
  int i;
  double max=-1E100;
  for(i=n;i-->0;) if(fabs(a[i])>max) max=fabs(a[i]);
  return max;
}

/*
Adapted from: http://www.mymathlib.com/diffeq/bulirsch_stoer.html
 */

static int Rational_Extrapolation_to_Zero(double *fzero,double tableau[],
					  double x[],double f,int n) 
{
  double t, up, across, denominator, dum;
  int col;

  if (n==0) {  *fzero = f; tableau[0] = f; return 0; }
  if ( x[n] == 0.0 ) { *fzero = f; return -2; }
   
  across = 0.0;                                                        
  up = tableau[0];                                                    
  tableau[0] = f;                                               

  for (col = 1; col <= n; col++) {
    if(tableau[col-1]==0 && across==0){t=0;break;}
    denominator = tableau[col-1] - across;
    if (denominator == 0.0) return -1;
    dum = 1.0 - (tableau[col-1] - up) / denominator;
    denominator = (x[n - col] / x[n]) * dum - 1.0;
    if (denominator == 0.0) return -1;
    t = tableau[col-1] + ( tableau[col-1] - up ) / denominator;
    across = up;
    up = tableau[col];
    tableau[col] = t;
  }
  *fzero = t;
  return 0;
}

static int Polynomial_Extrapolation_to_Zero(double *fzero,double tableau[],
					    double x[], double f, int n )
{
  double back_two_columns;    //  T[row,col-2];
  double old_aux;             //  T[row-1,col];
  double new_value;           //  T[row,col];
  double vertical_diff;       //  T[row,col]-T[row-1,col]
  double backslant_diff;      //  T[row,col]-T[row,col-1]
  double forwardslant_diff;   //  T[row,col]-T[row-1,col-1];
  double denominator;        
  int i;

  if (n == 0) { tableau[0] = f; return 0; }
  if ( x[n] == 0.0 ) { *fzero = f; return -2; }

  back_two_columns = 0.0;
  old_aux = tableau[0];
  tableau[0] = f;
  for (i = 0; i < n; i++) {
    if(tableau[i]==0 && old_aux==0){tableau[n]=0.0;break;}
    vertical_diff = tableau[i] - old_aux;
    backslant_diff = tableau[i] - back_two_columns;
    forwardslant_diff = backslant_diff - vertical_diff;
    denominator = (x[n-i-1]/x[n]) * forwardslant_diff - backslant_diff;
    if (denominator == 0.0) return -1;
    back_two_columns = old_aux;
    old_aux = tableau[i+1];
    tableau[i+1] = tableau[i] + vertical_diff * backslant_diff / denominator;
  }
  *fzero = tableau[n];
  return 0;
}

static int Graggs_Method(int (*f)(double,double*,double*,void*),
			 double y0[],
			 double t0,double t,
			 int NUMBER_OF_STEPS,
			 void *params,
			 double yres[]) {
  
  double* pars=(double*)params;
  int order=(int)pars[0],i;
  double y1[order],dydt[order],y2[order],yaux[order];
  double h=(t-t0)/(double)NUMBER_OF_STEPS;
  double h2=h+h;

  copyVec(yaux,y0,order);
  (*f)(t0,yaux,dydt,params);
  sumVec(y1,1,yaux,h,dydt,order);

  while(--NUMBER_OF_STEPS) {
    t0+=h;
    (*f)(t0,y1,dydt,params);
    sumVec(y2,1,yaux,h2,dydt,order);
    copyVec(yaux,y1,order);
    copyVec(y1,y2,order);
  } 

  (*f)(t,y1,dydt,params);

  sumVec(yres,0.5,yaux,0.5,y1,order);
  sumVec(yres,1,yres,0.5*h,dydt,order);
  return 0;
}

int Gragg_Bulirsch_Stoer(int (*f)(double,double*,double*,void*), 
			 double y0[], double y1[],
			 double t, double h, double *h_new, 
			 double epsilon, double yscale, 
			 int rational_extrapolate,
			 void *params)
{
  double* pars=(double*)params;
  int order=(int)pars[0];
  double step_size2[ATTEMPTS];
  double tableau[order][ATTEMPTS+1];
  double dum;
  double est[order],dest[order],destmax;
  double old_est[order];
  
  int (*Extrapolate)(double*,double*,double*,double,int);
  int i,j;
  int err;

  if(yscale==0.0) return -3;
  if(rational_extrapolate) Extrapolate=Rational_Extrapolation_to_Zero;
  else Extrapolate=Polynomial_Extrapolation_to_Zero;
 
  Graggs_Method(f,y0,t,t+h,NUMBER_OF_STEPS[0],params,est);
  step_size2[0]=(dum=h/(double)NUMBER_OF_STEPS[0],dum*dum);
  
  copyVec(y1,est,order);
  
  for(i=order;i-->0;){
    err=Extrapolate(&y1[i],tableau[i],step_size2,est[i],0);
    if(err<0) return err-1;
  }

  for(i = 1; i < ATTEMPTS; i++) {
    copyVec(old_est,y1,order);
    Graggs_Method(f,y0,t,t+h,NUMBER_OF_STEPS[i],params,est);
    step_size2[i]=(dum=h/(double)NUMBER_OF_STEPS[i],dum*dum);

    for(j=order;j-->0;){
      err=Extrapolate(&y1[j],tableau[j],step_size2,est[j],i);
      if(err<0) return err-1;
    }
    
    sumVec(dest,1.0/yscale,y1,-1.0/yscale,old_est,order);
    destmax=maxAbsVec(dest,order);

    if(destmax<epsilon){
      if(i>1) *h_new=8.0*h/(double)NUMBER_OF_STEPS[i-1];
      else *h_new=h;
      return 0;
    }
  }
  return -1;
}

double energy2B(double X[])
{
  double v=vnorm_c(X+3);
  double r=vnorm_c(X);
  double E=0.5*v*v-1/r;
  return E;
}

int EoM(double t,double y[],double dydt[],void *params) 
{ 
  //COMPUTE THE CONTRIBUTION OF EVERY OBJECT
  int i;
  double r,object[6],R[3],Rmag,tmp,GM,fac;
  double ab;
  int* ps=(int*)params;
  
  fac=UT*UT/(UL/1E3*UL/1E3*UL/1E3);
  dydt[CX]=y[CVX];
  dydt[CY]=y[CVY];
  dydt[CZ]=y[CVZ];
  dydt[CVX]=0.0;
  dydt[CVY]=0.0;
  dydt[CVZ]=0.0;

  for(i=NUMOBJS;i-->0;){
    if(!ACTIVE[i]) continue;
    spkezr_c(OBJS[i],t*UT,ECJ2000,"NONE",SSB,object,&tmp);
    vscl_c(1E3/UL,object,object);
    sumVec(R,1.0,y,-1.0,object,3);
    Rmag=vnorm_c(R);
    if(Rmag*UL/1e3<=RADII[i]){
      fprintf(stderr,"\t\tObject has collided with %s at t = %e days (Rmag = %e, RADII = %e)\n",
	     LABELS[i],t*UT/DAY,Rmag,RADII[i]);
      throw(1);
    }else{
      GM=GMASSES[i]*fac;
      ab=-GM/(Rmag*Rmag*Rmag);
      if(ps[1]==1){
	fprintf(stdout,"\t\tAcceleration body %d: %.17e\n",i,ab);
      }
      sumVec(dydt+3,1.0,dydt+3,ab,R,3);
    }
  }
  return 0;
}

int initObserver(SpiceDouble t,struct ObserverStruct* observer)
{
  SpiceDouble rho,vcirc,vrot[3];
  SpiceDouble lt;
  
  observer->t=t;
  
  //CHECK DATES
  SpiceDouble tref;
  if(t>=ETINI && t<=ETEND){
    tref=t;
  }else{
    SpiceChar UTC[100];
    SpiceDouble dt;
    deltet_c(t,"et",&dt);
    et2utc_c(t+dt,"ISOC",2,100,UTC);
    if(t<ETINI){
      UTC[0]='1';UTC[1]='9';UTC[2]='6';UTC[3]='3';
      str2et_c(UTC,&tref);
      fprintf(stdout,"ETINI = %.10e, t = %.10e, tref = %.10e\n",ETINI,t,tref);
    }
    else{
      UTC[0]='2';UTC[1]='0';UTC[2]='3';UTC[3]='2';
      str2et_c(UTC,&tref);
      fprintf(stdout,"ETEND = %.10e, t = %.10e, tref = %.10e\n",ETEND,t,tref);
    }
  }
  //DEBUGGING
  //printf("t = %e, tref = %e\n",t,tref);

  //CONVERSION FROM EARTH SYSTEM TO ECLIPTIC SYSTEM AT TIME T
  pxform_c("ITRF93",ECJ2000,tref,observer->MEJ);
  /*
  printf("%e,%e,%e\n%e,%e,%e\n%e,%e,%e\n",
	 observer->MEJ[0][1],observer->MEJ[0][2],observer->MEJ[0][3],
	 observer->MEJ[1][1],observer->MEJ[1][2],observer->MEJ[1][3],
	 observer->MEJ[2][1],observer->MEJ[2][2],observer->MEJ[2][3]);
  exit(0);
  */

  //CONVERSION FROM EARTH SYSTEM TO ECLIPTIC SYSTEM AT TIME T
  pxform_c("ECLIPJ2000","EARTHTRUEEPOCH",tref,observer->MEE);

  //LOCATE OBSERVER 
  georec_c(D2R(observer->lon),D2R(observer->lat),observer->alt/1000.0,
	   REARTH,FEARTH,observer->posearth);

  //DEBUGGING
  //printf("Observer = %s\n",vec2str(observer->posearth));

  //TOPOCENTRIC CONVERSION MATRICES
  horgeo(observer->lat,observer->lon,observer->hm,observer->hi);

  //VELOCITY OF OBSERVER DUE TO EARTH ROTATION
  rho=sqrt(observer->posearth[0]*observer->posearth[0]+
	   observer->posearth[1]*observer->posearth[1]);
  vcirc=2*M_PI*rho/GSL_CONST_MKSA_DAY;
  vpack_c(0.0,-vcirc,0.0,vrot);
  mxv_c(observer->hi,vrot,observer->v);

  //POSITION OF THE EARTH
  spkezr_c(EARTH_ID,t,ECJ2000,"NONE","SOLAR SYSTEM BARYCENTER",
	   observer->earth,&lt);

  //DEBUGGING
  //printf("Earth ECJ2000 = %s\n",vec2str(observer->earth,"%.17e "));

  //POSITION WITH RESPECT TO SSB IN ECLIPJ2000
  mxv_c(observer->MEJ,observer->posearth,observer->posj2000);
  vadd_c(observer->earth,observer->posj2000,observer->posabs);

  //DEBUGGING
  //printf("Observer ECJ2000 = %s\n",vec2str(observer->posabs,"%.17e "));

  //POSITION WITH RESPECT TO SSB IN ECLIPEPOCH
  //This is not working
  mxv_c(observer->MEE,observer->posabs,observer->posepoch);
  mxv_c(observer->MEE,(observer->posabs)+3,(observer->posepoch)+3);

}

int observerVelocity(struct ObserverStruct *observer,
		     SpiceDouble elev,SpiceDouble Az,SpiceDouble v)
{
  ////////////////////////////////////////////////////////////// 
  //DETERMINE THE OBSERVER VELOCITY WITH RESPECT TO SSB
  ////////////////////////////////////////////////////////////// 

  //VELOCITY OF OBSERVER IN SPACE W.R.T. TO LOCAL REFERENCE
  SpiceDouble cA=cos(D2R(Az)),sA=sin(D2R(Az)),ch=cos(D2R(elev)),sh=sin(D2R(elev));
  SpiceDouble vloc[3];
  if(Az==0 && elev==0 && v==0){
    vpack_c(0,0,0,vloc);
  }else{
    vpack_c(v*ch*cA,-v*ch*sA,v*sh,vloc);
  }

  //IMPACT VELOCITY IS THE INVERSE
  vscl_c(-1,vloc,vloc);

  //VELOCITY OF OBSERVER IN SPACE W.R.T. TO ITRF93
  SpiceDouble vmot[3];
  mxv_c(observer->hi,vloc,vmot);

  //TOTAL VELOCITY WITH RESPECT ITRF93
  vadd_c(observer->v,vmot,observer->posearth+3);

  //VELOCITY W.R.T. EARTH CENTER IN ECLIPJ2000 RF
  mxv_c(observer->MEJ,observer->posearth+3,observer->posj2000+3);

  //VELOCITY W.R.T. SOLAR SYSTEM BARYCENTER IN J2000 RF
  vadd_c(observer->earth+3,observer->posj2000+3,observer->posabs+3);

  /*NEW*/
  //************************************************************
  //COMPUTING DIRECTION OF INCOMING VELOCITY IN ECLIPJ2000
  /*
    It does not take into account Earth rotation effect on velocity
   */
  //************************************************************
  SpiceDouble uv[3],nuv;
  if(Az==0 && elev==0 && v==0){
    vpack_c(0,0,0,vloc);
  }else{
    vpack_c(ch*cA,-ch*sA,sh,uv);
  }
  //IMPACT DIRECTION IS THE INVERSE
  vscl_c(-1,uv,uv);
  //DIRECTION IN SPACE W.R.T. TO ITRF93
  mxv_c(observer->hi,uv,vmot);
  //DIRECTION IN SPACE W.R.T. ECLIPJ2000
  mxv_c(observer->MEJ,vmot,observer->uv);

  return 0;
}

int rayPropagation(struct ObserverStruct *observer,
		   SpiceDouble deltat,
		   SpiceDouble elements[6])
{
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //INITIAL CONDITIONS FOR PROPAGATION
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SpiceDouble x=observer->posabs[0];
  SpiceDouble y=observer->posabs[1];
  SpiceDouble z=observer->posabs[2];
  SpiceDouble vx=observer->posabs[3];
  SpiceDouble vy=observer->posabs[4];
  SpiceDouble vz=observer->posabs[5];

  deltat*=365.25*GSL_CONST_MKSA_DAY;
  SpiceDouble tini=observer->t;
  double direction=deltat/abs(deltat);
  double params[]={6};

  //UNITS
  UL=GSL_CONST_MKSA_ASTRONOMICAL_UNIT;
  UM=MSUN;
  GGLOBAL=1.0;
  UT=sqrt(UL*UL*UL/(GCONST*UM));
  UV=UL/UT;

  //INITIAL CONDITIONS
  double X0[6],X[6],Xu[6],E[8],a;
  vpack_c(x*1E3/UL,y*1E3/UL,z*1E3/UL,X0);
  vpack_c(vx*1E3/UV,vy*1E3/UV,vz*1E3/UV,X0+3);

  //DYNAMICAL TIMESCALE
  a=vnorm_c(X0);
  double tdyn=2*M_PI*sqrt(a*a*a/(GGLOBAL*MSUN/UM));

  //TIME LIMITS
  deltat/=UT;

  //GUESS TIME-STEP AS 1/1000 OF THE CHARACTERISTIC DYNAMICAL TIME
  double h=direction*tdyn/1000.0,h_used,h_next,h_adjust,delt;

  double t_start=tini/UT;
  double t_step=deltat;
  double tend=t_start+deltat;
  double t_stop=tend;

  double t=t_start;

  //INTEGRATION
  int status;

  t_stop = t_start + t_step;

  h_used = h;
  int nstall=0;
  do {
    //ADJUST H UNTIL OBTAINING A PROPER TIMESTEP
    while(1){
      status=Gragg_Bulirsch_Stoer(EoM,X0,X,t,h_used,&h_next,1.0,TOLERANCE,EXTMET,params);
      if(status) h_used/=4.0;
      else break;
    }
    if(fabs(h_used/t_step)<HTOL) nstall++;
    else nstall=0;
    if(nstall>MAXSTALL){
      fprintf(stderr,"\t\tIntegration has stalled at t = %e days with h/DT = %e\n",
	      t*UT/DAY,h_used/t_step);
      throw(1);
    }
    t+=h_used;
    copyVec(X0,X,6);
    if(direction*(t+h_next-t_stop)>0) h_used=t+h_next-t_stop;
    else h_used=h_next;
  }while(direction*(t-(t_stop-direction*1.e-10))<0);

  //PREVIOUS PROCEDURE WILL LEAVE YOU STILL APART FROM FINAL TIME, SO ADJUST
  if(direction*(t-t_stop)>0){
    h_adjust=(t_stop-t);
    status=Gragg_Bulirsch_Stoer(EoM,X0,X,t,h_adjust,&h_next,1.0,TOLERANCE,EXTMET,params);
    copyVec(X0,X,6);
    t=t_stop;
  }

  //CONVERTING TO CLASSICAL ELEMENTS IN KM AND KM/S
  vscl_c(UL/1E3,X0,Xu);vscl_c(UV/1E3,X0+3,Xu+3);
  oscelt_c(Xu,t*UT,GKMS*MSUN,E);

  vsclg_c(180/M_PI,E+2,4,E+2);
  vsclg_c(1E3/UL,E,1,E);

  //STORE THE ELEMENTS
  copyVec(elements,E,6);

  return 0;
}

int argsError(char* pname,char* msg="Bad options.")
{
    char cmd[1000];
    sprintf(cmd,"cat .help/%s.help",pname);
    fprintf(stderr,msg);
    system(cmd);
    exit(1);
}

char *str_replace(char *orig, char *rep, char *with) 
{
  /*
    Source:
    http://stackoverflow.com/questions/779875/what-is-the-function-to-replace-string-in-c
    You must free the result if result is non-NULL.
  */
  char *result; // the return string
  char *ins;    // the next insert point
  char *tmp;    // varies
  int len_rep;  // length of rep
  int len_with; // length of with
  int len_front; // distance between rep and end of last rep
  int count;    // number of replacements

  if (!orig)
    return NULL;
  if (!rep)
    rep = "";
  len_rep = strlen(rep);
  if (!with)
    with = "";
  len_with = strlen(with);

  ins = orig;
  for (count = 0; tmp = strstr(ins, rep); ++count) {
    ins = tmp + len_rep;
  }

  // first time through the loop, all the variable are set correctly
  // from here on,
  //    tmp points to the end of the result string
  //    ins points to the next occurrence of rep in orig
  //    orig points to the remainder of orig after "end of rep"
  tmp = result = (char*) malloc(strlen(orig) + (len_with - len_rep) * count + 1);

  if (!result)
    return NULL;

  while (count--) {
    ins = strstr(orig, rep);
    len_front = ins - orig;
    tmp = strncpy(tmp, orig, len_front) + len_front;
    tmp = strcpy(tmp, with) + len_with;
    orig += len_front + len_rep; // move to next "end of rep"
  }
  strcpy(tmp, orig);
  return result;
}

//PARSE LINES READ FROM A COMMA SEPARATED VALUES FILE
int parseLine(char line[],char** cols,int *ncols,char sep[]=",")
{
  int i=0;
  char *found;
  while((found=strsep(&line,sep))!=NULL) strcpy(cols[i++],found);
  *ncols=i;
  return 0;
}

//CALCULATE THE EQUATORIAL->GALACTIC MATRIX 
//ASSUMED CONSTANTS
#define RANGP 192.85948 //RA NORTH GALACTIC POLE
#define DECNGP 27.12825 //DEC NORTH GALACTIC POLE
#define THETA0 122.93192 //theta_0 GALACTIC SYSTEM
#define KC1 9.7779e8 //PC KM^-1 YR^-1, C1 CONVERSION UNIT CONSTANT IN BAILER-JONES(2015)
#define KC2 4.74047 //AU KM^-1 YR^-1, C2 CONVERSION UNIT CONSTANT IN BAILER-JONES(2015)

//GALACTIC COORDINATES OF THE SUN AND GALACTIC POTENTIAL
/*
  See
  http://docs.astropy.org/en/stable/api/astropy.coordinates.Galactocentric.html
  for details
 */
//TRUNCATION RADIUS OF THE SOLAR SYSTEM
#define RTRUNC 5E4 //AU
//GIVEN BY https://ui.adsabs.harvard.edu/#abs/2009ApJ…692.1075G/abstract
//REPLACED BY https://academic.oup.com/mnras/article-abstract/465/1/472/2417491?redirectedFrom=PDF
#define ROSUN 8.2E3 //PC
//#define ROSUN 8.0E3 //PC
//GIVEN BY https://ui.adsabs.harvard.edu/#abs/2001ApJ…553..184C/abstract
//REPLACED BY: https://academic.oup.com/mnras/article-abstract/465/1/472/2417491?redirectedFrom=PDF
#define ZSUN 17.0 //PC, 
//#define ZSUN 10.0 //PC, 
#define PHISUN 0.0 //DEG
//GIVEN BY: https://ui.adsabs.harvard.edu/#abs/2010MNRAS.403.1829S/abstract
#define USUN 11.1 //km/s, U IN DIRECTION OF -r 
#define VSUN 12.24 //km/s V IN DIRECTION OF -phi
#define WSUN 7.25 //km/s, W IN DIRECTION OF z
//GIVEN BY: https://ui.adsabs.harvard.edu/#abs/2015ApJS..216...29B/abstract
#define VCIRC 220.0 //km/s

//GIVEN BY BOILER-JONES (2017)
#define MDISK 7.91e10 //Msun
#define ADISK 3500.0 //PC
#define BDISK 250.0 //PC
#define MBULGE 1.40e10 //Msun
#define BBULGE 350.0 //PC
#define MHALO 6.98e11 //Msun
#define BHALO 24000.0 //PC

/*
  Gradient of the Kuzmin Potential

  phi(R,z) = -mu / sqrt( R^2 + (a+ sqrt(z^2+b^2))^2)

  Let's call: 
     u = sqrt(b^2+z^2)
     v = a + u

  dphi/dz = mu z v / [ u (v^2 + R^2)^(3/2)]

  dphi/dR = mu R / ((v^2 + R^2)^(3/2)

 */
int gradKuzmin(double R,double z,double mu,double a,double b,
		double *dphidR,double *dphidz)
{
  double u,v,D;
  
  u=sqrt(b*b+z*z);
  v=a+u;
  D=pow(v*v+R*R,3./2);
  *dphidR=mu*R/D;
  *dphidz=mu*z*v/(u*D);

  VPRINT(stdout,"R=%e, z=%e, mu=%e, a=%e, b=%e, dphidR=%e,dphidz=%e\n",R,z,mu,a,b,*dphidR,*dphidz);

  return 0;
}

/*
  Three components potential of Miyamoto & Nagai (1975)
  
  Compute the derivatives of the potential
*/
int gradGalacticPotential(double R,double q,double z,
			  double *dphidR,double *dphidq,double *dphidz,
			  void *params)
{
  //Parameters: Component 0 is reserved
  double* pars=(double*)params;
  double mud=pars[1],ad=pars[2],bd=pars[3];
  double mub=pars[4],ab=pars[5],bb=pars[6];
  double muh=pars[7],ah=pars[8],bh=pars[9];
  double dfdR=0.0,dfdz=0.0;
  *dphidR=0.0;
  *dphidq=0.0;
  *dphidz=0.0;
    
  //Disk
  //printf("R=%e z=%e mu=%e a=%e b=%e\n",R,z,mud,ad,bd);
  gradKuzmin(R,z,mud,ad,bd,&dfdR,&dfdz);
  *dphidR+=dfdR;
  *dphidz+=dfdz;
  //printf("dPHIdz (disk) = %e\n",dfdz);
  //printf("dPHIdR (disk) = %e\n",dfdR);

  //Bulge
  gradKuzmin(R,z,mub,ab,bb,&dfdR,&dfdz);
  *dphidR+=dfdR;
  *dphidz+=dfdz;
  //printf("dPHIdz (bulge) = %e\n",dfdz);
  //printf("dPHIdR (bulge) = %e\n",dfdR);
  
  //dphi/dz
  gradKuzmin(R,z,muh,ah,bh,&dfdR,&dfdz);
  *dphidR+=dfdR;
  *dphidz+=dfdz;
  //printf("dPHIdz (halo) = %e\n",dfdz);
  //printf("dPHIdR (halo) = %e\n",dfdR);

  return 0;
}


/*
  Following prescription by Johnson & Soderblom, 1987
 */
void TGalacticMatrix(void)
{
  double A1[][3]={{+cos(THETA0*DEG),+sin(THETA0*DEG),0},
		  {+sin(THETA0*DEG),-cos(THETA0*DEG),0},
		  {0,0,+1}};
  double A2[][3]={{-sin(DECNGP*DEG),0,+cos(DECNGP*DEG)},
		  {0,-1,0},
		  {+cos(DECNGP*DEG),0,+sin(DECNGP*DEG)}};
  double A3[][3]={{+cos(RANGP*DEG),+sin(RANGP*DEG),0},
		  {+sin(RANGP*DEG),-cos(RANGP*DEG),0},
		  {0,0,1}};
  double C1[3][3],C[3][3];

  mxm_c(A2,A3,C1);
  mxm_c(A1,C1,C);
  
  fprintf(stdout,"Galactic matrix (manual):\n");
  fprintf(stdout,"\t%s\n",vec2str(C[0],"%.5e "));
  fprintf(stdout,"\t%s\n",vec2str(C[1],"%.5e "));
  fprintf(stdout,"\t%s\n",vec2str(C[2],"%.5e "));

}

/*
  Following prescription by Johnson & Soderblom, 1987

  Test code:
  
  //TEST UVW
  ra=(1+49./60+23.35579/3600)*15;
  dec=-(10+42./60+12.8593/3600.);
  mura=-144;dmura=3.4;
  mudec=-88;dmudec=3.2;
  par=46.4;dpar=6.7;
  vr=-1.5;dvr=1.6;
  calcUVW(ra,dec,par,dpar,mura,dmura,mudec,dmudec,vr,dvr,UVW,dUVW);
  fprintf(stdout,"UVW = %s +/- %s\n",vec2str(UVW,"%.5lf "),vec2str(dUVW,"%.5lf "));
  exit(0);
  
 */
int calcUVW(double ra,double dec,
	    double par,double dpar,
	    double mura,double dmura,
	    double mudec,double dmudec,
	    double vr,double dvr,
	    double UVW[3],double dUVW[3])
{
  //CORRECTION BY DECLINATION
  //MURA SHOULD BE ALREADY CORRECTED BY DECLINATION
  //mura=mura*cos(dec*DEG);

  //COMPUTE VSKY
  double kc2=KC2;
  //kc2=4.74057;

  double vsky[]={
    vr,/*RADIAL, km/s*/
    kc2*mura/par,/*RA, km/s*/
    kc2*mudec/par/*DEC, km/s*/
  };
  
  //TRANSFORM TO LSR

  //ORIGINAL
  /*
  double TM[][3]={{-0.06699,-0.87276,-0.48354},
		  {+0.49273,-0.45035,+0.74458},
		  {-0.86760,-0.18837,+0.46020}};
  //*/
  //*
  double TM[3][3];
  pxform_c("J2000","GALACTIC",0,TM);
  //*/
	  
  double AM[][3]={{+cos(ra*DEG)*cos(dec*DEG),-sin(ra*DEG),-cos(ra*DEG)*sin(dec*DEG)},
		  {+sin(ra*DEG)*cos(dec*DEG),+cos(ra*DEG),-sin(ra*DEG)*sin(dec*DEG)},
		  {+sin(dec*DEG),0,+cos(dec*DEG)}}; 
  double BM[3][3];
  mxm_c(TM,AM,BM);
  mxv_c(BM,vsky,UVW);

  //COMPUTE ERRORS
  double BM2[][3]={{BM[0][0]*BM[0][0],BM[0][1]*BM[0][1],BM[0][2]*BM[0][2]},
		   {BM[1][0]*BM[1][0],BM[1][1]*BM[1][1],BM[1][2]*BM[1][2]},
		   {BM[2][0]*BM[2][0],BM[2][1]*BM[2][1],BM[2][2]*BM[2][2]}};
  double aux1=(KC2/par)*(KC2/par),aux2=(dpar/par)*(dpar/par);
  double dmu[]={
    dvr*dvr,
    aux1*(dmura*dmura+(mura*mura)*aux2),
    aux1*(dmudec*dmudec+(mudec*mudec)*aux2)
  };
  double bcr[]={
    BM[0][1]*BM[0][2],
    BM[1][1]*BM[1][2],
    BM[2][1]*BM[2][2]
  };

  double C[3][3],V[3];
  mxv_c(BM2,dmu,V);
  vscl_c(2*mura*mudec/(par*par)*KC2*aux2,bcr,bcr);
  
  vadd_c(V,bcr,dUVW);
  dUVW[0]=sqrt(dUVW[0]);dUVW[1]=sqrt(dUVW[1]);dUVW[2]=sqrt(dUVW[2]);
  
  return 0;
}

/*
  Transform XYZ,UVW w.r.t. LSR to X,Y,Z,VX,VY,VZ w.r.t. Galactic Center
  r: XYZ (Cartesian Galactic Coordinates), UVW 
  Units: km,km/s

  Test:
  double xLSR_s[]={-3.16028e+02*PARSEC/1e3,
		   1.84156e+01*PARSEC/1e3,
		   -3.58527e+02*PARSEC/1e3, 
                   1.33746e+01, -1.72967e+01, -1.54271e+01};
  VPRINT(stdout,"LSR (star) = %s\n",vec2strn(xLSR_s,6,"%.17e "));
  double xGC_s[6];
  
  //CONVERT FROM LSR TO GC
  LSR2GC(xLSR_s,xGC_s);
  vscl_c(1e3/PARSEC,xGC_s,xGC_s);
  VPRINT(stdout,"GC (star) = %s (pc,km/s)\n",vec2strn(xGC_s,6,"%.17e "));

  Compared to AstroPy:

  import astropy.coordinates as coord
  import astropy.units as u

  c1 = coord.ICRS(ra=45.1128*u.degree, dec=0.380844*u.degree,
    distance=(2.09081*u.mas).to(u.pc, u.parallax()),
    pm_ra_cosdec=-1.57293*np.cos(0.380844*np.pi/180)*u.mas/u.yr,
    pm_dec=-11.6616*u.mas/u.yr,
    radial_velocity=2.061*u.km/u.s)
  
  print(c1.transform_to(coord.Galactocentric)
    (x, y, z) in pc
    (-8617.14887136,  18.41629297, -330.49785565)
    (v_x, v_y, v_z) in km / s
    ( 24.42433714,  214.94327839, -8.22049138)>
 */
int LSR2GC(double LSR[6],double GC[6])
{
  //CONSTANT
  double UVW_Sun[]={USUN,(VSUN+VCIRC),WSUN};
  double theta=asin(ZSUN/ROSUN);
  double MGC[][3]={{cos(theta),0,sin(theta)},
		   {0,1.0,0},
		   {-sin(theta),0,cos(theta)}};

  VPRINT(stdout,"theta = %.3lf degrees\n",theta*RAD);
  VPRINT(stdout,"UVW (SUN) = %s\n",vec2str(UVW_Sun,"%.5e "));

  //Input
  VPRINT(stdout,"UVW (LSR) = %s\n",vec2str(LSR+3,"%.5e "));

  //Add Solar System Barycenter Velocity
  vadd_c(LSR+3,UVW_Sun,GC+3);
  VPRINT(stdout,"GC (inclined) = %s\n",vec2str(GC+3,"%.5e "));

  //Rotate for zsun
  mxv_c(MGC,GC+3,GC+3);
  VPRINT(stdout,"GC (ICRS) = %s\n",vec2str(GC+3,"%.5e "));

  //Compute cartesian coordinates respect to galactic center
  double r_Sun[]={-ROSUN*PARSEC/1e3,0,0};
  double rup[3];
  VPRINT(stdout,"r (SUN) = %s\n",vec2str(r_Sun,"%.5e "));
  VPRINT(stdout,"r (LSR) = %s\n",vec2str(LSR,"%.5e "));
  vadd_c(LSR,r_Sun,rup);
  VPRINT(stdout,"r (inclined) = %s\n",vec2str(rup,"%.5e "));

  //Rotate for zsun
  mxv_c(MGC,rup,GC);
  VPRINT(stdout,"r (GC) = %s\n",vec2str(GC,"%.5e "));

  return 0;
}

/*
  Transform state vector from cartesian (x,y,z,dx/dt,dy/dt,dz/dt) to
  polar (R,q,z,dR/dt,dq/dt,dz/dt).

  "units" is required to reconcile units of distance with units of
  velocity if they are different.  For instance if you use [R]=pc, and
  [v]=km/s, dvq/dt=R dqdt and dvqdt and R must be converted to the
  same units.  Since dqdt=(dvq/dt)/R, asuming than UV is the factor to
  convert dvq/dt in a common set of units and UL is the analogous for
  radius, then dqdt=[(dvq/dt)/R]_SI * UV/UR.

  Units = (units of velocity)/(unit of distance)

  Example: if rho is in parsec and vrho in km/s units=1e3/PARSEC

  Test:
  double xGC_p[6],units=1e3/PARSEC;
  cart2polar(xGC_s,xGC_p,units);
  VPRINT(stdout,"GC (polar) = %s (pc,km/s)\n",vec2strn(xGC_p,6,"%.17e "));

  polar2cart(xGC_p,xGC_s,units);
  VPRINT(stdout,"GC (cart.) = %s (pc,km/s)\n",vec2strn(xGC_s,6,"%.17e "));
  VPRINT(stdout,"GC (polar) = %s (pc,km/s)\n",vec2strn(xGC_p,6,"%.17e "));
 */
int cart2polar(double x[6],double p[6],double units=1)
{
  double *r=x,*v=x+3;
  double *rho=p,*vrho=p+3;

  //State vector in cartesian
  reccyl_c(r,&rho[0],&rho[1],&rho[2]);
  double cosq=cos(rho[1]),sinq=sin(rho[1]);
  double IM[][3]={{cosq,sinq,0},{-sinq,cosq,0},{0,0,1}};
  mxv_c(IM,v,vrho);
  vrho[1]=vrho[1]/rho[0]*units;
  
  return 0;
}
/*
  Transform state vector from polar (R,q,z,dR/dt,dq/dt,dz/dt) to
  cartesian (x,y,z,dx/dt,dy/dt,dz/dt).

  "units" is required to reconcile units of distance with units of
  velocity if they are different.  For instance if you use [R]=pc, and
  [v]=km/s, dvq/dt=R dqdt and dvqdt and R must be converted to the
  same units.  Since dqdt=(dvq/dt)/R, asuming than UV is the factor to
  convert dvq/dt in a common set of units and UL is the analogous for
  radius, then dqdt=[(dvq/dt)/R]_SI * UV/UR.

  Units = (units of velocity)/(unit of distance)

  Example: if rho is in parsec and vrho in km/s units=1e3/PARSEC
 */
int polar2cart(double p[6],double x[6],double units=1.0)
{
  double *r=x,*v=x+3;
  double *rho=p,*vrho=p+3;
  double cosq=cos(rho[1]),sinq=sin(rho[1]);
  double M[][3]={{cosq,-sinq,0},{sinq,cosq,0},{0,0,1}};

  cylrec_c(rho[0],rho[1],rho[2],r);

  vrho[1]=vrho[1]*rho[0]/units;
  mxv_c(M,vrho,v);
  vrho[1]=vrho[1]/rho[0]*units;

  return 0;
}

int integrateEoM(double tini,double X0[],double h,int npoints,double duration,
		 int nsys,int eom(double,double*,double*,void*),void *params,
		 double *ts,double** X)
{
  VPRINT(stdout,"Entre\n");
  //INTEGRATE GC
  double direction=duration/fabs(duration);
  h*=direction;
  double t_start=tini;
  double t_step=duration/(npoints-1);
  double tend=t_start+duration;
  double t_stop=tend;
  double t=t_start;
  double h_used=h;
  double deltat,h_next,h_adjust;
  int status;
  double *x0=(double*)malloc(nsys*sizeof(double));
  double *x=(double*)malloc(nsys*sizeof(double));

  VPRINT(stdout,"Integration parameters:\n");
  VPRINT(stdout,"\tnsys = %d\n",nsys);
  VPRINT(stdout,"\ttini = %.5e UT\n",tini);
  VPRINT(stdout,"\tt_start = %.5e UT\n",t_start);
  VPRINT(stdout,"\tt_stop = %.5e UT\n",t_stop);
  VPRINT(stdout,"\tt_step = %.5e UT\n",t_step);
  VPRINT(stdout,"\ttend = %.5e UT\n",tend);
  VPRINT(stdout,"\th = %.5e UT\n",h);
  VPRINT(stdout,"\tt = %.5e UT\n",t);
  VPRINT(stdout,"\th_used = %.5e UT\n",h_used);
  VPRINT(stdout,"\tDirection = %.0lf UT\n",direction);
  if(VERBOSE) getchar();
  
  //INITIAL CONDITIONS
  copyVec(x0,X0,nsys);
  copyVec(x,x0,nsys);
  VPRINT(stdout,"\ty = %s\n",vec2strn(x0,nsys,"%.7e "));

  //INTEGRATE
  char qpause=0;
  for(int i=0;i<npoints;i++) {
    ts[i]=t;
    copyVec(X[i],x0,nsys);
    deltat=t-tini;
    if(direction*((t_start+t_step)-tend)>0) t_step=(tend-t_start);
    t_stop=t_start+t_step;
    h_used=h;
    VPRINT(stdout,"Step %d: t = %.5e UT\n",i,t);
    VPRINT(stdout,"\tt_start = %.5e UT\n",t_start);
    VPRINT(stdout,"\tt_stop = %.5e UT\n",t_stop);
    VPRINT(stdout,"\th_used = %.5e UT\n",h_used);
    VPRINT(stdout,"\ty0 = %s\n",vec2strn(x0,nsys,"%.7e "));
    VPRINT(stdout,"\ty = %s\n",vec2strn(x,nsys,"%.7e "));

    int nmax=(int)(10*fabs(t_step/h));
    int nint=0;
    do{
      while(1){
	status=Gragg_Bulirsch_Stoer(eom,x0,x,t,h_used,&h_next,1.0,
				    TOLERANCE,EXTMET,params);

	VPRINT(stdout,"\t\tStatus = %d, h_used = %e, h_next = %e, x = %s\n",
	       status,h_used,h_next,vec2strn(x,nsys/2,"%.17e,"));
	if(VERBOSE) getchar();
	
	if(status){
	  //if(fabs(h_next)>fabs(h_used)) h_used=h_next;
	  //else 
	  h_used/=4.0;
	  VPRINT(stdout,"\t\t\tAdaptando paso h_used = %e, h_next = %e\n",h_used,h_next);
	  if(direction*h_used<0){
	    fprintf(stdout,"************ERROR*********** TIME STEP HAS THE WRONG SIGN\n");
	    h_used=h;
	    //exit(0);
	  }
	  if(VERBOSE) getchar();
	  qpause=1;
	}
	else{
	  break;
	}
      }
      VPRINT(stdout,"\t\tSali del condenado ciclo con h_used = %e, h_next = %e\n",h_used,h_next);

      t+=h_used;
      copyVec(x0,x,nsys);

      VPRINT(stdout,"\t\tLlegue a t = %.17e (t_stop = %.17e)\n",t,t_stop);

      if(direction*(t+h_next-t_stop)>0){
	//h_used=t+h_next-t_stop;
	h_used=h_next+(t_stop-(t+h_next));
	VPRINT(stdout,"\t\tUso lo que me falta: %e\n",h_used);
	if(VERBOSE) getchar();
      }
      else{
	h_used=h_next;
	VPRINT(stdout,"\t\tUso h_next\n");
      }

      VPRINT(stdout,"\t\tEn el siguiente paso usaré h = %e\n",h_used);
      if(VERBOSE && qpause) getchar();
      
      nint++;
    }while(direction*(t-(t_stop-direction*fabs(t_step)*1.e-7))<0 && nint<nmax);
    if(nint==nmax){
      free(x0);
      free(x);
      throw(1);
    }

    VPRINT(stdout,"\tComplete ese condenado intervalo\n");

    if(direction*(t-t_stop)>0){
      h_adjust=(t_stop-t);
      status=Gragg_Bulirsch_Stoer(eom,x0,x,t,h_adjust,&h_next,1.0,
				  TOLERANCE,EXTMET,params);
      copyVec(x0,x,nsys);
      t=t_stop;
    }
    t_start = t;
    if(direction*(t_start-tend)>0) break;
  }
  
  free(x0);
  free(x);
  return 0;
}

/*
  Equations as given by:
  https://www.aanda.org/articles/aa/pdf/2001/44/aah2819.pdf

  Variables are:

  y[0]:R
  y[1]:q (q angular coordinate, theta, phi)
  y[2]:z

  (if there are more than 1 particle, y[6],y[7],y[8] will be the
  corresponding coordinates of particle 2 and so on)

  y[3]:dRDt
  y[4]:dqdt 
  y[5]:dzdt

  (same as before if there is several particles)
 */
int EoMGalactic(double t,double y[],double dydt[],void *params) 
{ 
  double* ps=(double*)params;
  int nsys=(int)ps[0];
  int Npart=(int)(nsys/6);

  double R,q,z,vR,uq,vz,dphidR,dphidq,dphidz;
  
  //VPRINT(stdout,"Number of particles: %d (sys. %d)\n",Npart,nsys);
  //VPRINT(stdout,"\ty = %s\n",vec2strn(y,nsys,"%.17e "));
  //fprintf(stdout,"\ty = %s\n",vec2strn(y,nsys,"%.17e "));

  //dydt=vy
  for(int i=Npart;i-->0;){
    int ip=6*i;
    int j=ip;
    dydt[j]=y[3+j];j++;
    dydt[j]=y[3+j];j++;
    dydt[j]=y[3+j];j++;
  }

  //dvdt=f(y,vy,t)
  for(int i=Npart;i-->0;){
    int ip=6*i;

    //POTENTIAL
    int j=ip;
    R=y[j];j++;
    q=y[j];j++;
    z=y[j];j++;
    vR=y[j];j++;
    uq=y[j];j++;
    vz=y[j];j++;

    gradGalacticPotential(R,q,z,&dphidR,&dphidq,&dphidz,ps);
    /*
    fprintf(stdout,"\ty = %s\n",vec2strn(y+ip,6,"%.17e "));
    fprintf(stdout,"Particle %d:\n\tR = %.17e, z = %.17e, dphidR = %.17e, dphidq = %.17e, dphidz = %.17e\n",
	    i,R,z,dphidR,dphidq,dphidz);
    exit(0);
    //*/

    int k=ip+3;
    /*
    //LMA
    dydt[k]=0;k++;
    dydt[k]=0;k++;
    dydt[k]=0;k++;
    //*/
    //*
    //LMA POLAR
    //dphidR=0.0;dphidq=0.0;dphidz=0.0;
    dydt[k]=-dphidR+R*uq*uq;k++;
    dydt[k]=-dphidq-2*vR*uq/R;k++;
    dydt[k]=-dphidz;k++;
    //*/

    //*
    VPRINT(stdout,"\tR = %e, vR = %e, q = %e, uq = %e, z = %e, vz = %e, dphidR = %e, dphidq = %e, dphidz = %e, nsys = %d\n",
	   R,vR,q,uq,z,vz,dphidR,dphidq,dphidz,nsys);
    VPRINT(stdout,"\ty = %s\n\tdydt = %s\n",
	   vec2strn(y,nsys,"%.17e "),
	   vec2strn(dydt,nsys,"%.17e "));
    if(VERBOSE) getchar();
  }

  //*/
  //exit(0);
  //*/
  return 0;
}

int findTime(double t,double tsp[],int Ntimesp)
{
  double sgn=t/fabs(t);
  for(int i=Ntimesp-1;i-->0;){
    if((t-tsp[i])*(t-tsp[i+1])<0) return i;
  }
  return -1;
}

/*
  Compute the distance from a point to an interval

  Perpendicular distance :
  Formula in: 
  https://math.stackexchange.com/questions/1300484/distance-between-line-and-a-point
 */
double distancePointLine(double p[],double p1[],double p2[],double *dtfrac)
{
  double p2mp1[3],pmp1[3],pmp2[3],pm[3],d[3],p2mp1a,pma,dmin,d1,d2;

  //Useful vectors
  vsub_c(p2,p1,p2mp1);
  vsub_c(p,p1,pmp1);
  vsub_c(p,p2,pmp2);
  p2mp1a=vnorm_c(p2mp1);

  //Compute perpendicular distance
  vscl_c(vdot_c(pmp1,p2mp1)/(p2mp1a*p2mp1a),p2mp1,pm);
  vsub_c(pmp1,pm,d);
  pma=vnorm_c(d);
  *dtfrac=pm[0]/p2mp1[0];

  if(*dtfrac>=0 && *dtfrac<=1){
    VPRINT(stdout,"Inside distance (dtfrac = %e)\n",*dtfrac);
    dmin=pma;
  }
  else{
    VPRINT(stdout,"Extreme distance (d1=%e,d2=%e)\n",vnorm_c(pmp1),vnorm_c(pmp2));
    if(*dtfrac<0){
      dmin=vnorm_c(pmp1);
      *dtfrac=0;
    }else{
      dmin=vnorm_c(pmp2);
      *dtfrac=0;
    }
  }

  return dmin;
}

int generateMultivariate(double** cov,double *mu,double **x,int n,int Nobs)
{
  gsl_matrix *C=gsl_matrix_alloc(n,n);
  gsl_vector *muv=gsl_vector_alloc(n);
  gsl_vector *xv=gsl_vector_alloc(n);

  for(int i=n;i-->0;){
    gsl_vector_set(muv,i,mu[i]);
    for(int j=n;j-->0;)
      gsl_matrix_set(C,i,j,cov[i][j]);
  }

  //PREPARE
  gsl_linalg_cholesky_decomp1(C);
  
  //GENERATE
  for(int i=0;i<Nobs;i++){
    gsl_ran_multivariate_gaussian(RAND,muv,C,xv);
    for(int j=n;j-->0;) x[i][j]=gsl_vector_get(xv,j);
  }
    
  //RETURN
  gsl_matrix_free(C);
  gsl_vector_free(muv);
  gsl_vector_free(xv);
  
  return 0;
}

double terminalDistance(double t,void *params)
{
  double* ps=(double*)params;
  double **x1=matrixAllocate(2,6),**x2=matrixAllocate(2,6);
  double* ipars;
  double ts[2];
  double h=fabs(t)/100;
  double dx[6],d,dv;
  ipars=ps+12;

  copyVec(x1[0],ps,6);
  copyVec(x2[0],ps+6,6);

  VPRINT(stdout,"Attempting to integrate with tint = %.6e\n",t);
  VPRINT(stdout,"Parameters: %s\n",vec2strn(ipars,10,"%.5e "));

  VPRINT(stdout,"Initial conditions:\n\t%s\n\t%s\n",
	  vec2strn(x1[0],6,"%.5e "),vec2strn(x2[0],6,"%.5e "));

  //Integrate star
  integrateEoM(0,ps,h,2,t,6,EoMGalactic,ipars,ts,x1);
  VPRINT(stdout,"Integration result star: %s\n",vec2strn(x1[1],6,"%.5e "));
  
  //Integrate particle
  integrateEoM(0,ps+6,h,2,t,6,EoMGalactic,ipars,ts,x2);
  VPRINT(stdout,"Integration result particle: %s\n",vec2strn(x2[1],6,"%.5e "));

  //Distance
  vsubg_c(x1[1],x2[1],6,dx);
  d=vnorm_c(dx);
  dv=vnorm_c(dx+3);
  VPRINT(stdout,"Distance = %e, Velocity difference = %e\n",d,dv);
  
  if(VERBOSE) getchar();
  
  if((int)ps[22]==1){
    /*
    for(int i=6;i-->0;){
      ps[i]=x1[1][i];
      ps[6+i]=x2[1][i];
    }
    */
    copyVec(ps,x1[1],6);
    copyVec(ps+6,x2[1],6);
  }
  
  return d;
}

#define MVERBOSE 0
int minDistanceDiscrete(double *xs,double **xp,double *tsp,int Ntimes,double duration,double *params,
			double *dmin,double *tmin)
{
  //INTEGRATION PARAMETERS
  double hstep=fabs(duration)/(10*Ntimes);
  double *ts=(double*)malloc(Ntimes*sizeof(double));
  double **xInt=matrixAllocate(Ntimes,6);
  double x[6],*xp1,*xp2,dx[6],xpmin[6];
  double t1,t2,ftmin,dyn_vrel;
  
  //INTEGRATE STAR
  if(MVERBOSE) fprintf(stdout,"Integrating star until %e:\n",duration);
  integrateEoM(0,xs,hstep,Ntimes,duration,
	       6,EoMGalactic,params,
	       ts,xInt);

  //FIND TIME OF MINIMUM DISTANCE
  double dyn_dmin=1.0e+100,dyn_tmin;
  for(int i=1;i<Ntimes;i++){

    //STAR POSITION 
    polar2cart(xInt[i],x,1.0);
    if(MVERBOSE) fprintf(stdout,"\tPosition %d @ t=%e:%s\n",i,ts[i],vec2strn(x,6,"%.5e,"));

    //FIND TIME IN PARTICLE INTEGRATION
    double t=ts[i];
    int it=findTime(t,tsp,Ntimes);
    if(MVERBOSE) fprintf(stdout,"\t\tCorrespond to interval %d: [%e,%e]\n",it,tsp[it],tsp[it+1]);

    //FIND POSITION OF NOMINAL PARTICLE AT INITIAL POINT OF INTERVAL
    xp1=xp[it];
    xp2=xp[it+1];
    if(MVERBOSE) fprintf(stdout,"\t\tNominal particle position 1 (%e):%s\n",tsp[it],vec2strn(xp1,6,"%.5e "));
    if(MVERBOSE) fprintf(stdout,"\t\tNominal particle position 2 (%e):%s\n",tsp[it+1],vec2strn(xp2,6,"%.5e "));

    //DIFFERENCE
    vsubg_c(x,xp1,6,dx);
    if(MVERBOSE) fprintf(stdout,"\t\tDifference 1: [%s] (pos:%.5e,vel:%.5e]\n",
	   vec2strn(dx,6,"%.5e,"),vnorm_c(dx),vnorm_c(dx+3)*UV);
    vsubg_c(x,xp2,6,dx);
    if(MVERBOSE) fprintf(stdout,"\t\tDifference 2: [%s] (pos:%.5e,vel:%.5e]\n",
	   vec2strn(dx,6,"%.5e,"),vnorm_c(dx),vnorm_c(dx+3)*UV);

    //CALCULATE DISTANCE FROM PARTICLE TO INTERVAL
    *dmin=distancePointLine(x,xp1,xp2,&ftmin);
    if(MVERBOSE) fprintf(stdout,"\t\tPosition factor: %.5e\n",ftmin);
    if(MVERBOSE) fprintf(stdout,"\t\tDistance at minimum: %.5e\n",*dmin);
      
    //CORRECTED TIME
    t1=tsp[it]<tsp[it+1]?tsp[it]:tsp[it+1];
    t2=tsp[it]>tsp[it+1]?tsp[it]:tsp[it+1];
    *tmin=t1+ftmin*(t2-t1);
    if(MVERBOSE) fprintf(stdout,"\t\tCorrected time at minimum:%e\n",*tmin);

    //INTERPOLATED POSITION OF NOMINAL PARTICLE
    vsubg_c(xp2,xp1,6,dx);
    vscl_c((*tmin-tsp[it])/(tsp[it+1]-tsp[it]),dx,dx);
    vscl_c((*tmin-tsp[it])/(tsp[it+1]-tsp[it]),dx+3,dx+3);
    vaddg_c(xp1,dx,6,xpmin);

    //DISTANCE TO STAR
    vsubg_c(xpmin,x,6,dx);
    if(MVERBOSE) fprintf(stdout,"\t\tDifference at minimum: [%s] (pos:%.5e,vel:%.5e]\n",
	   vec2strn(dx,6,"%.5e,"),vnorm_c(dx),vnorm_c(dx+3)*UV/1e3);

    //COMPUTE MINIMUM DISTANCE
    if(*dmin<=dyn_dmin){
      dyn_dmin=*dmin;
      dyn_tmin=*tmin;
      dyn_vrel=vnorm_c(dx+3)*UV/1e3;
    }
  }

  fprintf(stdout,"\tDynamical minimum distance: %e\n",dyn_dmin);
  fprintf(stdout,"\tDynamical minimum time: %e\n",dyn_tmin);
  return 0;
}


int minDistanceBad(double *xs,double *xp,double mint,double maxt,double tmin0,
		double *dmin,double *tmin,double *params)
{
  int status;
  double ps[23];

  copyVec(ps,xs,6);
  copyVec(ps+6,xp,6);
  copyVec(ps+12,params,10);
  //LAST PARAMETERS INSTRUCT THE ROUTINES TO RETURN THE FINAL STATE
  ps[22]=0.0;

  VPRINT(stdout,"Initial conditions:%s\n",vec2strn(ps,12,"%.5e "));
  VPRINT(stdout,"Other parameters:%s\n",vec2strn(ps+12,10,"%.5e "));

  *dmin=terminalDistance(*tmin,ps);

  gsl_function F={.function=&terminalDistance,.params=ps};
  gsl_min_fminimizer *s=gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
  gsl_min_fminimizer_set(s,&F,tmin0,mint,maxt);

  //MINIMIZE
  int i=0;
  do{
    i++;
    status=gsl_min_fminimizer_iterate(s);
    *tmin=gsl_min_fminimizer_x_minimum(s);
    mint=gsl_min_fminimizer_x_lower(s);
    maxt=gsl_min_fminimizer_x_upper(s);
    status=gsl_min_test_interval(mint,maxt,0.0,1e-5);
    if(status==GSL_SUCCESS)
      VPRINT(stdout,"Minimization converged after %d iterations\n",i);
  }while(status==GSL_CONTINUE && i<100);

  //COMPUTE STATE AT MINIMIZE
  ps[22]=1;
  *dmin=terminalDistance(*tmin,ps);
  VPRINT(stdout,"Final conditions:%s\n",vec2strn(ps,12,"%.5e "));
  copyVec(xs,ps,6);
  copyVec(xp,ps+6,6);
  
  return 0;
}

/*
  Use multidimensional minimization methods to avoid guessing
 */
double terminalDistance2(const gsl_vector *x,void *params)
{
  double t;
  double* ps=(double*)params;
  double **x1=matrixAllocate(2,6),**x2=matrixAllocate(2,6);
  double *x1c=vectorAllocate(6);
  double *x2c=vectorAllocate(6);
  double *dx=vectorAllocate(6);
  double* ipars;
  double ts[2];
  double h;
  double d,dv;
  ipars=ps+12;

  t=gsl_vector_get(x,0);
  //fprintf(stdout,"%e\n",t);
  h=fabs(t)/100;
  copyVec(x1[0],ps,6);
  copyVec(x2[0],ps+6,6);

  VPRINT(stdout,"Attempting to integrate with tint = %.6e\n",t);
  VPRINT(stdout,"Parameters: %s\n",vec2strn(ipars,10,"%.5e "));

  VPRINT(stdout,"Initial conditions:\n\t%s\n\t%s\n",
	  vec2strn(x1[0],6,"%.5e "),vec2strn(x2[0],6,"%.5e "));

  //Integrate star
  integrateEoM(0,ps,h,2,t,6,EoMGalactic,ipars,ts,x1);
  VPRINT(stdout,"Integration result star: %s\n",vec2strn(x1[1],6,"%.5e "));

  //Convert to cartesian
  polar2cart(x1[1],x1c);

  //Integrate particle
  h=1e4;
  integrateEoM(0,ps+6,h,2,t,6,EoMGalactic,ipars,ts,x2);
  VPRINT(stdout,"Integration result particle: %s\n",vec2strn(x2[1],6,"%.5e "));

  //Convert to cartesian
  polar2cart(x2[1],x2c);

  //Distance
  vsubg_c(x1c,x2c,6,dx);
  d=vnorm_c(dx);
  dv=vnorm_c(dx+3);
  VPRINT(stdout,"Distance = %e, Velocity difference = %e\n",d,dv);

  if(VERBOSE) getchar();
  
  if((int)ps[22]==1){
    copyVec(ps,x1[1],6);
    copyVec(ps+6,x2[1],6);
  }
  
  freeMatrix(x1,2,6);
  freeMatrix(x2,2,6);
  free(x1c);
  free(x2c);
  free(dx);
  return d;
}

int minDistance(double *xs,double *xp,double tmin0,
		 double *dmin,double *tmin,double *params)
{
  int status;
  double *ps=vectorAllocate(23);
  double size;
  double error=fabs(tmin0)/10;
  double tolerance=error/10;

  copyVec(ps,xs,6);
  copyVec(ps+6,xp,6);
  copyVec(ps+12,params,10);
  //LAST PARAMETERS INSTRUCT THE ROUTINES TO RETURN THE FINAL STATE
  ps[22]=0.0;

  VPRINT(stdout,"Initial conditions:%s\n",vec2strn(ps,12,"%.5e "));
  VPRINT(stdout,"Other parameters:%s\n",vec2strn(ps+12,10,"%.5e "));
  gsl_multimin_function F;
  F.n=1;
  F.f=&terminalDistance2;
  F.params=ps;
  gsl_multimin_fminimizer *s=
    gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,1);
  gsl_vector *x=gsl_vector_alloc(1);
  gsl_vector *dx=gsl_vector_alloc(1);
  
  gsl_vector_set(x,0,tmin0);
  gsl_vector_set_all(dx,fabs(tmin0)/10);
  gsl_multimin_fminimizer_set(s,&F,x,dx);

  //MINIMIZE
  int i=0;
  do{
    i++;
    VPRINT(stdout,"Iteration %d\n",i);
    status=gsl_multimin_fminimizer_iterate(s);
    if(status) break;
    size=gsl_multimin_fminimizer_size(s);
    status=gsl_multimin_test_size(size,tolerance);
    if(status==GSL_SUCCESS)
      VPRINT(stdout,"Minimization converged after %d iterations\n",i);
  }while(status==GSL_CONTINUE && i<100);

  //GET RESULT
  *tmin=gsl_vector_get(s->x,0);
  *dmin=s->fval;

  //COMPUTE STATE AT MINIMIZE
  ps[22]=1;
  gsl_vector_set(x,0,*tmin);
  *dmin=terminalDistance2(x,ps);
  VPRINT(stdout,"Final conditions:%s\n",vec2strn(ps,12,"%.5e "));
  copyVec(xs,ps,6);
  copyVec(xp,ps+6,6);
  
  free(ps);
  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(x);
  gsl_vector_free(dx);

  return 0;
}

double wFunction2(double d,void *params)
{
  double h=*(double*)params;
  double q=d/h;
  double w;

  w=1/(h*h*h)*exp(-d*d/(h*h));

  return w;
}

double wFunction(double d,void *params)
{
  double h=*(double*)params;
  double q=d/h;
  double w;
  if(q<1)
    w=0.25*(2-q)*(2-q)*(2-q)-(1-q)*(1-q)*(1-q);
  else if(q<2)
    w=0.25*(2-q)*(2-q)*(2-q);
  else
    w=0;
  return w;
}

double wNormalization(double h)
{
  double norm,error;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  gsl_function F={.function=&wFunction,.params=&h};
  gsl_integration_qags(&F,0.0,2*h,0.0,1e-7,1000,w,&norm,&error);
  return 1/norm;
}

double wNormalization2(double h)
{
  double norm,error;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  gsl_function F={.function=&wFunction2,.params=&h};
  gsl_integration_qagiu(&F,0.0,0.0,1e-7,1000,w,&norm,&error);
  return 1/norm;
}

/*
  Minimum and maximum using the Tournament Method

  Adapted from: http://www.geeksforgeeks.org/maximum-and-minimum-in-an-array/

  Number of comparisons: 3n/2 - 2
  (Linear search takes 2n-3)
 */
int getMinMax(double array[],int left,int right,double *min,double *max,
	      int stride=0,int period=1)
{
  int mid;
  double minleft,maxleft;
  double minright,maxright;
    
  *max = array[left*period+stride];
  *min = array[left*period+stride];

  if(right == left)
    return 0; 
  mid = (left + right)/2;  

  getMinMax(array, left, mid, &minleft, &maxleft, stride, period);
  getMinMax(array, mid+1, right, &minright, &maxright, stride, period);  
    
  if (maxleft > maxright)
    *max = maxleft;
  else
    *max = maxright;    
     
  /* Take the minimum of both sub array */
  if (minleft < minright)
    *min = minleft;
  else
    *min = minright;     

  return 0;
}

int centerCloud(double **xp,int Npart,int Ncoord,
		double *xc,double **xg,double *r)
{
  for(int j=0;j<Ncoord;j++){
    double xm=0.0;
    for(int i=0;i<Npart;i++){
      xm+=xp[i][j];
    }
    xc[j]=xm/Npart;
  }
  return 0;
}

//Taken from: https://en.wikiversity.org/wiki/C_Source_Code/Find_the_median_and_mean
int getPercentile(double x[],int n,double p/*p-value:0-1*/,double nominal,
		  double *min,double *max,double *f,
		  double *ql,double *qm,double *qu) 
{
  double temp;
  int i, j;

  for(i=0; i<n-1; i++) {
    for(j=i+1; j<n; j++) {
      if(x[j] < x[i]) {
	temp = x[i];
	x[i] = x[j];
	x[j] = temp;
      }
    }
  }
  
  //FRACTION BELOW NOMINAL
  for(i=0;i<n;i++){
    if(x[i]>nominal) break;
  }
  *f=(1.0*i)/n;

  //MEDIAN
  if(n%2==0) {
    *qm=((x[n/2] + x[n/2 - 1]) / 2.0);
  } else {
    *qm=x[n/2];
  }

  //LOWER
  int nl=(int)floor((1-p)*n);
  *ql=(x[nl]+x[nl+1])/2;
  nl=nl==(n-1)?(n-2):nl;
  
  //UPPER
  int nu=(int)ceil(p*n);
  nu=nu==0?1:nu;
  *qu=(x[nu]+x[nu-1])/2;

  //MINIMUM AND MAXIMUM
  *min=x[0];
  *max=x[n-1];
}

int cloudVolume(double *x,int Npart,double *dV)
{
  /*
    Better implementation is a convex hull

    See:
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
  */

  int ip,jp;
  double *xi,*xj;

  double *xci=vectorAllocate(6);
  double *xcj=vectorAllocate(6);
  double *dx=vectorAllocate(3);
  double *ds=vectorAllocate(Npart);
  double **dd=matrixAllocate(Npart,Npart);
  double dmin;

  //Calculate the average distance
  double ad=0.0;
  for(int i=0;i<Npart;i++){
    ip=6*i;
    xi=x+ip;
    polar2cart(xi,xci);
    for(int j=0;j<Npart;j++){
      jp=6*j;
      xj=x+jp;
      ds[j]=1e100;
      if(i==j) continue;
      if(j<i) ds[j]=dd[j][i];
      else{
	polar2cart(xj,xcj);
	vsub_c(xcj,xci,dx);
	ds[j]=dd[i][j]=vnorm_c(dx);
      }
    }
    dmin=gsl_stats_min(ds,1,Npart);
    ad+=dmin;
  }
  ad/=Npart;

  free(xci);
  free(xcj);
  free(dx);
  free(ds);
  freeMatrix(dd,Npart,Npart);

  *dV=(4*M_PI/3*ad*ad*ad);
  return 0;
}

int cloudProperties2(double *x,int Npart,
		     double *r90,double *v90)
{
  int ip,is;
  double *xadius=vectorAllocate(6);
  double *center=vectorAllocate(6);
  double *radius=vectorAllocate(Npart);
  double *vadius=vectorAllocate(Npart);
  double *xc=vectorAllocate(6);
  double *xp;
  double min,max,ql,qm,qu;
  double nom,f;

  //Get the center
  for(int j=0;j<6;j++){
    double xm=0.0;
    is=0;
    for(int i=Npart;i-->0;){
      ip=6*i;
      //Particle
      xp=x+ip;
      if(xp[0]==99.99) continue;
      //Convert to cartesian
      polar2cart(xp,xc);
      xm+=xc[j];
      is++;
    }
    center[j]=xm/is;
  }
  
  //Calculate the radius-vector
  is=0;
  for(int i=Npart;i-->0;){
    ip=6*i;
    xp=x+ip;
    if(xp[0]==99.99) continue;
    //Convert to cartesian
    polar2cart(xp,xc);
    //Compute radius
    vsubg_c(xc,center,6,xadius);
    radius[is]=vnorm_c(xadius);
    vadius[is]=vnorm_c(xadius+3);
    is++;
  }

  //Calculate the percentiles
  nom=0;
  getPercentile(radius,is,0.9,nom,&min,&max,&f,&ql,&qm,&qu);
  *r90=qu;
  getPercentile(vadius,is,0.9,nom,&min,&max,&f,&ql,&qm,&qu);
  *v90=qu;
  
  free(radius);
  free(vadius);
  free(xadius);
  free(center);
  free(xc);
}

#define CVERBOSE 0
int cloudProperties(double *xp,int Npart,
		    double *radius,double *vradius,
		    double *rinter,double *vinter)
{
  double Rmin,Rmax,qmin,qmax,zmin,zmax;
  double *x=vectorAllocate(6);
  double *xg=vectorAllocate(6);
  double *dx=vectorAllocate(6);
  double *xc=vectorAllocate(6*Npart);

  //Get ranges of variables
  getMinMax(xp,0,Npart-1,&Rmin,&Rmax,0,6);
  getMinMax(xp,0,Npart-1,&qmin,&qmax,1,6);
  getMinMax(xp,0,Npart-1,&zmin,&zmax,2,6);

  //Convert ranges to variables
  x[0]=Rmin;x[1]=Rmin*cos(qmin);x[2]=zmin;
  xg[0]=Rmax;xg[1]=Rmax*cos(qmax);xg[2]=zmax;
  if(CVERBOSE) fprintf(stdout,"\t\tR=(%e,%e),q=(%e,%e),z=(%e,%e)\n",
	  Rmin,Rmax,qmin,qmax,zmin,zmax);
  vsub_c(x,xg,dx);
  *radius=vnorm_c(dx)/2;
  
  //CONVERSION TO CARTESIAN COORDINATES
  int ip;
  for(int i=0;i<Npart;i++){
    ip=6*i;
    polar2cart(xp+ip,xc+ip);
  }
  if(CVERBOSE) fprintf(stdout,"%s\n",vec2strn(xc,6*Npart,"%e "));

  //VELOCITY RADIUS
  double vxmin,vxmax,vymin,vymax,vzmin,vzmax;
  getMinMax(xc,0,Npart-1,&vxmin,&vxmax,3,6);
  getMinMax(xc,0,Npart-1,&vymin,&vymax,4,6);
  getMinMax(xc,0,Npart-1,&vzmin,&vzmax,5,6);
  if(CVERBOSE) fprintf(stdout,"\t\tvx=(%e,%e),vy=(%e,%e),vz=(%e,%e)\n",
	  vxmin,vxmax,vymin,vymax,vzmin,vzmax);
  x[0]=vxmin;x[1]=vymin;x[2]=vzmin;
  xg[0]=vxmax;xg[1]=vymax;xg[2]=vzmax;
  vsub_c(x,xg,dx);
  *vradius=vnorm_c(dx)/2;

  //INTEROBJECT SEPARATION IN SPACE
  int jp;
  double *xi,*xj;
  double dmean=0,dstd=0,vmean=0,vstd=0,dminmean=0,vminmean=0,dint,vint;
  int k=0;
  for(int i=Npart;i-->0;){
    ip=6*i;
    xi=xc+ip;
    double dmin=1e100,vmin=1e100;
    for(int j=Npart;j-->0;){
      if(j==i) continue;
      jp=6*j;
      xj=xc+jp;
      vsubg_c(xi,xj,6,dx);
      dint=vnorm_c(dx);
      vint=vnorm_c(dx+3);
      dmean+=dint;
      vmean+=vint;
      dstd+=dint*dint;
      vstd+=vint*vint;
      dmin=dint<dmin?dint:dmin;
      vmin=vint<vmin?vint:vmin;
      k++;
    }
    dminmean+=dmin;
    vminmean+=vmin;
    if(CVERBOSE) fprintf(stdout,"Minimum distance to %d:%e\n",i,dmin);
  }
  dminmean/=Npart;
  vminmean/=Npart;
  if(CVERBOSE) fprintf(stdout,"Average distance between particles:%e\n",dminmean);

  *rinter=dminmean;
  *vinter=vminmean;

  free(x);
  free(xg);
  free(dx);
  free(xc);
  return 0;
}

int printHeader(FILE* stream,const char *msg,char mark='*',int ntabs=0,int nmax=60)
{
  char *bar=charVectorAllocate(1000);
  char *tabs=charVectorAllocate(100);
  int n=MAX(strlen(msg),nmax);

  sprintf(tabs,"");
  if(ntabs==1) sprintf(tabs,"\t");
  if(ntabs==2) sprintf(tabs,"\t\t");
  if(ntabs==3) sprintf(tabs,"\t\t\t");
  
  sprintf(bar,"%s",tabs);
  for(int i=n;i-->0;) sprintf(bar,"%s%c",bar,mark);

  fprintf(stream,"%s\n%s%s\n%s\n",bar,tabs,msg,bar);

  free(bar);
  free(tabs);
  return 0;
}

//Taken from: https://en.wikiversity.org/wiki/C_Source_Code/Find_the_median_and_mean
int quantilesVector(double x[],int n,double *min,double *max,
		    double *q050,double *q005,double *q095) 
{
  double temp;
  int i, j;

  for(i=0; i<n-1; i++) {
    for(j=i+1; j<n; j++) {
      if(x[j] < x[i]) {
	temp = x[i];
	x[i] = x[j];
	x[j] = temp;
      }
    }
  }

  //MIN PERC

  //MEDIAN
  if(n%2==0) {
    *q050=((x[n/2] + x[n/2 - 1]) / 2.0);
  } else {
    *q050=x[n/2];
  }

  //LOWER
  int nl=(int)floor(0.05*n);
  *q005=(x[nl]+x[nl+1])/2;
  nl=nl==(n-1)?(n-2):nl;
  
  //UPPER
  int nu=(int)ceil(0.95*n);
  nu=nu==0?1:nu;
  *q095=(x[nu]+x[nu-1])/2;

  //MINIMUM AND MAXIMUM
  *min=x[0];
  *max=x[n-1];
}

struct vinfpar {
  double xmin,xmax;
  gsl_interp_accel *a;
  gsl_spline *s;
};

double vinfPosterior(double x,void *params)
{
  struct vinfpar *pars=(struct vinfpar*)params;
  double p;
  if(x<pars->xmin) x=1.01*pars->xmin;
  if(x>pars->xmax){
    p=12*pow(10.0,-3*x);
    /*
    VPRINT(stdout,"Extrapolation for %e: %e\n",x,p);
    p=gsl_spline_eval(pars->s,pars->xmax,pars->a);
    VPRINT(stdout,"Compare with maximum: %e\n",p);
    */
  }else{
    p=gsl_spline_eval(pars->s,x,pars->a);
  }
  return p;
}

double vinfProbability(double v1,double v2,void *params)
{
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
  gsl_function F={.function=vinfPosterior,.params=params};
  double vint,vint_error;
  gsl_integration_qags(&F,v1,v2,0,1e-5,1000,w,&vint,&vint_error);
  gsl_integration_workspace_free(w);
  return vint;
}

double starDistance(double par)
{
  double d;
  d=AU/tan(par/(60*60*1000.0)*DEG)/PARSEC;
  return d;
}

/*
  Compute the solid angle subtended by a set of vectors in space.
 */
double solidAngle(double** vecs,int nvec,double *vmed,double *vstd)
{
  double fill=0.7861;//Experimental
  double q,f;
  double v,d,ad,dmin;

  double *vs=vectorAllocate(nvec);
  double *ds=vectorAllocate(nvec);
  double *ls=vectorAllocate(nvec);
  double *fs=vectorAllocate(nvec);
  double **dd=matrixAllocate(nvec,nvec);

  //Convert to spherical
  for(int i=0;i<nvec;i++){
    recsph_c(vecs[i],&vs[i],&q,&f);
    fs[i]=f<0?2*M_PI+f:f;
    ls[i]=M_PI/2-q;
  }
  *vmed=gsl_stats_mean(vs,1,nvec);
  *vstd=gsl_stats_sd(vs,1,nvec);

  //Compute the average distance between points
  ad=0.0;
  for(int i=0;i<nvec;i++){
    for(int j=0;j<nvec;j++){
      ds[j]=1e100;
      if(i==j) continue;
      if(j<i) ds[j]=dd[j][i];
      else ds[j]=dd[i][j]=greatCircleDistance(fs[i],fs[j],ls[i],ls[j]);
    }
    dmin=gsl_stats_min(ds,1,nvec);
    ad+=dmin;
  }
  ad/=nvec;

  //Solid angle approximation
  double sa=nvec*M_PI*ad*ad/fill;
  
  free(ds);
  free(vs);
  free(ls);
  free(fs);
  freeMatrix(dd,nvec,nvec);

  return sa;
}

