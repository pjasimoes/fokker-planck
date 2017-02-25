/*

  .:: Time-dependent Fokker-Planck code ::.

  Solves the electron distribution function in energy, cosine of
  pitch-angle, position and time f(E,u,s,t) considering collisional 
  energy loss and pitch-angle scattering and magnetic mirroring.

  For details on input and output see IDL wrapper procedure
  fp.pro

  Reference: Hamilton, Lu and Petrosian,  Apj 354, 726, 1990
             (see for details about the equation and numerical schemes)

  History:
  - original FORTRAN code (by R. Hamilton and E. Lu) (pre-1997)
  - speed improvment by J. McTiernan (1997) and L. Ofman (2001)
  - translated to C/C++, IDL wrapper, improved injection functions, 
    new models for plasma density and magnetic field, allow calls 
    for specific terms of the equation, binary output, progress bar
    by Paulo Simoes (2011-2012)
  - added pitch-angle scattering term ruled by scattering mean free 
    path, P.Simoes, 2012
     
  Files:
  fp.cpp (this file)
  gauleg.h (Gauss-Legendre integration function)
  makefile (compiles with g++)
  fp.pro (IDL wrapper)

  Paulo Simoes, Glasgow, 2011-2012
    paulo.simoes@glasgow.ac.uk
    pjasimoes@gmail.com

*/
#include <iostream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <string>
#include <fstream>
#include <cstdlib>
#include <time.h>

#ifdef _OPENMP
 #include <omp.h>
#endif

#include "gauleg.h"

#define C   2.9979246e+10
#define E0  510.99891
#define EPS 1e-200
#define PI  3.141592653589793 
#define LN0 5.0e22   //1.0/(4*PI*r0*r0*lnC)=1e24/20.  where  lnC=20.0 and r0=2.82e-13
#define LAMBDA 1.0e8 // lambda for turbulence

///////////////////////////////////////////////////////////////////////////////
// PROTOTYPES
///////////////////////////////////////////////////////////////////////////////

//void inj(double t, double t0, double dt, double tau, double d0, double dp, double p0, double dx, double x0, double s[], double e[], double u[], int ne, int nu, int ns, double *fptr);

void tinj(double t, double dt, double t0, double tau, double *fptr, double *pinj, int n);
void initInj(double *pinj, double d0, double dp, double p0, double dx, double x0, double s[], double e[], double u[], int ne, int nu, int ns);
void sup(double dt, double ds, double beta[], double u[], int ne, int nu, int ns, double *fptr);
void bup(double dt, double dbds[], double onmu22[], double dmu[], double beta[], double u[], int ne, int nu, int ns, double *fptr);
void eup(double dt, double delog, double emin, double beta2[],double np[], double e[], int ne, int nu, int ns, double *fptr);
void mup(double dt, double onmu2[], double rb3g2[], double np[], double dmu[], double u[], int ne, int nu, int ns, double *fptr);
void scup(double dt, double onmu2[], double beta[], double lambda, double dmu[], double u[], int ne, int nu, int ns, double *fptr);
void bfield(double s, double mr, double l, double n, double &b, double &dbds, int model);
double density(double s, double l, double nmin, double nmax);

void readBfield(double s[], double dbds[], int ns);
double linterpol(double x, double xm[], double ym[], int n);

void progressbar( int percent , double time, double tnow, double estim);

using namespace std;

/********************************************************************************
 * MAIN
 ********************************************************************************/

int main(int argc, char** argv)
{

  clock_t pStart = clock();

  //bool output = 1;

  string outputFile="fp.out";
  string inputFile="fp.in";

  int optind=1;
  bool dt_override=0;
  double dt_val=0.0;;
  /* decode arguments */
  while ((optind < argc) && (argv[optind][0]=='-')) 
    {
      string sw = argv[optind];
     
      /* use switch '-o filename' to name the output file */
      if (sw=="-o") {
	optind++;
	outputFile = argv[optind];
      }
      /* use switch '-i filename' to name the input file */
      if (sw=="-i") {
	optind++;
	inputFile = argv[optind];
      }
      /* use switch '-dt value' to override dt */
      if (sw=="-t") {
	optind++;
	dt_val = atof(argv[optind]);
	dt_override=1;
      }
    }

 /* handling input file: */
  ifstream fin;
  fin.open(inputFile.c_str());
  if (!fin.is_open())
    {
      cerr << "ERROR: Could not open " << inputFile << " for reading." << endl;
      fin.clear();
      fin.close();
      exit(EXIT_FAILURE);
    }
  
  /* handling output file: */
  ofstream fout;
  fout.open(outputFile.c_str(), ios::out | ios::binary);
  if (!fout.is_open())
    {
      cerr << "ERROR: Could not open " << outputFile << " for writing." <<endl;
      fout.clear();
      fout.close();
      exit(EXIT_FAILURE);
    }

  int ne,nu,ns,nt,narray;
  fin >> ne; // energy (30)
  fin >> nu; // pitch angle (43)
  fin >> ns; // position (101)
  fin >> nt; // time (output)
  narray=ne*nu*ns;

   // grid resolutions
  double dt,ds,du,delog;

  double * e      = new double [ne]; // energy in units of mc2
  double * beta   = new double [ne]; // beta = velocity / c 
  double * rb3g2  = new double [ne]; // beta^3 * gamma^2 (for diffusion coef.) 
  double * beta2  = new double [ne]; // beta at 1/2 grid step 
  double * u      = new double [nu]; // pitch-angle cosine (u=cos(phi)) 
  double * onmu2  = new double [nu]; // 1-u*u 
  double * onmu22 = new double [nu]; 
  double * dmu    = new double [nu]; // u grid step (in case du is not uniform)   
  double * s      = new double [ns]; // position 
  double * dbds   = new double [ns]; // d(lnB)/ds 
  double * np     = new double [ns]; // plasma density array 
  double * b      = new double [ns]; // magnetic field (just for output) 
  double * tr     = new double [nt]; // array of time points to write output
  double * f      = new double [narray]; // distribution function f(E,u,s)
  double * finj   = new double [narray]; // injection function S(E,u,s)
  double * fdata  = new double [narray*nt]; // storage

  ///////////////////////////////////////////////////////////////////////////////  

  //time_t start = time(NULL); // start time (evaluating performance)
  int bmodel;
  double emax,emin,smax,smin,mr,nexp,nmin,nmax,tmax,d0,dp,p0,dx,x0,t0,tau,lambda;
  bool run[6];

  ///////////////////////////////////////////////////////////////////////////////
  // INPUT PARAMETERS (read from file)
  ///////////////////////////////////////////////////////////////////////////////

  fin >> emin;
  fin >> emax;
  fin >> smin;
  fin >> smax;
  fin >> mr;
  fin >> nexp;
  fin >> bmodel;
  fin >> nmin;
  fin >> nmax;
  fin >> tmax;
  fin >> d0;
  fin >> dp;
  fin >> p0;
  fin >> dx;
  fin >> x0;
  fin >> t0;
  fin >> tau;
  fin >> lambda;
  fin >> run[0];
  fin >> run[1];
  fin >> run[2];
  fin >> run[3];
  fin >> run[4];
  fin >> run[5];

  /* closing input file */
  fin.clear();
  fin.close();

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  /*
    cout << lambda << endl;
    lambda = LAMBDA;
    cout << lambda << endl;
  */
 
  /*
   * Initializing: grid arrays
   */

  double gammhalf,gamm;
  //  double em1;
  delog = (log10(emax)-log10(emin))/float(ne-1);
 // create energy grid
  for (int i=0;i<ne;i++)
    {
      e[i]=pow(10.0,float(i)*delog+log10(emin))/E0;
      //em1=pow(10.0,float(i-1)*delog+log10(emin))/E0;
      //      de[i] = e[i] - em1;
      gamm = e[i]+1.0;
      beta[i] = sqrt(1.0-1.0/(gamm*gamm));
      gammhalf = pow(10.0,float(i-0.5)*delog+log10(emin))/E0 + 1.0;
      beta2[i]=sqrt(1.0-1.0/(gammhalf*gammhalf));
      rb3g2[i] = 1.0/(pow(beta[i],3.0)*gamm*gamm);
    }

  // create pitch-angle grid
  du = 2.0/float(nu-1); // (+1.0 - (-1.0))/(nu-1)
  for (int i=1;i<nu;i++)
    {
      u[i]=float(i)*du-1.0;
      dmu[i] = du; // cte, change later if needed
      onmu2[i] = 1.0 - u[i]*u[i];
      onmu22[i] = 1.0 - 0.25*(u[i]+u[i-1])*(u[i]+u[i-1]);
      // onmu22 is (1-mu^2) for mu[i+1/2]
    }
  u[0]=-1.0;
  dmu[0]=du;
  onmu2[0]=0.0; // 1.0-(-1)^2=0;

  // create position grid 
  // create magnetic field array
  double sj,nj;
  double ib,idbds;       // temp variables for the output of bfield()  
  ds=(smax-smin)/float(ns-1);
  for (int i=0;i<ns;i++)
    { 
      sj = i*ds+smin;
      s[i]=sj;
      nj = density(sj,smax,nmin,nmax);
      np[i] = nj;
      //  if (s[i] < 0) mr=2.0; else mr=4.0;
      bfield(sj,mr,smax,nexp,ib,idbds,bmodel);
      b[i] = ib;
      dbds[i] = idbds;//*LN0/nj;
    }
  if (bmodel == 2)
    {
      cout << "BFIELD: INPUT BY USER" << endl;
      readBfield(s,dbds,ns);
    }

  /********************************************************************************** 
  // DEFINITION OF TIME STEP
  ***********************************************************************************/
  //  dycoul = 0.9*min(0.5*du/rb3g2[0], beta[0]*de[1]);
  //double dycoul = beta[0]*de[0];
  double dycoul = beta[0]*(emin-pow(10.0,(-1.0)*delog+log10(emin)))/E0;

  //  double dtcoul = beta[0]*de[0] / nmax * LN0/C;
  double dydiff = 0.5*dmu[nu-1]/rb3g2[0];
  dycoul = min(dydiff,dycoul);
  double dtcoul = dycoul /C/nmax*LN0; //1e24/20.;

  // advective term
  //  v=ds/dt=mu*beta*c -> dt=ds/(mu*beta*c)
  double dtadv =  ds/(beta[ne-1]*C); // mu=1 max(beta) (worst case)

  // mirroring term
  double dymir = 1e20;
  for (int i=1;i<nu-1;i++)
    dymir = min(dymir,du/onmu2[i]);
  double dymir0=dymir;
  for (int i=0;i<ns;i++)
    dymir = min(dymir,dymir0*nmax/LN0/abs(dbds[i]));
  // dymir = dymir/(0.5*beta[ne-1]*C);
  dymir = dymir/(0.5*beta[ne-1]);
  //  double dtmir = LN0/nmax/C*dymir;
  double dtmir = LN0/nmax/C*dymir;

  // dt for turbulence scattering (need to confirm the relativistic correction)
  double dtturb = 2.0*lambda*dmu[nu-1]/C/beta[ne-1]/pow(sqrt(1.0/(1-beta[ne-1]*beta[ne-1])),2./3.);

  // find minimum dt (all the 'ifs' are there so we just use dt for the processes we asked for.)
  dt = 0.1; // whatever value to start things
  if (run[0]) dt = min(dt,dtadv);
  if (run[2]) dt = min(dt,dtmir);
  if (run[1] or run[3]) dt = min(dt,dtcoul);
  if (run[4]) dt = min(dt,dtturb);
 
  cout << "TIME SCALES\n";
  cout << "dt[adv]  = " << dtadv << endl;
  cout << "dt[mir]  = " << dtmir << endl;
  cout << "dt[coul] = " << dtcoul << endl;
  cout << "dt[turb] = " << dtturb << endl;
  cout << "dt[min]  = " << dt    << endl;

  if (dt_override)
    {
      dt=dt_val;
      cout << "USER: dt[override]  = " << dt    << endl;
    }

  // INITIALIZING TIME VARIABLES
  double t=0.0;       // time [sec]
  int itr=0,tstep=0;  // time step (for integration and output) counters

  double dtw=tmax/float(nt-1);
  dtw=max(dt,dtw);
  for (int i=0;i<nt;i++)
    tr[i]=float(i)*dtw;
  cout << "dt[write]  = " << dtw    << endl;
///////////////////////////////////////////////////////////////////////////////
// INIT: DONE.
///////////////////////////////////////////////////////////////////////////////
// MAIN
///////////////////////////////////////////////////////////////////////////////
 
// WRITING

  fout.write(reinterpret_cast<const char*> (&ne), sizeof(int));
  fout.write(reinterpret_cast<const char*> (&nu), sizeof(int));
  fout.write(reinterpret_cast<const char*> (&ns), sizeof(int));
  fout.write(reinterpret_cast<const char*> (&nt), sizeof(int));
  fout.write(reinterpret_cast<const char*> (&delog), sizeof(double));
  fout.write(reinterpret_cast<const char*> (&du), sizeof(double));
  fout.write(reinterpret_cast<const char*> (&ds), sizeof(double));
  fout.write(reinterpret_cast<const char*> (&dt), sizeof(double));
  fout.write(reinterpret_cast<const char*> (e), sizeof(double)*ne);
  fout.write(reinterpret_cast<const char*> (u), sizeof(double)*nu);
  fout.write(reinterpret_cast<const char*> (s), sizeof(double)*ns);
  fout.write(reinterpret_cast<const char*> (b), sizeof(double)*ns);
  fout.write(reinterpret_cast<const char*> (dbds), sizeof(double)*ns);
  fout.write(reinterpret_cast<const char*> (np), sizeof(double)*ns);
  
  /*
  // WRITE GRID SIZES
  fout << ne << "\t" << nu << "\t" << ns<<"\t" << nt << endl;
  // WRITE GRID RESOLUTIONS
  fout << delog << "\t" << du << "\t" << ds << "\t" << dt << endl;
  
  for (int i=0;i<ne;i++)
    fout << e[i]*E0 << "\t";   // WRITE ENERGY ARRAY
  fout << endl;
  
  for (int i=0;i<nu;i++)
    fout << u[i] << "\t";      // WRITE PITCH-ANGLE ARRAY
  fout << endl;
  
  for (int i=0;i<ns;i++)
    fout << s[i] << "\t";      // WRITE POSITION ARRAY
  fout << endl;
  
  for (int i=0;i<ns;i++)
    fout << b[i] << "\t";      // WRITE b(s) ARRAY
  fout << endl;
  
  for (int i=0;i<ns;i++)
    fout << dbds[i] << "\t";   // WRITE dbds(s) ARRAY
  fout << endl;
  
  for (int i=0;i<ns;i++)
    fout << np[i] << "\t";     // WRITE DENSITY ARRAY
  fout << endl;
  */
///////////////////////////////////////////////////////////////////////////////
// MAIN LOOP
///////////////////////////////////////////////////////////////////////////////

  double *ptr = &f[0];
  double *pinj = &finj[0];

  // initializing f
  //for (int i=0;i<(narray);i++)
  //  f[i]=0.0;

  // computes the injection function finj(E,u,s)
  initInj(pinj,d0,dp,p0,dx,x0,s,e,u,ne,nu,ns);

  // computes the initial injection into distribution
  tinj(t,dt,t0,tau,ptr,pinj,narray);

  //inj(t,t0,dt,tau,d0,dp,p0,dx,x0,s,e,u, ne,nu,ns, ptr);

  float percent;
  int ttotal = tmax/dt;
  clock_t tStart = clock();
  double tNow=0.0,tMean=0.0,tEstim=0.0;

  while (t < tmax)
    {
      // progress bar  //////////////////////////////
      tNow = (double)(clock() - tStart)/CLOCKS_PER_SEC;
      percent = floor(float(tstep)/ttotal*100.0);
      progressbar(percent,tEstim,tNow,tMean*ttotal);
      tstep++; // increment time step counter  
      tNow = (double)(clock()-tStart)/CLOCKS_PER_SEC;
      tMean = (tMean+tNow/(double)(tstep))/2.0;
      tEstim = tMean * (tmax/dt-tstep);
      //////////////////////////////////////////////
      
      // if t steps over tr for output, write output
      if (t >= tr[itr]) 
	{
	  // cout << t << "  "<< tr[itr] << "  "<< t-tr[itr] <<endl;
	  // write output
	  /*
	    fout << t << endl;
	    for (int i=0;i<(narray);i++)
	    fout << f[i] << "\n";
	  */
	  tr[itr]=t;
	  for (int i=0;i<narray;i++)
	    fdata[i+narray*itr]=f[i];

	  //cout << narray*itr << "\t" << (narray-1)+narray*itr << endl;

	  //fout.write(reinterpret_cast<const char*> (f), sizeof(double)*narray);

	  itr++; // advances output time step
	}

      // advective term
      if (run[0]) 
	sup(dt,ds,beta,u, ne,nu,ns, ptr);
      
      // collisional loss term
      if (run[1]) 
	eup(dt,delog,emin,beta2,np,e, ne,nu,ns, ptr);
      
      // mirroring term
      if (run[2]) 
	bup(dt, dbds, onmu22, dmu, beta, u, ne,nu,ns, ptr);
	   
      // collisional diffusion term
      if (run[3]) 
	mup(dt, onmu2, rb3g2, np, dmu, u, ne,nu,ns, ptr);

      // general diffusion term
      if (run[4]) 
      scup(dt, onmu2, beta, lambda, dmu, u, ne,nu,ns, ptr);

      // injection term
      if (run[5])
	tinj(t,dt,t0,tau,ptr,pinj,narray);
     
      t += dt; // step time
    }
  /* END OF MAIN LOOP */

  fout.write(reinterpret_cast<const char*> (&itr), sizeof(int));
  fout.write(reinterpret_cast<const char*> (tr) , sizeof(double)*itr);
  fout.write(reinterpret_cast<const char*> (fdata), sizeof(double)*narray*itr);

  /* closing output file */
  fout.clear();
  fout.close();

  delete [] fdata;
  delete [] f;
  delete [] finj;
  delete [] e; 
  delete [] beta;
  delete [] rb3g2;
  delete [] beta2;
  delete [] u;
  delete [] onmu2;
  delete [] onmu22;
  delete [] dmu;
  delete [] s;
  delete [] dbds;
  delete [] np;
  delete [] b;
  delete [] tr;

  cout << "\nprogram done.\n";
  cout << "elapsed time: " << (double)(clock() - pStart)/CLOCKS_PER_SEC <<endl;
//cout << "elapsed time (sec): " << difftime(time(NULL),start) << endl;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// SUP
///////////////////////////////////////////////////////////////////////////////

void sup(double dt, double ds, double beta[], double u[], int ne, int nu, int ns, double *fptr)
{
 
  /*
    variables for 
    f_m^n = phic (center)
    f_{m+1}^n = phin (next)
    f_{m-1}^n = phip (previous)
    delt = (f_m^n-f_{m-1}^n)*(f_{m+1}^n-f_m^n)
    a1ma = a(1-a)
    a1pa = a(1+a) // note: typo in the paper: it is "a(1+a)", not "a(1-a)" for a<0.
    dphi, dphip are the 
  */

  double a0 = C*dt/ds, a01;

  double tfin[nu][ns];  
  /* temporary array for input distribution 
     f^n (n = time, m = position)*/
  double tfout[nu][ns]; 
  /* temporary array for output distribution f_m^{n+1}*/

  double a,a1pa,phic,dphi,delt,si,dphip,deltp,phip,a1ma,phin,delto;
  int index;

  //#pragma omp parallel for schedule(dynamic,4) private(tfin,tfout,a,a01,a1pa,phic,dphi,delt,si,dphip,deltp,phip,a1ma,phin,delto,index)

  for (int k=0;k<ne;k++)
    {
      a01 = a0 * beta[k];
          
      for (int l=0;l<nu;l++)
	for (int m=0;m<ns;m++)
	  {
	    index = k*nu*ns + l*ns + m;	
	    tfin[l][m] = fptr[index];
	  }
    
      for (int l=0;l<nu;l++)
	{
	  a = a01*u[l];  /* electron longitudinal speed ('a' 
			    in Hamilton et al, 1990)    */
	  
	  if (a < 0)       /* NEGATIVE VELOCITIES (Eq. 8 
			      in Hamilton et al, 1990, for a < 0)*/
	    {
	      a1pa = a*(1.0+a);
	    
	      phic = tfin[l][ns-1];
	      phin = tfin[l][ns-2];

	      dphi = phic - phin;   
	      delt = -dphi*phic;
	      si = copysign(1.0,delt);
	      delt = -delt/(phin+EPS)*(si+1.0)/2.0;

	      tfout[l][ns-1] = phic*(1.0+a)-a1pa*delt*(si+1.0)/2.;
	      
	      for (int m=ns-2;m>0;m--)
		{
		  phip = phic; 
		  phic = phin;
		  phin = tfin[l][m-1];
		  
		  dphip = dphi;
		  dphi = phic - phin;
		  deltp = delt;
		  delt = dphi*dphip;
		  si = copysign(1.0,delt);
		  delt=delt/(phip - phin+EPS)*(si+1.0)/2.;
		  tfout[l][m] = phic - a*dphip + a1pa*(deltp - delt);
		}
	      phip = phic;
	      phic = phin;
	      dphip = dphi;
	      deltp = delt;
	      delt = tfin[l][0]*dphip;
	      si=copysign(1.,delt);
	      delt= delt/(tfin[l][1]+EPS)*(si+1.)/2.;     
	      tfout[l][0] = phic - a*dphip + a1pa*(deltp - delt);
	    }
	  else 	      // POSITIVE VELOCITIES (Eq. 8 in Hamilton et al, 1990, for a > 0)
	    {      
	      a1ma = a*(1.0-a);
	      phic = tfin[l][0];
	      phin = tfin[l][1];
	      dphip = phin - phic;
	      delt = phic*dphip;
	      si=copysign(1.,delt);
	      delt= delt/(phin+EPS)*(si+1.0)/2.;
	      tfout[l][0] = phic*(1.0-a) - a1ma*delt;
	      for (int m=1;m<ns-1;m++)
		{
		  phip = phic;
		  phic = phin;
		  phin = tfin[l][m+1];
		  delto = delt;
		  dphi = dphip;
		  dphip = phin - phic;
		  delt = dphi*dphip;
		  si=copysign(1.0,delt);
		  delt=delt/(phin - phip+EPS)*(si+1.0)/2.;
		  tfout[l][m] = phic - a*dphi - a1ma*(delt - delto);
		}
	      phip = phic;
	      phic = phin;
	      delto = delt;
	      dphi = dphip;
	      delt = -phic*dphi;
	      si=copysign(1.0,delt);
	      delt= -delt/(phip+EPS)*(si+1.0)/2.;
	      tfout[l][ns-1]=phic- a*dphi - a1ma*(delt-delto);
	    }
	
	  for (int l=0;l<nu;l++)
	    for (int m=0;m<ns;m++)
	      {
		index = k*nu*ns + l*ns + m; 
		fptr[index] = tfout[l][m];
		//f[k][l][m] = tfout[l][m];
	      }
	}
    }
}
///////////////////////////////////////////////////////////////////////////////

void bup(double dt, double dbds[], double onmu22[], double dmu[], double beta[], double u[], int ne, int nu, int ns, double *fptr)
{

  double tfin[nu][ns];  // temporary array for input distribution f_m^n
  double tfout[nu][ns]; // temporary array for output distribution f_m^{n+1}
  
  double a,aa,psf,fl,flp,d3,d2,delt,d1,si,dk;
  int lm1,lp1,index;

  //#pragma omp parallel for schedule(dynamic,4) private(a,aa,psf,fl,flp,d3,d2,delt,d1,si,dk,lm1,lp1,tfin,tfout,index)
  for (int k=0;k<ne;k++)
    {  
      a = dt*C*beta[k];
      for (int l=0;l<nu;l++)
	for (int m=0;m<ns;m++)
	  {
	    index = k*nu*ns + l*ns + m;	
	    tfin[l][m] = fptr[index];
	  }   
      for (int m=0;m<ns;m++) 
	{
	  aa = a*dbds[m];
	  if (aa == 0.0)
	    {
	      for (int l=0;l<nu;l++)
		tfout[l][m] = tfin[l][m];
	    }
	  else
	    {
	      psf = exp(-aa);
	      tfout[nu-1][m] = psf*tfin[nu-1][m];
	      tfout[0][m] = tfin[0][m]/psf; 
	      
	      if (aa > 0.0) // converging field. Eqs. for V<0 (IMPORTANT)
		{		  
		  // fl = tfin[nu-1][m]*onmu22[nu-1]; //line below subst this
		  fl =  tfin[nu-1][m]*(1.0-0.25*(u[nu-1]+u[nu-2])*(u[nu-1]+u[nu-2]));
		  for (int l=nu-2;l>0;l--)
		    {
		      lm1 = l-1;
		      lp1 = l+1;
		      flp = fl;
		      d3 = tfin[l][m];
		      d2 = 0.5*(d3+tfin[lm1][m]);
		      delt = 0.5*dmu[l]/dmu[lp1];
		      d1 = d3*(1.+delt)-tfin[lp1][m]*delt;
		      si=copysign(1.,d3-tfin[lm1][m]);
		      dk=-max(d3,min(d1,d2))*(si-1.)/2.
			+min(d3,max(d1,d2))*(si+1.)/2.;
		      // fl = dk*onmu22[l];  //line below subst this
		      fl = dk*(1.0-0.25*(u[l]+u[lm1])*(u[l]+u[lm1]));
		      tfout[l][m]=d3+aa*(flp-fl)/(u[lp1]-u[lm1]); 
		      // aa sign is switched multiplying by -1, that's why d3+aa
		    }
		  
		}
	      else // also eqs for V<0
		{
		  //fl = tfin[0][m]*onmu22[1]; //line below subst this
		  fl =  tfin[0][m]*(1.0-0.25*(u[1]+u[0])*(u[1]+u[0]));
		  for (int l=1;l<nu-1;l++) 
		    {
		      lm1 = l-1;
		      lp1 = l+1;
		      flp = fl;
		      d3 = tfin[l][m];
		      d2 = 0.5*(tfin[lp1][m]+d3);
		      delt = 0.5*dmu[lp1]/dmu[l];
		      d1 = d3*(1.+delt)-tfin[lm1][m]*delt;
		      si=copysign(1.,d3-tfin[lp1][m]); 
		      dk=-max(d3,min(d1,d2))*(si-1.)/2.
			+min(d3,max(d1,d2))*(si+1.)/2.;
		      //fl = dk*onmu22[lp1]; //line below subst this
		      fl = dk*(1.0-0.25*(u[lp1]+u[l])*(u[lp1]+u[l])) ;
		      tfout[l][m]=d3+aa*(fl-flp)/(u[lp1]-u[lm1]);
		    }
		}
	    }
	}
      
      for (int l=0;l<nu;l++)
	for (int m=0;m<ns;m++)
	  {
	    index = k*nu*ns + l*ns + m;	
	    fptr[index] = tfout[l][m];
	  }
    }
  
}

///////////////////////////////////////////////////////////////////////////////

void eup(double dt, double delog, double emin,double beta2[], double np[], double e[], int ne, int nu, int ns, double *fptr)
{
 
  //  double tin[ne][nu][ns];
  //  double tout[ne][nu][ns];
  //  double * tin = new double [ne*nu*ns];

  //double * tout = new double [ne*nu*ns];
  double tout[ne*nu*ns];
  
  double a,d3,d2,d1,dk,si,si0,f1,f1p,delt;
  int km1,kp1,index,indexm1,indexp1; 

  // temporary array
  /*  for (int k=0;k<ne;k++)
    for (int l=0;l<nu;l++) 
      for (int m=0;m<ns;m++)   
	{
	  index = k*nu*ns + l*ns + m;	
	  //tin[k][l][m] = fptr[index];
	  tin[index] = fptr[index];
	}
  */
  double acte = 2.0*dt*C/LN0;

  //#pragma omp parallel for schedule(dynamic,4) private(a,d3,d2,d1,dk,si,si0,f1,f1p,delt,km1,kp1,indexm1,indexp1) shared(fptr,tout,acte)
  for (int m=0;m<ns;m++)   
    for (int l=0;l<nu;l++) 
      {   
	a=acte*np[m];
	//k = ne-1;
	km1 = ne-2;
	d3 = fptr[(ne-1)*nu*ns+l*ns+m];
	d1 = fptr[(ne-2)*nu*ns+l*ns+m];
	//d3 = tin[ne-1][l][m];
	//d1 = tin[km1][l][m]; 
	/* here d1 is not the same as d1 inside 
		    the 'for' loop (just used d1 for U_{m-1}) */
	d2 = 0.5*(d3+d1);
	si0 = copysign(1.0,d3-d1);
	dk = max(d3,d2)*(1.0-si0)/2.0+min(d3,d2)*(si0+1.)/2.0;
	f1 = dk/beta2[ne-1];
	//tout[ne-1][l][m] = d3 - a*f1/(pow(10.0,delog*ne+log10(emin))/E0 - e[km1]);
	tout[(ne-1)*ns*nu+l*ns+m] = d3 - a*f1/(pow(10.0,delog*ne+log10(emin))/E0 - e[km1]);
	  
	/* pow(10.0,delog*ne+log10(emin)) gives exactly the 'next' 
	   grid point e[ne], as it would exist.
	   this works because the grid is regular in log10(e),
	   with step 'delog'.*/
	for (int k=ne-2;k>=1;k--)
	  {	
	    kp1 = k+1;
	    km1 = k-1;
	    index=k*nu*ns+l*ns+m;
	    indexm1=km1*nu*ns+l*ns+m;
	    indexp1=kp1*nu*ns+l*ns+m;

	    f1p = f1;
	    //d3 = tin[k][l][m];
	    d3 = fptr[index];
	    //d2 = 0.5*(d3 + tin[km1][l][m]);
	    d2 = 0.5*(d3 + fptr[indexm1]);
	    delt=0.25*(e[k+1]-e[k-1])/(e[k+1]-e[k]);
	    //d1 = d3*(1.0+delt)-tin[kp1][l][m]*delt;
	    d1 = d3*(1.0+delt)-fptr[indexp1]*delt;
	    //si = copysign(1.0,d3-tin[km1][l][m]);
	    si = copysign(1.0,d3-fptr[indexm1]);
	    dk = max(d3,min(d1,d2))*(1.0-si)/2.+min(d3,max(d1,d2))*(si+1.)/2.;
	    f1 = dk/beta2[k]; /* beta2[k] stands for V_{m-1/2} 
				 (the constant part is a)*/
	    //tout[k][l][m] = d3 + a*(f1p-f1)/(e[kp1]-e[km1]);
	    tout[index] = d3 + a*(f1p-f1)/(e[kp1]-e[km1]);
	  }
	/*k = 0;*/
	kp1 = +1;
	km1 = -1;
	f1p = f1;
	//d3 = tin[0][l][m];
	d3 = fptr[l*ns+m];
	delt=0.25*(e[kp1]-pow(10.0,delog*(-1.0)+log10(emin))/E0)/(e[kp1]-e[0]);
	//d1 = d3*(1.0+delt)-tin[kp1][l][m]*delt;
	d1 = d3*(1.0+delt)-fptr[kp1*nu*ns+l*ns+m]*delt;
	dk = min(d3,d1);
	f1 = dk/beta2[0];
	//tout[0][l][m] = d3 + a*(f1p-f1)/(e[kp1]-pow(10.0,delog*(-1.0)+log10(emin))/E0);
	tout[l*ns+m] = d3 + a*(f1p-f1)/(e[kp1]-pow(10.0,delog*(-1.0)+log10(emin))/E0);
      }

  for (int k=0;k<ne;k++)
    for (int l=0;l<nu;l++) 
      for (int m=0;m<ns;m++) 
	{  
	  index = k*nu*ns + l*ns + m;
	  //fptr[index] = tout[k][l][m];
	  fptr[index] = tout[index];
	}
  
  //  delete [] tin;
  //delete [] tout;  
}  

///////////////////////////////////////////////////////////////////////////////
void mup(double dt, double onmu2[], double rb3g2[], double np[], double dmu[], double u[], int ne, int nu, int ns, double *fptr)
{
  /* f[k][l][m] is updated for the pitch angle 
     diffusion term.  The Crank-Nicholson method is employed and 
     subroutine cntridag is called to solve the tridiagonal matrix.
     (from original code)  */

  int lp1,lm1,index,indexp1,indexm1,index0,indexn;
  double dy,b4,dydsum,b5,aa,bb1,bet;

  double * a = new double[ne*nu*ns];
  double * b = new double[ne*nu*ns];
  double * c = new double[ne*nu*ns];
  double * r = new double[ne*nu*ns];
  double * gam = new double[ne*nu*ns];
  double * tin = new double[ne*nu*ns];
  double * tout = new double[ne*nu*ns];

  /*
  double a[ne*nu*ns];
  double b[ne*nu*ns];
  double c[ne*nu*ns];
  double r[ne*nu*ns];
  double gam[ne*nu*ns];
  double tin[ne*nu*ns];
  double tout[ne*nu*ns];
  */

  /* for (int k=0;k<ne;k++)
    for (int l=0;l<nu;l++) 
      for (int m=0;m<ns;m++)   
	{
	  index = k*nu*ns + l*ns + m;	
	  tin[index] = fptr[index];
	  }*/

  for (int i=0;i<(ne*nu*ns);i++)
    tin[i] = fptr[i];
	
  dy =  C/LN0 * dt;
  

  //#pragma omp parallel shared(a,b,c,r,gam,tout)
  //#pragma omp for schedule(dynamic,4) private(lp1,lm1,index,indexp1,indexm1,dydsum,b4,b5,aa,bb1)
  for (int k=0;k<ne;k++)
    for (int m=0;m<ns;m++)
      for (int l=1;l<nu-1;l++)
	{  
	  b4 = np[m]*rb3g2[k];
	  lp1 = l+1;
	  lm1 = l-1;

	  index = k*nu*ns + l*ns + m;	
	  indexp1 = k*nu*ns + lp1*ns + m;	
	  indexm1 = k*nu*ns + lm1*ns + m;	

	  dydsum = dy/(dmu[lp1]+dmu[l]);
	  b5 = b4*onmu2[l];
	  aa = b5*dydsum;
	  bb1 = (-b4*u[l])*dydsum;
	  /* 
	     a, b, c are the elements of the tridiagonal matrix. 
	     r is the vector product of the matrix and the vector 
	     of updated phi.  (note from original code) 
	  */
	  a[index] = bb1-aa/dmu[l];
	  b[index] = 1.0+b5*dy/(dmu[l]*dmu[lp1]);
	  c[index] = -bb1-aa/dmu[lp1];
	  r[index] = aa*((tin[indexp1]-tin[index])/dmu[lp1]-
			   (tin[index]-tin[indexm1])/dmu[l])+bb1*(tin[indexp1]-tin[indexm1])+tin[index];
	}
  //#pragma omp parallel shared(a,b,c,r,gam,tout)
  //#pragma omp for schedule(dynamic,4) private(index,indexp1,indexm1,b4,index0,indexn,bet)
  for (int k=0;k<ne;k++) 
    for (int m=0;m<ns;m++) 
      {
	  //index = k*nu*ns + l*ns + m;	
	  index0 = k*nu*ns + 0*ns + m;	
	  indexn = k*nu*ns + (nu-1)*ns + m;
	  b4 = np[m]*rb3g2[k];
	  /* at amu=-1, the only term remaining is the first derivative 
	     from the pitch angle diffusion term which we approximate by 
	     a forward difference (note from original code) */
	  c[index0] = -b4*dy/dmu[0];
	  b[index0] = 1.0-c[index0];
	  r[index0] = -c[index0]*(tin[k*nu*ns + 1*ns + m]-tin[index0])+tin[index0];
	  /* at amu=1, use the backward difference  (note from original code) */
	  a[indexn] = -b4*dy/dmu[nu-1];
	  b[indexn] = 1.0-a[indexn];
	  r[indexn] = a[indexn]*(tin[indexn]-tin[k*nu*ns+(nu-2)*ns+m])+tin[indexn];
	  /* solving tridiagonal matrix */
	  bet = b[index0];
	  tout[index0] = r[index0]/bet;
	  for (int j=1;j<nu;j++)
	    {
	      index = k*nu*ns + j*ns + m;	
	      indexm1 = k*nu*ns + (j-1)*ns + m;	
	      gam[index] = c[indexm1]/bet;
	      bet = b[index]-a[index]*gam[index];
	      tout[index] = (r[index]-a[index]*tout[indexm1])/bet;
	    }
	  for (int j=nu-2;j>=0;j--)
	    {
	      index = k*nu*ns + j*ns + m;
	      indexp1 = k*nu*ns + (j+1)*ns + m;
	      tout[index] = tout[index]-gam[indexp1]*tout[indexp1];
	    }
      }
  
    for (int i=0;i<(ne*nu*ns);i++)
      fptr[i] = tout[i];

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] r;
    delete [] gam;
    delete [] tin;
    delete [] tout;

    /*    for (int k=0;k<ne;k++)
      for (int l=0;l<nu;l++) 
	for (int m=0;m<ns;m++)   
	  {
	    index = k*nu*ns + l*ns + m;
	    fptr[index] = tout[index]; 
	    // f[k][l][m] = tout[k][l][m];
	    }*/
    
}

///////////////////////////////////////////////////////////////////////////////
void scup(double dt, double onmu2[], double beta[], double lambda, double dmu[], double u[], int ne, int nu, int ns, double *fptr)
{
  /* f[k][l][m] is updated for the pitch angle 
     diffusion term.  The Crank-Nicholson method is employed and 
     subroutine cntridag is called to solve the tridiagonal matrix.
     (from original code)  

     new subroutine: general pitch-angle scattering, controlled by lambda=mean free path for scattering

*/

  int lp1,lm1,index,indexp1,indexm1,index0,indexn;
  double dy,b4,dydsum,b5,aa,bb1,bet;

  double * a = new double[ne*nu*ns];
  double * b = new double[ne*nu*ns];
  double * c = new double[ne*nu*ns];
  double * r = new double[ne*nu*ns];
  double * gam = new double[ne*nu*ns];
  double * tin = new double[ne*nu*ns];
  double * tout = new double[ne*nu*ns];

  for (int i=0;i<(ne*nu*ns);i++)
    tin[i] = fptr[i];
	
  //dy =  C/LN0 * dt;
  dy = C*dt/2./lambda;
  
  //#pragma omp parallel shared(a,b,c,r,gam,tout)
  //#pragma omp for schedule(dynamic,4) private(lp1,lm1,index,indexp1,indexm1,dydsum,b4,b5,aa,bb1)
  for (int k=0;k<ne;k++)
    for (int m=0;m<ns;m++)
      for (int l=1;l<nu-1;l++)
	{  
	  b4 = beta[k]/pow(sqrt(1.0/(1-beta[k]*beta[k])),2./3.);// beta[k]*C/2.0/lambda //np[m]*rb3g2[k]; // speed (space and energy dependent)
	  lp1 = l+1;
	  lm1 = l-1;

	  index = k*nu*ns + l*ns + m;	
	  indexp1 = k*nu*ns + lp1*ns + m;	
	  indexm1 = k*nu*ns + lm1*ns + m;	

	  dydsum = dy/(dmu[lp1]+dmu[l]);

	  b5 = b4*onmu2[l];
	  
	  aa = b5*dydsum;
	  bb1 = (-b4*u[l])*dydsum;
	  /* 
	     a, b, c are the elements of the tridiagonal matrix. 
	     r is the vector product of the matrix and the vector 
	     of updated phi.  (note from original code) 
	  */
	  a[index] = bb1-aa/dmu[l];
	  b[index] = 1.0+b5*dy/(dmu[l]*dmu[lp1]);
	  c[index] = -bb1-aa/dmu[lp1];
	  r[index] = aa*((tin[indexp1]-tin[index])/dmu[lp1]-
			   (tin[index]-tin[indexm1])/dmu[l])+bb1*(tin[indexp1]-tin[indexm1])+tin[index];
	}
  //#pragma omp parallel shared(a,b,c,r,gam,tout)
  //#pragma omp for schedule(dynamic,4) private(index,indexp1,indexm1,b4,index0,indexn,bet)
  for (int k=0;k<ne;k++) 
    for (int m=0;m<ns;m++) 
      {
	  //index = k*nu*ns + l*ns + m;	
	  index0 = k*nu*ns + 0*ns + m;	
	  indexn = k*nu*ns + (nu-1)*ns + m;

	  	  b4 = beta[k]/pow(sqrt(1.0/(1-beta[k]*beta[k])),2./3.);//*C/2.0/lambda; //speed (space and energy dependent)

	  /* at amu=-1, the only term remaining is the first derivative 
	     from the pitch angle diffusion term which we approximate by 
	     a forward difference (note from original code) */
	  c[index0] = -b4*dy/dmu[0];
	  b[index0] = 1.0-c[index0];
	  r[index0] = -c[index0]*(tin[k*nu*ns + 1*ns + m]-tin[index0])+tin[index0];
	  /* at amu=1, use the backward difference  (note from original code) */
	  a[indexn] = -b4*dy/dmu[nu-1];
	  b[indexn] = 1.0-a[indexn];
	  r[indexn] = a[indexn]*(tin[indexn]-tin[k*nu*ns+(nu-2)*ns+m])+tin[indexn];
	  /* solving tridiagonal matrix */
	  bet = b[index0];
	  tout[index0] = r[index0]/bet;
	  for (int j=1;j<nu;j++)
	    {
	      index = k*nu*ns + j*ns + m;	
	      indexm1 = k*nu*ns + (j-1)*ns + m;	
	      gam[index] = c[indexm1]/bet;
	      bet = b[index]-a[index]*gam[index];
	      tout[index] = (r[index]-a[index]*tout[indexm1])/bet;
	    }
	  for (int j=nu-2;j>=0;j--)
	    {
	      index = k*nu*ns + j*ns + m;
	      indexp1 = k*nu*ns + (j+1)*ns + m;
	      tout[index] = tout[index]-gam[indexp1]*tout[indexp1];
	    }
      }
  
    for (int i=0;i<(ne*nu*ns);i++)
      fptr[i] = tout[i];

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] r;
    delete [] gam;
    delete [] tin;
    delete [] tout;
    
}
///////////////////////////////////////////////////////////////////////////////

double density(double s, double l, double nmin, double nmax)
{
  // my model (PJ)
  //double k=100;
  //return nmin*exp( pow(s/l,k) * log(nmax/nmin));

  /* exponential model (Burge et al 2012) FIXED for loop with
     L=[-8.32e8,8.32e8], to get density=1e12 at each footpoint
   */
  //  double nc=5e10;
  //  double np=1e14;
  double nc=nmin;
  double np=nmax;
  double h0=3e7; //3e8
  if (nmin == nmax) 
    {return nmax;}
  else 
    {return nc+np*exp(-(l-abs(s))/h0);}
}

///////////////////////////////////////////////////////////////////////////////

void bfield(double si, double mr, double l, double n, double &b, double &dbds, int model)
{
  /* **************************************************************
     returns the magnetic field b and its log derivative d(ln B)/ds
     for parabolic model (1) or exponential model (0)
     ************************************************************** */

  float b0 = 1.0;

  switch (model)
    {   
    case 2:  /* read bfield from file */
      {
	b=0.0;
	dbds=0.0;
	break;
      }
    case 1:
      {
	/* **********************************
	   Parabolic model: 
	   B(s) = B0 * (1+s^n/LB^n)  
	   where: LB = l / sqrt(mr-1.0)
	   l: half loop length (cm)
	   B0: magnetic field at looptop B(l/2)=B0
	   mr: mirror ratio B(s)/B0
	   n: n=2 (parabolic), higher value can be used (n=4)
	   ************************************* */
	double LB = l / sqrt(mr-1.0);
	b = b0 * (1.0 + pow(si,n) / pow(LB,n));
	dbds = n * pow(si,n-1.0) / (pow(LB,n) + pow(si,n)); //dbds = n*(s)**(n-1.0) / (LB**n + (s)**n)
	break;
      }
    case 0:
      {
	/* **********************************
	   Exponential model: 
	   B(s) = B0 * exp(log(mr)*(s/l)^n)
	   where:
	   L: half loop length (cm)
	   B0: magnetic field at looptop B(L/2)=B0
	   mr: mirror ratio B(s)/B0
	   n: higher values cause stronger convergence near the footpoints
	   ************************************* */
	dbds = log(mr) / pow(l,n) * n * pow(si,n-1.0); // d (ln B) / ds
	b = b0 * exp(log(mr) * pow(si/l,n));           // b(s)
	break;
      }
    }
}

///////////////////////////////////////////////////////////////////////////////

void initInj(double *pinj, double d0, double dp, double p0, double dx, double x0, double s[], double e[], double u[], int ne, int nu, int ns)
{

  int index;
  double ae,au,as,x,p;
  double fe[ne],fu[nu],fs[ns];

  /* normalization */
  ae = (d0-1.0)/(pow(e[0]*E0,(1.-d0))-pow(e[ne-1]*E0,(1.-d0)));
  // au = 1.0/(sqrt(2.0*PI)*dp);
  // as = 1.0/(sqrt(2.0*PI)*dx);

  for (int i=0;i<ne;i++)
    if (e[i]*E0 >= 10.) 
      fe[i] = ae*pow(e[i]*E0,-d0);
    else 
      fe[i]=0.0;

  // normalization of fu
  int nint=100;
  double iu[nint],wu[nint],fiu=0.0;
  gauleg(-1.0,1.0,iu,wu,nint);
  for (int i=0;i<nint;i++)
    {
      if (p0 == 2) 
	{
	  fiu+=wu[i]*(exp(-0.5*((iu[i]-1.0)*(iu[i]-1.0))/(dp*dp)));
	}
      else 
	{
	  p=iu[i]-p0;
	  fiu+=wu[i]*exp(-0.5*(p*p)/(dp*dp));
	}
    }
  au=1.0/fiu;
  // done.

  // normalization of fs

  gauleg(s[0],s[ns-1],iu,wu,nint);
  fiu=0.0;
  for (int i=0;i<nint;i++)
    {
      x=iu[i]-x0;
      fiu+=wu[i]*exp(-0.5*(x*x)/(dx*dx));
    }
  as=1.0/fiu;
  // done.

  if (dp == 0) // only beam streaming with mu=1
    {
      for (int j=1;j<nu;j++)
	fu[j]=0.0;

      if (p0 == 2) // both directions
	{
	  fu[0]=0.5/(u[1]-u[0]);
	  fu[nu-1]=0.5/(u[1]-u[0]);
	}
	  else fu[0]=1.0/(u[1]-u[0]);
    }
  else  //beam with width dp centered at p0
    for (int j=0;j<nu;j++)
      {
	if (p0 == 2)
	  {
	    fu[j] = au/2.0*exp(-0.5*((u[j]-1.0)*(u[j]-1.0))/(dp*dp)) 
	      + au/2.0*exp(-0.5*((u[j]+1.0)*(u[j]+1.0))/(dp*dp));
	  }
	else 
	  {
	    p = u[j] - p0;
	    fu[j] = (dp >= 2) ? 0.5 : au*exp(-0.5*(p*p)/(dp*dp));
	  }
      }

  // isotropic 
  //  if (dp >= 2)
  //  for (int j=0;j<nu;j++)
  //    fu[j]=0.5;

  if (dx == 0) 
    {
      for (int k=0;k<ns;k++)
	fs[k]=0.0;
      fs[(ns-1)/2]=1.0/abs(s[1]-s[0]);
    }
  else 
    {
      for (int k=0;k<ns;k++)
	{
	  x = s[k] - x0;
	  fs[k] = as*exp(-0.5*x*x/(dx*dx)); 
	}
    }
  for (int i=0;i<ne;i++)
    for (int j=0;j<nu;j++)
      for (int k=0;k<ns;k++)
	{
	  index = i*nu*ns + j*ns + k;
	  pinj[index] += fe[i]*fu[j]*fs[k];
	}
}
///////////////////////////////////////////////////////////////////////////////

void tinj(double t, double dt, double t0, double tau, double *fptr, double *pinj, int n)
{

  double ft,tt; 

  tt=(t-t0);
  ft=(tau == 0) ? 1.0 : 1.0/(sqrt(2.0*PI)*tau)*exp(-0.5*(tt*tt)/(tau*tau))*dt;
    
  for (int i=0;i<n;i++)
    fptr[i] += pinj[i]*ft;

}

///////////////////////////////////////////////////////////////////////////////

void readBfield(double s[], double dbds[], int ns)
{
  string bfieldFileName="bfielddata.txt";
  ifstream bfieldFile;
  bfieldFile.open(bfieldFileName.c_str());
  if (!bfieldFile.is_open())
    {
      cerr << "ERROR: Could not open " << bfieldFile << " for reading." << endl;
      bfieldFile.clear();
      bfieldFile.close();
      exit(EXIT_FAILURE);
    }
  
  int npts;
  bfieldFile >> npts;
  double * posData = new double [npts];
  double * bfieldData = new double [npts];
  for (int i=0;i<npts;i++)
    bfieldFile >> posData[i];
  for (int i=0;i<npts;i++)
    bfieldFile >> bfieldData[i];
  
  for (int i=0;i<ns;i++)
    dbds[i]=linterpol(s[i],posData,bfieldData,npts);
    
  delete [] bfieldData;
  bfieldFile.clear();
  bfieldFile.close();
}

///////////////////////////////////////////////////////////////////////////////
double linterpol(double x, double xm[], double ym[], int n)
{
  /* linear interpolation */
  double r;

  for (int i=0; i < n; i++)
    {
      if ((x >= xm[i]) && (x <= xm[i+1]))
	{
	  r = ym[i]+(x-xm[i])*(ym[i+1]-ym[i])/(xm[i+1]-xm[i]); 
	  return r;
	}
    }
 
  r = ym[0];
  if (x >= xm[n-1])
    r = ym[n-1];

  return r;
}
////////////////////////////////////////////////////////////////////////////////
void progressbar( int percent , double time, double tnow, double estim)
{
/* ********************************************************************************
   progress bar code from http://nakkaya.com/2009/11/08/command-line-progress-bar/ 
   (many thanks)
   ******************************************************************************** */

  string bar;

  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }
  
  cout<< "\r" "[" << bar << "] ";
  cout.width( 3 );
  cout<< percent << "% "; 
  cout << "| secs to finish: ";
  cout.width(6);
  cout << time << "  | elapsed: "; 
  cout.width(6);
  cout << tnow << " | total estimated: ";
  cout.width(6);
  cout << estim;
  cout << std::flush;
  
}

/*
For an array of [ A ][ B ][ C ], we can index it thus:
(a * B * C) + (b * C) + (c * 1)
*/



