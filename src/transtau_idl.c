#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "export.h"
#define csign(x) (x < 0.0 ? -1 : 1)
#define sqr(x) ((x)*(x))
#define BIGNEG -1000000
#define		Ls		1

/******************************************************************************/
IDL_LONG transtau_idl(int argc,void *argv[])
{

  IDL_LONG retval=1;
  float aa,t0,eqstate,boxside,omegam,omegav,weff,rmax;
  IDL_LONG Ng,Nmatch;
  IDL_LONG *iskew;
  float *ddat,*vdat,*tau,*tau_real,*tau_nov,*tau_cold;
  float *yskew,*zskew; 
  float *xqso,*yqso,*zqso; 
  float *theta,*phi,*psi,*tqso;
  float *wqso;

  // allocate internal variables 
  int ix,ik,ip,isim,im;
  float	u,u0,du,u1,du1,du_cold,sig_cold,bcold,tcold,t,tau_const;
  float rho,vel,b,sig,sig1,nh1,omegak;
  float r2,dx,dy,dz,adx,denom,rperp2,norm,ndotq,cos_psi;
  float dot_x,dot_y,dot_z,time_delay;
  float	Hz,mpart,t2v,vside;
  float	tmin=10.,bmin2=0.1,dmax=1e3,r2min=1e-6;
  double dNg;
  
  // Allocate pointers from IDL
  aa         =    *(float *)argv[0];  // scale factor
  t0         =    *(float *)argv[1];  // temperature in K
  eqstate    =    *(float *)argv[2];  // equantion of state
  boxside    =    *(float *)argv[3];  // size of box
  omegam     =    *(float *)argv[4];  // omega matter
  omegav     =    *(float *)argv[5];  // omega vacuum
  weff       =    *(float *)argv[6];  // w dark energy eq. state
  rmax       =    *(float *)argv[7]; // maximum distance to evalute signal
  Ng         = *(IDL_LONG *)argv[8];  // total number of spectral pixels
  ddat       =     (float *)argv[9];  // density array (Nspec X Ng)
  vdat       =     (float *)argv[10]; // velocity array
  tau        =     (float *)argv[11]; // tau array
  tau_nov    =     (float *)argv[12]; // tau array
  tau_real   =     (float *)argv[13]; // tau array
  tau_cold   =     (float *)argv[14]; // tau array
  Nmatch     = *(IDL_LONG *)argv[15]; // number of halo-skewer matches
  iskew      =  (IDL_LONG *)argv[16]; // Skewer index
  yskew      =     (float *)argv[17]; // Skewer y position
  zskew      =     (float *)argv[18]; // Skewer z position
  xqso       =     (float *)argv[19]; // halo x position (line-of-sight)
  yqso       =     (float *)argv[20]; // halo y position
  zqso       =     (float *)argv[21]; // halo z position
  theta      =     (float *)argv[22]; // quasar emission polar angle
  phi        =     (float *)argv[23]; // quasar emission azimuthal angle
  psi        =     (float *)argv[24]; // quasar emission cone half-angle
  tqso       =     (float *)argv[25]; // quasar time in units mpc/h proper
  wqso       =     (float *)argv[26]; // ratio gam_qso/gam_uvb @ 1 proper mpc/h

  omegak = 1.0 - omegam - omegav;
  Hz    = 100*sqrt(omegam/aa/aa/aa+omegav/pow(aa,3*(1+weff))+omegak/aa/aa);
  vside = aa*Hz*boxside;	/* Side length in km/s */
  mpart = 1.672634e-4;		/* Hydrogen */
  t2v   = 1.380658/mpart/1.e6;	/* Kelvin to km/s */
  dNg   = (double)Ng;
  // cold lines are for case of no thermal broadening
  tcold = 1.5e4; // temperature width of 100K or 1.28 km/s
  bcold = sqrt(2*t2v*tcold);
  tau_const = vside/dNg/sqrt(M_PI);
  
  //fprintf(stderr,"aa=%f t0=%f eqstate=%f boxside=%f\n",aa,t0,eqstate,boxside);
  //fprintf(stderr,"omegam=%f omegav=%f weff=%f Nspec=%d Ng=%d\n"
  //,omegam,omegav,weff,Nspec,Ng);
  //return 1;

  //fprintf(stderr,"iskew=%d im=%d x=%f xqso=%f dx=%f cos_psi=%f ndotq=%f r2=%f time_delay=%f tqso=%f\n",iskew[im]
  //,im,(double)ix*boxside/dNg,xqso[im],dx,cos_psi,ndotq,r2,time_delay,tqso[im]);

  for (ik=0; ik< Ng*Nmatch; ik++) {
    tau[ik]=0.0;
    tau_nov[ik]=0.0;
    tau_real[ik]=0.0;
    tau_cold[ik]=0.0;
  }
  
  for(im=0; im<Nmatch; im++) {
    // This counter does not reset the line. 
    fprintf(stderr,"Working on pair %6d of %6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",im,Nmatch);
    dy = yskew[im] - yqso[im];
    dz = zskew[im] - zqso[im]; 
    rperp2 = sqr(dy) + sqr(dz);
    dot_x = cos(theta[im]);
    dot_y = sin(theta[im])*cos(phi[im]);
    dot_z = sin(theta[im])*sin(phi[im]);
    // if psi is negative this means isotropic
    if(psi[im] < 0) cos_psi=BIGNEG;
    else cos_psi=cos(psi[im]);
    for (ix=0; ix<Ls*Ng; ix++) { /* Real space position */
      dx  = (double)ix*boxside/dNg - xqso[im];
      adx = fabs(dx); 
      if((boxside-adx)<adx) dx=-csign(dx)*(boxside-adx);
      if(fabs(dx) > rmax) continue; 
      // if out of range of interest don't do convo
      isim  = ix+iskew[im]*Ls*Ng;
      rho = ddat[isim]/(1.0+ddat[isim]/dmax);	/* Maxm density */
      vel = vdat[isim];
      if (rho>0) {
	t   = t0*pow(rho,eqstate-1.0)+tmin;	/* Minm T to avoid overflow */
	b   = sqrt(2*t2v*t + bmin2);		/* Minm b to avoid overflow */
	nh1 = rho*rho*pow(t/t0,-0.7);
	// transverse proximity effect stuff
	r2 = rperp2 + sqr(dx);
	norm = sqrt(r2);
	if(norm > 1.0e9) norm=1.0e9; // prevent overflow
	ndotq = dot_x*dx/norm + dot_y*dy/norm + dot_z*dz/norm;
	time_delay = -dx + sqrt(r2); 
	// positive dx means point in foreground of quasar toward observer
	if((ndotq >= cos_psi) && (time_delay <= tqso[im])) {
	  denom = (1.0 + wqso[im]/(r2 + r2min)); // rmin avoid overflow
	  if(denom > 1.0e9) denom=1.0e9; // prevent overflow
	  nh1 = nh1/denom;
	  //tau[Ls*Ng*im+ix] += 1.0; // TESTING!!!
	  //return 1;
	}
	tau_real[Ls*Ng*im+ix] += nh1;
	for (ip=-Ls*Ng; ip<=Ls*Ng; ip+=Ls*Ng) {
	  u0  = (double)(ix+ip)/dNg*vside;
	  u1  = u0; // u1 has no peculiar velocities
	  u0 += vel;
	  for (ik=0; ik<Ls*Ng; ik++) {  /* Spectral position. */
	    u   = (double)ik/dNg*vside;
	    du  = fabs(u-u0)/b;
	    if (du<6.0) {
	      sig = exp(-du*du)/b;
	      tau[Ls*Ng*im+ik] += sig*nh1;
	    }
	    // compute tau with peculiar vel set to zero
	    du1 = fabs(u-u1)/b;  
	    if (du1<6.0) {
	      sig1 = exp(-du1*du1)/b;
	      tau_nov[Ls*Ng*im+ik] += sig1*nh1;
	    }
	    // compute tau with thermal motions set to small number
	    du_cold  = fabs(u-u0)/bcold;
	    if (du_cold<6.0) {
	      sig_cold = exp(-du_cold*du_cold)/bcold;
	      tau_cold[Ls*Ng*im+ik] += sig_cold*nh1;
	    }
	  }
	}
      }
    }
  }
  
  // multiply all by tau_const
  for (ik=0; ik<Ng*Nmatch; ik++) {
    tau[ik]=tau[ik]*tau_const;
    tau_nov[ik]=tau_nov[ik]*tau_const;
    tau_real[ik]=tau_real[ik]*tau_const;
    tau_cold[ik]=tau_cold[ik]*tau_const;
  }

  return retval;
}
