#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "export.h"

#define		Ls		1

/******************************************************************************/
IDL_LONG maketau_idl(int argc,void *argv[])
{

  IDL_LONG retval=1;
  float aa,t0,eqstate,boxside,omegam,omegav,weff;
  IDL_LONG Nspec,Ng;
  float *ddat,*vdat,*tau,*tau_nov,*tau_real,*tau_cold;
  // allocate internal variables 
  int ix,im,is,ik,ip;
  float	u,u0,du,u1,du1,du_cold,bcold,tcold,sig_cold,tau_const;
  float t,rho,vel,b,sig,sig1,nh1,omegak;
  float	Hz,mpart,t2v,vside;
  float	tmin=10.,bmin2=0.1,dmax=1e3;
  double dNg;

  // Allocate pointers from IDL
  aa         =    *(float *)argv[0];  // scale factor
  t0         =    *(float *)argv[1];  // temperature in K
  eqstate    =    *(float *)argv[2];  // equantion of state
  boxside    =    *(float *)argv[3];  // size of box
  omegam     =    *(float *)argv[4]; // omega matter
  omegav     =    *(float *)argv[5]; // omega vacuum
  weff       =    *(float *)argv[6]; // w
  Nspec      = *(IDL_LONG *)argv[7]; // total num spectra is Nspec
  Ng         = *(IDL_LONG *)argv[8]; // total number of spectral pixels
  ddat       =     (float *)argv[9]; // density array (Nspec X Ng)
  vdat       =     (float *)argv[10]; // velocity array
  tau        =     (float *)argv[11]; // tau array
  tau_nov    =     (float *)argv[12]; // tau array
  tau_real   =     (float *)argv[13]; // tau array
  tau_cold   =     (float *)argv[14]; // tau array
  omegak = 1.0 - omegam - omegav;
  Hz    = 100*sqrt(omegam/aa/aa/aa+omegav/pow(aa,3*(1+weff))+omegak/aa/aa);
  vside = aa*Hz*boxside;	/* Side length in km/s */
  mpart = 1.672634e-4;		/* Hydrogen */
  t2v   = 1.380658/mpart/1.e6;	/* Kelvin to km/s */
  dNg   = (double)Ng;
  // cold lines are for case of no thermal broadening
  tcold = 1.5e4; // temperature width of 15,000 or 15.7 km/s, about pix scale
  bcold = sqrt(2*t2v*tcold);
  tau_const = vside/dNg/sqrt(M_PI);
  //fprintf(stderr,"aa=%f t0=%f eqstate=%f boxside=%f\n",aa,t0,eqstate,boxside);
  //fprintf(stderr,"omegam=%f omegav=%f weff=%f Nspec=%d Ng=%d\n",omegam,omegav,weff,Nspec,Ng);
  //return 1;

  for (ik=0; ik<Ls*Ng*Nspec; ik++) {
    tau[ik]=0.0;
    tau_nov[ik]=0.0;
    tau_real[ik]=0.0;
    tau_cold[ik]=0.0;
  }
  
  for (im=0; im<Nspec; im++) {
    // This counter does not reset the line. 
    fprintf(stderr,"Working on pair %6d of %6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",im,Nspec);
    for (ix=0; ix<Ls*Ng; ix++) { /* Real space position */
      is  = ix+im*Ls*Ng;
      rho = ddat[is]/(1.0+ddat[is]/dmax);	/* Maxm density */
      vel = vdat[is];
      if (rho>0) {
	t   = t0*pow(rho,eqstate-1.0)+tmin;	/* Minm T to avoid overflow */
	b   = sqrt(2*t2v*t + bmin2);		/* Minm b to avoid overflow */
	nh1 = rho*rho*pow(t/t0,-0.7);
	tau_real[is] += nh1;
	for (ip=-Ls*Ng; ip<=Ls*Ng; ip+=Ls*Ng) {
	  u0  = (double)(ix+ip)/dNg*vside;
	  u1  = u0; // u1 has no peculiar velocity
	  u0 += vel;
	  for (ik=0; ik<Ls*Ng; ik++) {  /* Spectral position. */
	    u   = (double)ik/dNg*vside;
	    du  = fabs(u-u0)/b;
	    if (du<6.0) {
	      sig = exp(-du*du)/b;
	      tau[Ls*Ng*im+ik] += sig*nh1;
	    }
	    // compute the tau with peculiar vel set to zero
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
  for (ik=0; ik<Ls*Ng*Nspec; ik++) {
    tau[ik]=tau[ik]*tau_const;
    tau_nov[ik]=tau_nov[ik]*tau_const;
    tau_real[ik]=tau_real[ik]*tau_const;
    tau_cold[ik]=tau_cold[ik]*tau_const;
  }
  
  
  return retval;
}
