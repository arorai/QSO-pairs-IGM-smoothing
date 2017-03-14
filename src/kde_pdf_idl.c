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
IDL_LONG kde_pdf_idl(int argc,void *argv[])
{

  int ii,ipdf;
  IDL_LONG retval=1;
  float hopt,xarg;
  IDL_LONG ni,npdf;
  double temp,norm;
  float *xpdf,*xi,*pdf;
  
  // Allocate pointers from IDL
  npdf       =    *(IDL_LONG *)argv[0];  // number of pdf points
  ni         =    *(IDL_LONG *)argv[1];  // number of data points
  hopt       =    *(float *)argv[2];     // smoothing length
  xpdf       =     (float *)argv[3]; // values to evaluate pdf at
  xi         =     (float *)argv[4]; // data to evaluate pdf for
  pdf        =     (float *)argv[5]; // pdf evaluated at xpdf

  norm=1.0/sqrt(2.0*M_PI)/(double)ni/(double)hopt;
  for(ipdf=0; ipdf<npdf; ipdf++) {
    temp=0.0;
    for (ii=0; ii<ni; ii++) { 
      xarg=(xpdf[ipdf]-xi[ii])/hopt;
      if(abs(xarg) < 10.0) temp += exp(-xarg*xarg/2.0);
    }
    pdf[ipdf]=temp*norm;
  }
  
  return retval;
}
