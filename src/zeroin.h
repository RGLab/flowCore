/*
Header file added to flowCore by John Ramey
The `R_zeroin` function used to be and is no longer apart of the R API.
See the following link for more information:
http://developer.r-project.org/blosxom.cgi/R-devel/2012/08/15#n2012-08-15
*/

#include <R_ext/Boolean.h>
#include <R_ext/RS.h>           /* F77_... */
#include <R_ext/BLAS.h>
#include <R.h>
#include <Rinternals.h>

double R_zeroin(double ax, double bx, double (*f)(double, void *), void *info, double *Tol, int *Maxit);
