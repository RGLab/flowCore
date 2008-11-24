#include <Rinternals.h>
#include <R.h>
#include <R_ext/Applic.h>


struct biexponential_info {
	double a,b,c,d,f,g,y;
};

double biexponential_fn(double x,void*info) {
	struct biexponential_info *p = (struct biexponential_info *)info;
	int    sig;
	if(x>p->g) { sig=1;x=x-p->g;} else {sig=0;x=p->g-x;}
	double B = p->a*exp(p->b*x)-p->c*exp(-p->d*x)+p->f;
	
	return ((sig==1) ? B : -B)-p->y;
}


SEXP biexponential_transform(SEXP input,SEXP A,SEXP B,SEXP C,SEXP D,SEXP F,SEXP G,SEXP tol,SEXP maxit) {
	SEXP    output;
	struct biexponential_info params;
	
	int     i,fail=0;
	
	params.a = asReal(A);
	params.b = asReal(B);
	params.c = asReal(C);
	params.d = asReal(D);
	params.f = asReal(F);
	params.g = asReal(G);
	
	if(TYPEOF(input) != REALSXP) error("Input must be real values.");
	PROTECT(output = duplicate(input));
	for(i=0;i<length(output);i++) {
		int    j;
		double Ax,Bx;
		double Tol   = asReal(tol);
		double x;
		int    MaxIt = asInteger(maxit);
		params.y = REAL(output)[i];
		if(params.y!= 0) {
		  x = params.y < 0 ? -log(-params.y)+params.g:log(params.y)-params.g;
		  //Find some starting parameters
		  for(j=0,Ax=.5;biexponential_fn(x-Ax,(void*)&params)>0;Ax*=10) if(++j > MaxIt) break;
		  for(j=0,Bx=.5;biexponential_fn(x+Bx,(void*)&params)<0;Bx*=10) if(++j > MaxIt) break;
		  REAL(output)[i] = R_zeroin(x-Ax,x+Bx,biexponential_fn,(void*)&params,&Tol,&MaxIt); 
		  if(MaxIt==0) 
		    fail++;
		}
		else{
		  REAL(output)[i] = 0;
		}
	}
	if(fail>0)
		warning("%d values of %d have not converged.",fail,length(output));
	UNPROTECT(1);
	return output;
}
