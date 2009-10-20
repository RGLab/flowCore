#include <Rinternals.h>
#include <R.h>
#include <R_ext/Applic.h>

struct biexponential_info {
	double a,b,c,d,f,w,y;
};

double biexponential_fn(double x,void*info) {
	struct biexponential_info *p = (struct biexponential_info *)info;
	double B = p->a*exp(p->b*(x-p->w))-p->c*exp(-p->d*(x-p->w))+p->f -p-> y;
    return (B);
}

struct sfun_info{
	double m,w,p,t,a,r;
};

/*
Logicle function calculation 
*/
double strans_fn(double y,void*info){
	struct sfun_info  *k = (struct sfun_info *)info;
	double tmp ;
	tmp = (y < (k->w +k->a)) ? (k->w + k->a - y) : (y - k->w - k->a) ; 
	tmp =  k->t*pow(10,-1*(k->m - k->w - k->a ))*(pow(10,tmp)-k->p*k->p*pow(10,-tmp/k->p)+ k->p*k->p-1);	
	return((y < (k->w +k->a)) ? -1*tmp -k->r  :tmp-k->r);
}

/*10/15/09 ngopalak: Updated biexponential transform so that it now calculates
the solution of a generic biexponential function 
S(x,a,b,c,d,f) = ae^(bx) -ce^(-dx) +f  instead of the logicle function
it was calculating earlier
*/
SEXP biexponential_transform(SEXP input,SEXP A,SEXP B,SEXP C,SEXP D,SEXP F,SEXP W,SEXP tol,SEXP maxit) {
	SEXP    output;
	struct biexponential_info params;
	
	int    i;
    int fail=0;
	double step;
	params.a = asReal(A);
	params.b = asReal(B);
	params.c = asReal(C);
	params.d = asReal(D);
	params.f = asReal(F);
    params.w = asReal(W);

	if(TYPEOF(input) != REALSXP) error("Input must be real values.");
	PROTECT(output = duplicate(input));
	for(i=0;i<length(output);i++) {
		int    j;
		double Tol   = asReal(tol);
		int    MaxIt = asInteger(maxit);
		params.y = REAL(output)[i];
                    for(j=0,step=0.5; biexponential_fn(-step,(void*)&params)*biexponential_fn(step,(void*)&params) >0;step*=1.5,j+=1){  
                        if(j > MaxIt){
                              break;
                          }
                    }
                    REAL(output)[i]=R_zeroin(-step,step,biexponential_fn,(void*)&params,&Tol,&MaxIt); 
                    if(MaxIt==-1){ 
                        fail=fail+1;
                    }
        }
	if(fail>0)
		warning("%d values of %d have not converged.",fail,length(output));
	UNPROTECT(1);
	return output;
}

/*
10/20/09 ngopalak: Updated the logicle transform so that it now uses a separate logicle function strans_info
instead of the generic biexponential function being used earlier. This change was necessary to 
simplify the implementation and improve code readability.

*/
SEXP logicle_transform(SEXP input,SEXP M,SEXP W,SEXP P,SEXP T,SEXP A,SEXP tol,SEXP maxit) {
	SEXP    output;
	struct sfun_info params;
	int     i,fail=0;
	double step;
	params.m = asReal(M);
	params.w = asReal(W);
	params.p = asReal(P);
	params.t = asReal(T);
	params.a = asReal(A);
	double Tol   = asReal(tol);
	int    MaxIt ;	
	if(TYPEOF(input) != REALSXP) error("Input must be real values.");
	PROTECT(output = duplicate(input));
	for(i=0;i<length(output);i++) {
		int    j;
		double res;
		double y = REAL(output)[i];
		MaxIt = asInteger(maxit);
      	params.r =y;
 		for(j=0,step=0.5; strans_fn(-step,(void*)&params)*strans_fn(step,(void*)&params) >0;step*=1.5,j+=1){  
            if(j > MaxIt){
                break;
            }
        }
        res=R_zeroin(-step,step,strans_fn,(void*)&params,&Tol,&MaxIt); 
	    REAL(output)[i] = res;
        if(MaxIt==-1){ 
            fail=fail+1;
        }
	
    }
	if(fail>0)
		warning("%d values of %d have not converged.",fail,length(output));
	UNPROTECT(1);
	return output;
}

/*10/06/09 ngopalak: Added a function to calculate the inverse of the logicle 
transformation, making use of the previously defined biexponential_fn to calculate 
the inverse
*/

SEXP invLogicle_transform(SEXP input,SEXP M,SEXP W,SEXP P,SEXP T,SEXP A){
	SEXP    output;
	struct sfun_info params;
	int     i;
	params.m = asReal(M);
	params.w = asReal(W);
	params.p = asReal(P);
	params.t = asReal(T);
	params.a = asReal(A);
	params.r = 0 ;
	if(TYPEOF(input) != REALSXP) error("Input must be real values.");
	PROTECT(output = duplicate(input));
	for(i=0;i<length(output);i++) {
		double x;
		x = REAL(output)[i];
		REAL(output)[i] =  strans_fn(x,(void*)&params);
	}
	UNPROTECT(1);
	return output;
}
