#include <Rinternals.h>
#include <R.h>
#include <R_ext/Applic.h>

struct biexponential_info {
	double a,b,c,d,f,w,y;
};

double biexponential_fn(double x,void*info) {
	struct biexponential_info *p = (struct biexponential_info *)info;
	double B = p->a*exp(p->b*(x-p->w))-p->c*exp(-p->d*(x-p->w))+p->f -p-> y;
	//return ((x>p->w) ? B : -B);
        return (B);
}

/*4/14/09 ngopalak: Updated biexponential transform so that it now calculates
the solution of a generic biexponential function 
S(x,a,b,c,d,f) = ae^(bx) -ce^(-dx) +f  instead of the logicle function
it was calculating earlier
*/
SEXP biexponential_transform(SEXP input,SEXP A,SEXP B,SEXP C,SEXP D,SEXP F,SEXP W,SEXP tol,SEXP maxit) {
	SEXP    output;
	struct biexponential_info params;
	
	int    i;
        int fail=0;
	double step,checkme=1;
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
		double x,y;
		int    MaxIt = asInteger(maxit);
		double res;
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

/*4/14/09 ngopalak: Updated Logicle transformation. the logicle transformation
 was implemented incorrectly.Additionally, the parameter argument f passed to 
this function was missing a multiplying
factor a
*/

SEXP logicle_transform(SEXP input,SEXP A,SEXP B,SEXP C,SEXP D,SEXP F,SEXP W,SEXP tol,SEXP maxit) {
	SEXP    output;
	struct biexponential_info params;
	int     i,fail=0;
	double step,checkme=1;
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
		double x,y;
		int    MaxIt = asInteger(maxit);
		double res;
		y = REAL(output)[i];
		params.y = (y >= params.w) ? y: 2*params.w-y;
		
                for(j=0,step=0.5; biexponential_fn(-step,(void*)&params)*biexponential_fn(step,(void*)&params) >0;step*=1.5,j+=1){  
                    if(j > MaxIt){
                        break;
                    }
                }
                res=R_zeroin(-step,step,biexponential_fn,(void*)&params,&Tol,&MaxIt); 
                REAL(output)[i] = (y>= params.w) ? res : -res;
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

SEXP invLogicle_transform(SEXP input,SEXP A,SEXP B,SEXP C,SEXP D,SEXP F,SEXP W){
	SEXP    output;
	struct biexponential_info params;
	int     i;
	params.a = asReal(A);
	params.b = asReal(B);
	params.c = asReal(C);
	params.d = asReal(D);
	params.f = asReal(F);
	params.w = asReal(W);
	params.y = 0 ;
	if(TYPEOF(input) != REALSXP) error("Input must be real values.");
	PROTECT(output = duplicate(input));
	for(i=0;i<length(output);i++) {
		double res;
		double x;
		x = REAL(output)[i];
		if( x >= 0){
			res = biexponential_fn(x,(void*)&params);
		
		}else{
			x = 2*params.w -x;
		    res = -biexponential_fn(x,(void*)&params);
		}
		REAL(output)[i] = res;
	
	}
	UNPROTECT(1);
	return output;
}
