#include <R.h>
#include <Rinternals.h>

SEXP biexponential_transform(SEXP input,SEXP A,SEXP B,SEXP C,SEXP D,SEXP F,SEXP G,SEXP tol,SEXP maxit);
SEXP logicle_transform(SEXP input,SEXP M,SEXP W,SEXP P,SEXP T,SEXP A,SEXP tol,SEXP maxit);
SEXP invLogicle_transform(SEXP input,SEXP M,SEXP W,SEXP P,SEXP T,SEXP A);
SEXP inPolygon(SEXP _data, SEXP _vertices);
SEXP inPolytope(SEXP _data, SEXP _A, SEXP _b);

