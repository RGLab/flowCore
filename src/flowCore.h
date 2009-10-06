#include <R.h>
#include <Rinternals.h>

SEXP biexponential_transform(SEXP input,SEXP A,SEXP B,SEXP C,SEXP D,SEXP F,SEXP G,SEXP tol,SEXP maxit);
SEXP logicle_transform(SEXP input,SEXP A,SEXP B,SEXP C,SEXP D,SEXP F,SEXP G,SEXP tol,SEXP maxit);
SEXP invLogicle_transform(SEXP input,SEXP A,SEXP B,SEXP C,SEXP D,SEXP F,SEXP G);
SEXP inPolygon(SEXP _data, SEXP _vertices);
SEXP inPolytope(SEXP _data, SEXP _A, SEXP _b);

//  SEXP biexponential_transform(SEXP *,SEXP *,SEXP *,SEXP *,SEXP *,SEXP *,SEXP *,SEXP *,SEXP *);
//  SEXP inPolygon(SEXP *, SEXP *);
//  SEXP inPolytope(SEXP *, SEXP *, SEXP *);
