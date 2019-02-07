#include <Rcpp.h>
#include <stdlib.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector poly_centroid(NumericMatrix verts) {
//Improvement -- Use RcppArmadillo::det (for cross) or even better just use boost::geometry::centroid
  int nrv = verts.nrow();
  if(nrv < 2 || verts.ncol() != 2)
    stop("Argument 'vertices' must be numeric matrix of two columns and at least\ntwo rows specifying vertices of a polygon on a two-dimensional plane");
  double area = 0.0, cx = 0.0, cy = 0.0, cross = 0.0;
  for(int i = 0; i < nrv; i++){
    cross = (verts(i, 0)*verts((i+1) % nrv, 1) - verts((i+1) % nrv, 0)*verts(i, 1));
    area += cross;
    cx += (verts(i, 0) + verts((i+1) % nrv, 0))*cross;
    cy += (verts(i, 1) + verts((i+1) % nrv, 1))*cross;
  }
  cx /= (3*area);
  cy /= (3*area);
  NumericVector centroid = NumericVector::create(cx, cy);
  return centroid;
}

