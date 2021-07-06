#include <cpp11.hpp>
#include <stdlib.h>

[[cpp11::register]] cpp11::doubles_matrix poly_centroid(
    cpp11::doubles_matrix verts) {
  int nrv = verts.nrow();
  if(nrv < 2 || verts.ncol() != 2)
    cpp11::stop("Argument 'vertices' must be numeric matrix of two columns and at least\ntwo rows specifying vertices of a polygon on a two-dimensional plane");
  double area = 0.0, cx = 0.0, cy = 0.0, cross = 0.0;
  for(int i = 0; i < nrv; i++){
    cross = (verts(i, 0)*verts((i+1) % nrv, 1) - verts((i+1) % nrv, 0)*verts(i, 1));
    area += cross;
    cx += (verts(i, 0) + verts((i+1) % nrv, 0))*cross;
    cy += (verts(i, 1) + verts((i+1) % nrv, 1))*cross;
  }
  cx /= (3*area);
  cy /= (3*area);
  cpp11::writable::doubles_matrix centroid(1, 2);  // a 1 x 2 matrix
  centroid(0, 0) = cx;
  centroid(0, 1) = cy;
  return centroid;
}

