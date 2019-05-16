/*
 * F. Hahne  12/10/2005
 */
#include <Rcpp.h>

#include <cytolib/in_polygon.hpp>

using namespace Rcpp;



/*-----------------------------------------------------------------
   interface to R with arguments:
     data :    matrix of numerics
     vertices:   matrix with vertices
------------------------------------------------------------------*/
//[[Rcpp::export]]
std::vector<bool> inPolygon(NumericMatrix data, NumericMatrix vertices)
{ 
  
  int nrd = data.nrow();
  int nrv = vertices.nrow();
  
  /* check input argument _data */
  if(nrd == 0)
  {
    std::vector<bool> res(nrd, false);
    return(res);
  }else
  {
    
    if(data.ncol() != 2)
      stop("Argument 'points' must be numeric matrix of two columns and at least\none row specifiying points on a two-dimensional plane");

    if(nrv < 2 || vertices.ncol() != 2)
      stop("Argument 'vertices' must be numeric matrix of two columns and at least\ntwo rows specifying vertices of a polygon on a two-dimensional plane");

     /* Do it! */
    vector<cytolib::POINT> points(nrv);
    for(int i = 0; i < nrv; i++)
    {
    	points[i].x = vertices[i];
    	points[i].y = vertices[i + nrv];
    }
    double * xdata = REAL(data);
    double * ydata = xdata + nrd;
    INDICE_TYPE parentInd(nrd);
    for(int i = 0; i < nrd; i++)
        parentInd[i] = i;
    INDICE_TYPE resInd;
    resInd.reserve(nrd);
    cytolib::in_polygon(xdata, ydata, points, parentInd, false, resInd);
    std::vector<bool> res(nrd);
    for(unsigned i = 0; i < resInd.size(); i++)
	   res[resInd[i]] = true;
    return(res);           
}
    
  
}
