/*
 * F. Hahne  12/10/2005
 */
#include "cpp11.hpp"

#include <cytolib/in_polygon.hpp>



/*-----------------------------------------------------------------
   interface to R with arguments:
     data :    matrix of numerics
     vertices:   matrix with vertices
------------------------------------------------------------------*/
[[cpp11::register]] std::vector<bool> inPolygon(
    cpp11::doubles_matrix data, cpp11::doubles_matrix vertices) {
 

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
      cpp11::stop("Argument 'points' must be numeric matrix of two columns and at least\none row specifiying points on a two-dimensional plane");

    if(nrv < 2 || vertices.ncol() != 2)
      cpp11::stop("Argument 'vertices' must be numeric matrix of two columns and at least\ntwo rows specifying vertices of a polygon on a two-dimensional plane");

     /* Do it! */
    vector<cytolib::CYTO_POINT> points(nrv);
    for(int i = 0; i < nrv; i++)
    {
      points[i].x = vertices(i, 0);  // vertices must have two columns
      points[i].y = vertices(i, 1);
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
