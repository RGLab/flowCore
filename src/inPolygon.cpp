/*
 * F. Hahne  12/10/2005
 */
#include <Rcpp.h>

#include <stdlib.h>
using namespace Rcpp;
/*-----------------------------------------------------------------
internal c function for calculation of pAUCs
-----------------------------------------------------------------*/

double min(double val1, double val2){
    if(val1 <= val2){
	return(val1);
    }
    else{
	return(val2);
    }
}

double max(double val1, double val2){
    if(val1 >= val2){
	return(val1);
    }
    else{
	return(val2);
    }
}



void inPolygon_c(double *data, int nrd, 
            double *vertices, int nrv, std::vector<bool> & res) {

  int i, j, counter;
  double xinters;
  double p1x, p2x, p1y, p2y;
  //find max py
  double p_y_max = vertices[nrv];
  for(i = 1; i < nrv; i++)
	  p_y_max = max(vertices[i + nrv], p_y_max);

  for(i=0; i<nrd; i++)
  {//iterate over points
    p1x=vertices[0];
    p1y=vertices[nrv];
    counter=0;
    for(j=1; j < nrv+1; j++)
    {// iterate over vertices 
      /*p1x,p1y and p2x,p2y are the endpoints of the current vertex*/
      if (j == nrv){//the last vertice must "loop around"
      	p2x = vertices[0];
      	p2y = vertices[0+nrv];
      }//if
      else
      {
      	p2x = vertices[j]; 
      	p2y = vertices[j+nrv];
      }//else
      /*if horizontal ray is in range of vertex find the x coordinate where
	ray and vertex intersect*/
      double p_y_bottom = p1y;
      double p_y_top = p2y;
      bool is_p2_top = true;
      if(p_y_bottom > p_y_top)
      {
    	  std::swap(p_y_bottom, p_y_top);
    	  is_p2_top = false;
      }
      double p_x_left = p1x;
      double p_x_right = p2x;
      if(p_x_left > p_x_right)
    	  std::swap(p_x_left, p_x_right);
      double data_x = data[i];
      double data_y = data[i + nrd];
      if(data_y >= p_y_bottom && data_y < p_y_top &&
         data_x <= p_x_right && p2y != p1y)
      {
  	      xinters = (data_y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x;
  	      /*if intersection x coordinate == point x coordinate it lies on the
  	      boundary of the polygon, which means "in"*/
        	if(xinters==data_x)
      	  {
        	  counter=1;
        	  break;
        	}
        	/*count how many vertices are passed by the ray*/
        	if (xinters > data_x)
      	  {
        	  counter++;
        	}
      }
      else if(data_y == p_y_max)//handle cell that is at the same y-level as the top vertex/edge
      {
    	  if(p_y_top == data_y)//one end of edge reach the same y as top
    	  {
    		  if(p_y_bottom == data_y)//horizontal top edge
    		  {
    			  counter = data_x >= p_x_left && data_x <= p_x_right;//whether on the edge

    		  }
    		  else//check if on the top vertex
    			  counter = data_x == (is_p2_top?p2x:p1x);
    	  }
    	  break;
      }
      p1x=p2x;
      p1y=p2y;
    }
    /*uneven number of vertices passed means "in"*/
    res[i] = counter % 2 > 0;
    
  }
}



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
    std::vector<bool> res(nrd);
    
    if(data.ncol() != 2)
      stop("Argument 'points' must be numeric matrix of two columns and at least\none row specifiying points on a two-dimensional plane");

    if(nrv < 2 || vertices.ncol() != 2)
      stop("Argument 'vertices' must be numeric matrix of two columns and at least\ntwo rows specifying vertices of a polygon on a two-dimensional plane");

     /* Do it! */
   inPolygon_c(REAL(data), nrd, REAL(vertices), nrv, res);
    return(res);           
}
    
  
}
