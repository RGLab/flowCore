/*
 * Gopalakrishnan N 10/14/2008
 */

#include <Rcpp.h>

#include <stdlib.h>
using namespace Rcpp;


void inPolytope_c(double *data, double *A, double *b, int nRowData, int nRowA, int nColA, 
                  std::vector<bool> & result) {

  int i, j, k;
  float sigma[nRowA*nRowData];
  float temp=0;

  for(i=0;i<nRowData;i++)
  {
    result[i]=1;
  }
  for ( k=0;k<nRowData;k++)
  {
    for( i=0;i<nRowA;i++)
    {      temp=0;
    for ( j=0;j<nColA;j++)
    {
      temp+=data[k+j*(nRowData)]*A[i+j*(nRowA)];
    }
    sigma[i+nRowA*k]=temp+b[i];
    }
  }

  for(k=0;k<nRowData;k++)
  {
    for(i=0;i<nRowA;i++)
    {  
       if(sigma[i+nRowA*k]>0)        //if(sigma[i][k]>0)
       {
        result[k]=0;
        break;
      }
 
    }
  }

 
}


/*-----------------------------------------------------------------
   interface to R with arguments:
     data :    matrix of numerics
     A:   matrix , number of columns = number of columns of data
    b: matrix with 1 column,number of rows = number of rows of A 
------------------------------------------------------------------*/
//[[Rcpp::export]]
std::vector<bool> inPolytope(NumericMatrix data, NumericMatrix A, NumericVector b)
{ 
  int nRowData = data.nrow();
  std::vector<bool> result(nRowData);   
  int nRowA = A.nrow();
  int nColA = A.ncol(); 

  
  if(b.length()!=nRowA)
      stop("Invalid argument 'b': must be a real vector of length 'nrow(A)'."); 
  
  inPolytope_c(REAL(data), REAL(A), REAL(b), nRowData, nRowA, nColA, result);


  return result;
}
