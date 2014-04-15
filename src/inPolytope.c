/*
 * Gopalakrishnan N 10/14/2008
 */

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h> 

#include <stdlib.h>


void inPolytope_c(double *data, double *A, double *b, int nRowData, int nRowA, int nColA, 
		 int *result) {

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

SEXP inPolytope(SEXP _data, SEXP _A, SEXP _b)
{ 
  SEXP result;      /* return value: a integer vector*/
  SEXP dimData;
  SEXP dimA;

  int dimB;
  double *data;
  double *A;
  double *b;

  int nRowData, nRowA, nColA;  /* dimensions of data    */

  /* check input argument _data */
  PROTECT(dimData = getAttrib(_data, R_DimSymbol));
  if(((!isReal(_data)) & !isInteger(_data)) | isNull(dimData) | (LENGTH(dimData)!=2))
     error("Invalid argument 'data': must be a real matrix.");
  data = REAL(AS_NUMERIC(_data));
  nRowData = INTEGER(dimData)[0];
  UNPROTECT(1);          
  /* done with dimData */

  /* check input argument _A */
  PROTECT(dimA = getAttrib(_A, R_DimSymbol));
   if((!isReal(_A)) | isNull(dimA) | (LENGTH(dimA)!=2))
      error("Invalid argument 'A': must be a real matrix."); 
  A = REAL(AS_NUMERIC(_A));
  nRowA  = INTEGER(dimA)[0];
  nColA  = INTEGER(dimA)[1];
  UNPROTECT(1);          
  /* done with A */

 /* check input argument _b */
  dimB = (LENGTH(_b));
  if((!isReal(_b)) | (dimB!=nRowA))
      error("Invalid argument 'b': must be a real vector of length 'nrow(A)'."); 
  b = REAL(AS_NUMERIC(_b));         
  /* done with b */

 
  /* allocate memory for return values */
  PROTECT(result = allocVector(INTSXP, nRowData));

  inPolytope_c(data, A, b, nRowData, nRowA, nColA, INTEGER(result));

  UNPROTECT(1); 

  return result;
}
