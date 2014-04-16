/**
  Hyperlog transformation added by Josef Spidlen.
  This hyperlog implementation is based on Java reference 
  implementation that is part of the full Gating-ML 2.0
  specification. The Java reference implementation has
  been provided by Wayne Moore, see hyperlog.notice.html
  for details. Josef Spidlen ported it to C/CPP and 
  integrated it with R/flowCore.
*/

#include <iostream>
#include "hyperlog.h"
extern "C" {
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

/**
 * Hyperlog tranform/inverse transform wrapper function, makes use of the Hyperlog
 * class adapted from Wayne Moore's Java Hyperlog implementation for the underlying
 * calculation of the transformation.
 **/

SEXP hyperlog_transform(SEXP input, SEXP T, SEXP W, SEXP M, SEXP A) {
    SEXP output;
    PROTECT(output = duplicate(input));
    try{
        Hyperlog *hplg = new Hyperlog(asReal(T), asReal(W), asReal(M), asReal(A));
        for (int i = 0; i < length(output); i++) {
            REAL(output)[i] = hplg->scale(REAL(output)[i]);
        }
        if (hplg != NULL) delete hplg;
    }
    catch(const char *str){
        Rf_error("Hyperlog Exception: %s \n", str);
    }

    UNPROTECT(1);
    return(output);
}

SEXP invHyperlog_transform(SEXP input, SEXP T, SEXP W, SEXP M, SEXP A) {
    SEXP output;
    PROTECT(output = duplicate(input));
    try{
        Hyperlog *hplg = new Hyperlog(asReal(T), asReal(W), asReal(M), asReal(A));
        for (int i = 0; i < length(output); i++) {
               REAL(output)[i] = hplg->inverse(REAL(output)[i])  ;
        }
        if (hplg != NULL) delete hplg;
    }
    catch(const char *str){
      Rf_error("Hyperlog Exception: %s \n", str) ;
    }
    UNPROTECT(1);
    return(output);
}

} // end of extern c 



