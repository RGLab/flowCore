#include <iostream>
#include "logicle.h"
extern "C" {
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
/**
 * Logicle tranform/inverse transform wrapper function, makes use of the Logicle
 *  class provided Wayne Moore for the underlying calculation of the transformation.
 *
 * */

SEXP logicle_transform(SEXP input, SEXP T, SEXP W, SEXP M, SEXP A) {

    SEXP output;
	PROTECT(output = duplicate(input));
    try{
        Logicle *lg = new Logicle(asReal(T), asReal(W), asReal(M), asReal(A));
        for (int i = 0; i < length(output); i++) {
            REAL(output)[i] = lg->scale(REAL(output)[i]) * asReal(M) ;
        }
    }
    catch(const char * str){
      Rf_error("Logicle Exception: %s \n", str) ;
    }

    UNPROTECT(1);
    return(output);
}

SEXP invLogicle_transform(SEXP input, SEXP T, SEXP W, SEXP M, SEXP A){
	SEXP output;
    PROTECT(output = duplicate(input));
    try{
        Logicle *lg = new Logicle(asReal(T), asReal(W), asReal(M), asReal(A));
        for (int i = 0; i < length(output); i++) {
               REAL(output)[i] = lg->inverse(REAL(output)[i]/asReal(M))  ;
        }
    }
    catch(const char * str){
      Rf_error("Logicle Exception: %s \n", str) ;
    }
    UNPROTECT(1);
    return(output);
}

}// end of extern c 



