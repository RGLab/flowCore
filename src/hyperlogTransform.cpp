/**
  Hyperlog transformation added by Josef Spidlen.
  This hyperlog implementation is based on Java reference 
  implementation that is part of the full Gating-ML 2.0
  specification. The Java reference implementation has
  been provided by Wayne Moore, see hyperlog.notice.html
  for details. Josef Spidlen ported it to C/CPP and 
  integrated it with R/flowCore.
*/

#include <Rcpp.h>
#include "hyperlog.h"


/**
 * Hyperlog tranform/inverse transform wrapper function, makes use of the Hyperlog
 * class adapted from Wayne Moore's Java Hyperlog implementation for the underlying
 * calculation of the transformation.
 **/
//[[Rcpp::export]]
std::vector<double> hyperlog_transform(std::vector<double> input, double T, double W, double M, double A, bool isInverse) {
	unsigned nLen = input.size();

	    try{
	    	Hyperlog lg = Hyperlog(T, W, M, A);
				for (unsigned i = 0; i < nLen; i++) {
					if(isInverse)
						input.at(i) = lg.inverse(input.at(i));
					else
						input.at(i) = lg.scale(input.at(i));
				}
	    	}
	    catch(const char * str){
	      std::string tmp= "Hyperlog Exception: ";

	    	Rcpp::stop(tmp.append(str));
	    }

        return(input);
}




