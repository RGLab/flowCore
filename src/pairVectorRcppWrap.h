/*
 * pairVectorRcppWrap.h
 *
 *  Created on: Feb 9, 2015
 *      Author: wjiang2
 */

#ifndef PAIRVECTORRCPPWRAP_H_
#define PAIRVECTORRCPPWRAP_H_
#include <cytolib/compensation.hpp>
using namespace cytolib;

#include <RcppArmadillo.h> //include this instead of Rcpp.h so that RcppArmadillo inclusion won't be preceded by Rcpp.h in RcppExport.cpp
#include <RcppCommon.h>

using namespace Rcpp;

typedef std::pair<std::string, std::string> myPair;
typedef std::vector<myPair> myPairs;


namespace Rcpp {

	template <> SEXP wrap(const myPairs & kw);
}



#endif /* PAIRVECTORRCPPWRAP_H_ */
