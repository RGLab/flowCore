/*
 * pairVectorRcppWrap.hpp
 *
 *  Created on: Feb 9, 2015
 *      Author: wjiang2
 */

#ifndef PAIRVECTORRCPPWRAP_HPP_
#define PAIRVECTORRCPPWRAP_HPP_

#include <RcppCommon.h>
#include <Rcpp.h>
using namespace Rcpp;

typedef std::pair<std::string, std::string> myPair;
typedef std::vector<myPair> myPairs;


namespace Rcpp {

	template <> SEXP wrap(const myPairs & kw);
}



#endif /* PAIRVECTORRCPPWRAP_HPP_ */
