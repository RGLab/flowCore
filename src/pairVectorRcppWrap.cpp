/*
 * pairVectorRcppWrap.cpp
 *
 *  Created on: Feb 9, 2015
 *      Author: wjiang2
 */
#include <Rcpp.h>
#include "pairVectorRcppWrap.hpp"

namespace Rcpp {
	template <> SEXP wrap(const myPairs & kw){
		unsigned nSize = kw.size();
		Rcpp::CharacterVector res(nSize);
		Rcpp::CharacterVector res_names(nSize);
		for(unsigned i = 0; i < nSize; i++){
			myPair thisKw = kw.at(i);
			res[i] = thisKw.second;
			res_names[i] = thisKw.first;
		}
		res.names() = res_names;
		return res;
	}
}





