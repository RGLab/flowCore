/*
 * readFCS.cpp
 *
 *  Created on: Feb 9, 2015
 *      Author: wjiang2
 */
#include <Rcpp.h>
#include <boost/regex.hpp>
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
// [[Rcpp::plugins(myRegEx)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::export]]
myPairs fcsTextParse(std::string txt, bool emptyValue){

		myPairs pairs;

		//remove trailing whitespace
//		txt = boost::regex_replace(txt,boost::regex("\s+$"), "");
		std::string div = txt.substr(0,1);

//		# regexes require double-escaping (*sigh*)
		if(div == "\\" || div == "|")
			div = "\\" + div;
//		Rcout << div << std::endl;

		boost::regex token_ex;
		unsigned i;
		myPair kw;
		if(emptyValue)
		{
//			when empty value is allowed, we have to take the assumption that there is no double delimiters in any keys or values,
			token_ex.set_expression(div);
			boost::sregex_token_iterator token_begin(txt.begin() + 1, txt.end(), token_ex, -1), token_end;
			i = 1;

			while(token_begin != token_end){

				std::string token = *token_begin++;
//				Rcout << i << ":" << std::endl;
				if((i++)%2 == 1)
					kw.first = token;//set key
				else{
					kw.second = token;//set value
					pairs.push_back(kw);//add the pair
				}


			}
		}else
		{
			 txt = txt + " "; //prepend a dummy space for token_iterater to work

//		  when empty value is not allowed, we safely parse the double delimiters as the valid values
//		  token_ex.set_expression("([^" + div + "]*)" + div ); //e.g. blahblah/

			div  =  "([^"+ div + "])" + div + "([^" + div + "])";
			token_ex.set_expression(div);
			int submatch[] = {-1,1,2};
			boost::sregex_token_iterator token_begin(txt.begin()+1,txt.end(), token_ex, submatch), token_end;

			i = 1;
			std::string prev, first, second, next;
			while(token_begin != token_end){
					prev = *token_begin++;
					if(token_begin == token_end)
							break;
					first = *token_begin++;
					second = *token_begin++;
					std::string token = next + prev + first;
					next = second;
//                      std::cout << i << ":" << std::endl;
					if((i++)%2 == 1)
						kw.first = token;//set key
//							std::cout << "K:" << token << std::endl;//set key
					else{
						kw.second = token;//set value
						pairs.push_back(kw);//add the pair
//							std::cout << "V:" << token << std::endl;//set key
						}

					}
		}





		return(pairs);
}
