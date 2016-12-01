/*
 * readFCS.cpp
 *
 *  Created on: Feb 9, 2015
 *      Author: wjiang2
 */

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include "pairVectorRcppWrap.h"
//#include <Rcpp.h>
// [[Rcpp::plugins(myRegEx)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::export]]
myPairs fcsTextParse(std::string txt, bool emptyValue){

		myPairs pairs;

		/*
		 * get the first character as delimiter
		 */
		std::string delimiter = txt.substr(0,1);

		/*
		 * check if string ends with delimiter
		 */
		bool isDelimiterEnd = txt.substr(txt.size()-1, 1) == delimiter;

//		regexes require double-escaping (*sigh*)
//		if(delimiter == "\\" || delimiter == "|")
			delimiter = "\\" + delimiter;



		std::string doubleDelimiter,magicString;
		doubleDelimiter = delimiter + delimiter;
		magicString = "\\0QuickAndDirty\\0";
//		std::cout << doubleDelimiter << ":" << magicString <<std::endl;
		unsigned i = 0; //counter
		myPair kw;
		/*
		 *	when empty value is allowed, we have to take the assumption that there is no double delimiters in any keys or values,
		 */
		if(!emptyValue)//replace the double delimiter with a magic strings
			txt = boost::regex_replace(txt, boost::regex(doubleDelimiter), magicString);//somehow boost::replace_all won't do the job for \\\\
		std::cout << txt << std::endl;

		/*
		 * then split by single delimiter
		 */
		boost::sregex_token_iterator token_begin(txt.begin() + 1, txt.end(), boost::regex(delimiter), -1), token_end;
		while(token_begin != token_end){
			i++;
			std::string token = *token_begin++;
//			std::cout << token << " ";
			if(!emptyValue){
				/*
				 * restore double delimiter when needed
				 * (this slows down things quite a bit, but still a lot faster than R version,
				 *  and this double delimiter logic is not normally invoked anyway)
				 */
				token = boost::regex_replace(token, boost::regex(magicString), doubleDelimiter);
//				std::cout << token;
			}
//			std::cout << std::endl;

			if((i)%2 == 1)
			{
				if(token.empty())
					// Rcpp::stop (temporarily switch from stop to range_error due to a bug in Rcpp 0.12.8)
					throw std::range_error("Empty keyword name detected!If it is due to the double delimiters in keyword value, please set emptyValue to FALSE and try again!");
				kw.first = token;//set key
			}
			else{
				kw.second = token;//set value
				pairs.push_back(kw);//add the pair
			}


		}

		/*
		 * check if kw and value are paired
		 */
		 if(i%2 == 1){
			 if(isDelimiterEnd){
			   // Rcpp::stop
			   std::string serror = "uneven number of tokens: ";
			   serror.append(boost::lexical_cast<std::string>(i-1));
			   throw std::range_error(serror.c_str());
			 }	 
			 else
				 Rcpp::Rcout << "the text section does not end with delimiter: " << delimiter << ". The last keyword is dropped." << std::endl;;
		 }



		return(pairs);
}
