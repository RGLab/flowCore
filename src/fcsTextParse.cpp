/*
 * readFCS.cpp
 *
 *  Created on: Feb 9, 2015
 *      Author: wjiang2
 */

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "pairVectorRcppWrap.h"
using namespace std;
//#include <Rcpp.h>
// [[Rcpp::plugins(myRegEx)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::export]]
myPairs fcsTextParse(std::string txt, bool emptyValue){

		myPairs pairs;

		/*
		 * get the first character as delimiter
		 */
		char delimiter = txt[0];

		/*
		 * check if string ends with delimiter
		 */
		bool isDelimiterEnd = txt[txt.size()-1] == delimiter;



		std::string doubleDelimiter,magicString;
		doubleDelimiter.push_back(delimiter);
		doubleDelimiter.push_back(delimiter);
		//search for the first unused odd char as replacememnt for double delimiter
		//FCS 3.1 states only 0-126 ASCII are legal delimiter, but we can't assume the file always follows the standard
		//also the TEXT main contain some special characters , thus we want to make sure the replacement char is not used anywhere in FCS TEXT
		unsigned char oddChar = 127;
		for(; oddChar < (unsigned char)256 ; oddChar++)
		{

			if(oddChar==delimiter||txt.find(oddChar)!=std::string::npos)
				continue;
			else
				break;
		}
		if(oddChar == (unsigned char) 256)
			Rcpp::stop("Can't find the unused odd character from ASCII(127-255) in FSC TEXT section!");

		std::string soddChar;
		soddChar.push_back(oddChar);
		myPair kw;
		/*
		 *	when empty value is allowed, we have to take the assumption that there is no double delimiters in any keys or values,
		 */
		if(!emptyValue)//replace the double delimiter with the odd char
			boost::replace_all(txt, doubleDelimiter, soddChar);

		std::vector<std::string> tokens;
		boost::split(tokens, txt, [delimiter](char c){return c == delimiter;});

		unsigned j = isDelimiterEnd?tokens.size()-2:tokens.size()-1;//last token, skip the last empty one when end with delimiter
		for(unsigned i = 1; i <= j; i++){//counter, start with 1 to skip the first empty tokens
			std::string token = tokens[i];
//			std::cout << token << " ";
			if(!emptyValue){
				/*
				 * restore double delimiter when needed
				 * (this slows down things quite a bit, but still a lot faster than R version,
				 *  and this double delimiter logic is not normally invoked anyway)
				 */
				boost::replace_all(token, soddChar, doubleDelimiter);
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
		 if(j%2 == 1){
			 std::string serror = "uneven number of tokens: ";
		     serror.append(boost::lexical_cast<std::string>(j));
		     Rcpp::Rcout << serror << std::endl;
			 Rcpp::Rcout << "The last keyword is dropped." << std::endl;
		 }



		return(pairs);
}
