// Copyright (c) 2021 Ozette Technologies
/*
 * readFCS.cpp
 *
 *  Created on: Feb 9, 2015
 *      Author: wjiang2
 */

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "cpp11.hpp"
#include <iostream>
#include <cytolib/compensation.hpp>
using namespace std;
using namespace cytolib;
[[cpp11::register]] cpp11::writable::doubles_matrix<> string_to_spill(string key){
	compensation comp(key);
 	arma::mat spillover = comp.get_spillover_mat();
cpp11::writable::doubles_matrix<> res(spillover.n_rows, spillover.n_cols);
  // copy spillover matrix.
  for (auto j = 0; j < spillover.n_cols; j++) {
    for (auto i = 0; i < spillover.n_rows; i++) {
      res(i, j) = spillover(i, j);
    }
  }
  
  cpp11::writable::strings markers(comp.marker);
  cpp11::writable::list_of<cpp11::writable::strings> mydims(
      {R_NilValue, markers});
  Rf_setAttrib(cpp11::as_sexp(res), cpp11::as_sexp({"dimnames"}),
               cpp11::as_sexp(mydims));
	return res;
}
[[cpp11::register]] std::string spill_to_string(
    cpp11::doubles_matrix<> rmat, std::vector<std::string> markers) {
  arma::Mat<double> mat(rmat.nrow(), rmat.ncol());
  for (auto j = 0; j < rmat.ncol(); j++) {
    for (auto i = 0; i < rmat.nrow(); i++) {
      mat(i, j) = rmat(i, j);
    }
  }
	compensation comp(mat, markers);
	return comp.to_string();

}

[[cpp11::register]] cpp11::sexp fcsTextParse(std::string txt, bool emptyValue){
typedef std::pair<std::string, std::string> myPair;
typedef std::vector<myPair> myPairs;

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
			cpp11::stop("Can't find the unused odd character from ASCII(127-255) in FSC TEXT section!");

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
				boost::replace_all(token, soddChar, string(1, delimiter));//unescape the double delimiter to single one
//				std::cout << token;
			}
//			std::cout << std::endl;

			if((i)%2 == 1)
			{
				if(token.empty())
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
		     serror.append("\n");
  			  Rprintf(serror.data());
			    Rprintf("The last keyword is dropped.\n");

		 }



		// Convert pairs to a named vector, which is what R expects.
  cpp11::writable::strings returned_named_vector(pairs.size());
  std::vector<std::string> keys(pairs.size());
  try {
    for (int i = 0; i < pairs.size(); i++) {
      keys[i] = pairs[i].first;
      returned_named_vector[i] = pairs[i].second;
    }
    returned_named_vector.names() = keys;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
  }
  return (returned_named_vector);

}
