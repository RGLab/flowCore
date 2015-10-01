/*
 * combine two uint16 (from readBin()) into one uint32 and convert it to double
 * so that the original uint32 value won't overflow R's INT_MAX 2^31-1
 *
 *  Created on: Sept 30, 2015
 *      Author: wjiang2
 */

#include <Rcpp.h>

//#include <vector>

// [[Rcpp::export]]
std::vector<double> uint2double(std::vector<unsigned> input, bool isBigEndian){
	unsigned nInput = input.size();
	unsigned nOut = nInput/2;
	std::vector<double> output(nOut);
	for(unsigned i = 0, j = 0; i < nInput - 1; i = i + 2, j++){
		unsigned left, right;
		if(isBigEndian)
		{
			left = input.at(i);
			right = input.at(i+1);


		}else
		{
			left = input.at(i+1);
			right = input.at(i);

		}
		left = left << 16;
		output.at(j) = left | right;
	}


	return (output);
}
