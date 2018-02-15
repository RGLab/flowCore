#include <Rcpp.h>
#include "convertRawBytes.h"
// [[Rcpp::plugins("cpp11")]]

/*
 * sort each element based on the byte order
 * The input is from readBin call
 */
// [[Rcpp::export]]
BYTES sortBytes(BYTES bytes, std::vector<unsigned short> byte_order)
{

  
  unsigned elementSize = byte_order.size();
  unsigned nTotalBytes = bytes.size();
  //how many element to return
  unsigned nElement = nTotalBytes / elementSize ;
  BYTES output(nTotalBytes);
  for(unsigned ind = 0; ind < nElement; ind++){
    for(unsigned i = 0; i < elementSize; i++){
      auto j = byte_order.at(i);
      
      auto pos_old = ind * elementSize + i;
      auto pos_new = ind * elementSize + j;
//       if(ind<=10)
//         Rcpp::Rcout << pos_old <<":" << pos_new << std::endl;
      output.at(pos_new) = bytes.at(pos_old);;
    }

  }
  return output;

}

