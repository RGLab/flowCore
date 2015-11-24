#include <Rcpp.h>
#include "convertRawBytes.h"



// [[Rcpp::plugins("cpp11")]]

/*
 * convert the raw vector to the respective type with the data that has different colSize across parameters
 * The input is from readBin call
 */
template <class T>
void convertRawBytes_impl( BYTES bytes,  std::vector<T> & output, std::vector<unsigned short> colSize
                    , unsigned short ncol, unsigned nrow, bool isBigEndian) 
  {
  
  
  
  auto pos = 0;
  // Rcpp::Rcout << "event=" << nrow << " ncol=" << ncol << std::endl;
  
  for(auto i = 0; i < nrow; i++){

    for(auto j = 0; j < ncol; j++){
      //convert each element
      auto thisStart = pos;
      auto thisSize = colSize.at(j);
      
      if(thisSize != 1 && thisSize != 2 && thisSize !=  4 && thisSize !=  8)
        Rcpp::stop("Multiple different odd bitwidths are not supported!");
      
      auto thisEnd = pos + thisSize - 1;
      auto ind = i * ncol + j;
      
      // Rcpp::Rcout << "start="  << thisStart << " " << "end=" << thisEnd << " " << "index=" << ind << std::endl;
      // Rcpp::Rcout << "bytes" << std::endl;
      
//       for(auto k = thisStart; k <= thisEnd; k++)
//         Rcpp::Rcout << std::hex << int(bytes.at(k));
//       Rcpp::Rcout << std::endl;
        
      if(isBigEndian)
      {
        BYTES tmp(thisSize);
        
        for(auto k = thisStart; k <= thisEnd; k++)
            tmp.at(k % thisSize) = bytes.at(k);
        std::reverse(tmp.begin(),tmp.end());
        
//         for(auto byte : tmp)
//           Rcpp::Rcout << std::hex << int(byte);
//         Rcpp::Rcout << std::endl;
        
        memcpy(&output.at(ind), &tmp.at(0), thisSize);
      }
      else
      {
        memcpy(&output.at(ind), &bytes.at(thisStart), thisSize);
      }
        
      pos = thisEnd + 1;	          
    } 
  }
  
}




// [[Rcpp::export]]
SEXP convertRawBytes(BYTES bytes, bool isInt, std::vector<unsigned short> colSize
                , unsigned short ncol, bool isBigEndian){
  
  // RCPP_RETURN_VECTOR does not work here since it only takes one function argument
  if(colSize.size()!=ncol)
    Rcpp::stop("The length of 'colSize' vector does not match to 'ncol'!");
  
  unsigned nBytes = bytes.size();
  //total bytes for each row
  unsigned short nRowSize = 0;
  for(auto size : colSize)
    nRowSize+=size;
  //how many rows
  unsigned nrow = nBytes/nRowSize;
  //how many element to return
  auto nElement = nrow * ncol;
  
  
  if(isInt)
  {
    std::vector<int> output(nElement);
    convertRawBytes_impl<int>(bytes, output, colSize, ncol,nrow, isBigEndian) ;
    return Rcpp::wrap(output);
  }
  else
  {
    //we only support the same byte size for numeric type
    switch(colSize.at(0))
    {
    case sizeof(float):
        {
            std::vector<float> output(nElement);
            convertRawBytes_impl<float>(bytes, output, colSize, ncol,nrow, isBigEndian);
            return Rcpp::wrap(output);
        }   
    case sizeof(double):
        {
          std::vector<double> output(nElement);
          convertRawBytes_impl<double>(bytes, output, colSize, ncol, nrow,isBigEndian);    
          return Rcpp::wrap(output);
        }
      
    }
    
  }
    
  
  
}
