#include <Rcpp.h>
#include "convertRawBytes.h"

template <class T>
T convertRaw_impl(const BYTES & src, unsigned short thisSize, bool isBigEndian, unsigned thisStart, unsigned thisEnd)
{
  T dest;
  if(isBigEndian)
  {
    BYTES tmp(thisSize);

    for(auto k = thisStart; k <= thisEnd; k++)
      tmp.at(k % thisSize) = src.at(k);
    std::reverse(tmp.begin(),tmp.end());

    //         for(auto byte : tmp)
    //           Rcpp::Rcout << std::hex << int(byte);
    //         Rcpp::Rcout << std::endl;

    memcpy(&dest, &tmp.at(0), thisSize);
  }
  else
  {
    memcpy(&dest, &src.at(thisStart), thisSize);
  }

  return dest;
}


// [[Rcpp::plugins("cpp11")]]

/*
 * convert the raw vector to the respective type with the data that has different colSize across parameters
 * The input is from readBin call
 */
// [[Rcpp::export]]
std::vector<double> convertRawBytes( BYTES bytes, bool isInt, std::vector<unsigned short> colSize
                    , unsigned short ncol, bool isBigEndian)
{

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
  //for the simplicity, we always convert the data to double regardless of datatype
  //because we need to do this anyway when datatype uint32
  std::vector<double> output(nElement);

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
      // Rcpp::Rcout << std::to_string(thisSize) << std::endl;
      if(isInt){
        switch(thisSize)
        {
        case sizeof(BYTE)://1 byte
            output.at(ind) = (double)convertRaw_impl<BYTE>(bytes, thisSize, isBigEndian, thisStart, thisEnd);
            break;
        case sizeof(unsigned short): //2 bytes
          output.at(ind) = (double)convertRaw_impl<unsigned short>(bytes, thisSize, isBigEndian, thisStart, thisEnd);
          break;
        case sizeof(unsigned)://4 bytes
          output.at(ind) = (double)convertRaw_impl<unsigned>(bytes, thisSize, isBigEndian, thisStart, thisEnd);
          break;
        default:
          Rcpp::stop("Unsupported bitwidths when performing channel-wise reading:" + std::to_string(thisSize));
        }

      }
      else
      {
        switch(thisSize)
        {
        case sizeof(float):
          output.at(ind) = (double)convertRaw_impl<float>(bytes, thisSize, isBigEndian, thisStart, thisEnd);
          break;
        case sizeof(double):
          output.at(ind) = (double)convertRaw_impl<double>(bytes, thisSize, isBigEndian, thisStart, thisEnd);
          break;
        default:
          Rcpp::stop("Unsupported bitwidths when performing channel-wise reading:" + std::to_string(thisSize));
        }
      }


      pos = thisEnd + 1;
    }
  }
  return output;
}

