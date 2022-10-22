// Copyright (c) 2021 Ozette Technologies
#include "convertRawBytes.h"
#include <string>
#include "cpp11.hpp"
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


    memcpy(&dest, &tmp.at(0), thisSize);
  }
  else
  {
    memcpy(&dest, &src.at(thisStart), thisSize);
  }

  return dest;
}



/*
 * convert the raw vector to the respective type with the data that has different colSize across parameters
 * The input is from readBin call
 */
[[cpp11::register]] cpp11::sexp convertRawBytes(
   std::vector<unsigned char> bytes, bool isInt, cpp11::integers colSize,
   int ncol, bool isBigEndian) {
 if (static_cast<std::size_t>(colSize.size()) !=
     static_cast<std::size_t>(ncol))
    throw std::range_error("The length of 'colSize' vector does not match to 'ncol'!");

  int nBytes = bytes.size();
  //total bytes for each row
  int nRowSize = 0;
  for(auto size : colSize)
    nRowSize+=size;
  //how many rows
  int nrow = nBytes/nRowSize;
  //how many element to return
  int nElement = nrow * ncol;
  //for the simplicity, we always convert the data to double regardless of datatype
  //because we need to do this anyway when datatype uint32
  cpp11::writable::doubles output(nElement);

  int pos = 0;

  for(int i = 0; i < nrow; i++){

    for(int j = 0; j < ncol; j++){
      //convert each element
      int thisStart = pos;
      int thisSize = colSize.at(j);

      if(thisSize != 1 && thisSize != 2 && thisSize !=  4 && thisSize !=  8)
        throw std::range_error("Multiple different odd bitwidths are not supported!");

      int thisEnd = pos + thisSize - 1;
      int ind = i * ncol + j;

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
        case sizeof(uint64_t)://8 bytes
          output.at(ind) = (double)convertRaw_impl<uint64_t>(bytes, thisSize, isBigEndian, thisStart, thisEnd);
          break;
        default:
          std::string serror = "Unsupported bitwidths when performing channel-wise reading:";
          serror.append(std::to_string(thisSize));
          throw std::range_error(serror.c_str());
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
          std::string serror ="Unsupported bitwidths when performing channel-wise reading:";
          serror.append(std::to_string(thisSize));
          throw std::range_error(serror.c_str());
        }
      }


      pos = thisEnd + 1;
    }
  }
  return output;
}

