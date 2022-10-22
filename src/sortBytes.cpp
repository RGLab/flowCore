#include <cpp11.hpp>
#include "convertRawBytes.h"

/*
 * sort each element based on the byte order
 * The input is from readBin call
 */
[[cpp11::register]] cpp11::raws sortBytes(cpp11::raws bytes,
                                         cpp11::doubles byte_order) 
{

  
  int elementSize = byte_order.size();
  int nTotalBytes = bytes.size();
  //how many element to return
  int nElement = nTotalBytes / elementSize ;
  cpp11::writable::raws output(nTotalBytes);
  for(int ind = 0; ind < nElement; ind++){
    for(int i = 0; i < elementSize; i++){
      auto j = byte_order.at(i);
      
      int pos_old = ind * elementSize + i;
      int pos_new = ind * elementSize + j;
      output.at(pos_new) = bytes.at(pos_old);;
    }

  }
  return output;

}

