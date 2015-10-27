#include "CImg.h"
using namespace cimg_library;
#include <iostream>
#include <ctime>

#include "ising.hpp"



int main() {
  const int L = 100;
  const long int seed = std::clock();
  const double beta = 4;
  
  Ising model(L, seed, beta);
  
  const int width = 1;
  
  // Getting the future buffer
  lat_t* image_buffer = model.buffer();

  //making a window
  CImgDisplay window(900,600,"Ising model",0);
  
  //binding a buffer to image
  CImg<lat_t> image(image_buffer,L,L,1,width,true);
  
  
  while (!window.is_closed()) { // event binding
    if (window.is_resized()) {
      window.resize();
    }
    window.display(image);
    
    
    // updating buffer live here
    model.try_flip();
    
    window.wait(5); // limits the framerate
  }
  
    
  
  
}
