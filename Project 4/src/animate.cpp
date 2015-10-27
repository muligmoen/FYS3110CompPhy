#include "CImg.h"
using namespace cimg_library;

#include <ctime>

#include "ising.hpp"

/*! \file animate.cpp
 * \brief Program to make an animation of the Ising system
 * 
 * This program sets up a window and prints the Ising system.
 * The program used is CImg, which is available from 
 * https://github.com/dtschump/CImg and http://cimg.eu/.
 * This program is installed as "sudo apt install cimg-dev" on
 * debian systems.
 * 
 * The system is set up randomly, and shows domains when enough
 * cycles has been completed.
 */ 

//! Main program of animate
int main() {
  const int L = 100;
  const long int seed = std::clock();
  const double beta = 1.0/2.4;
  const int wait_time = 1;
  
  Ising model(L, seed, beta);
  model.init_rand();
  
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
    
    window.wait(wait_time); // limits the framerate
  }
  
    
  
  
}
