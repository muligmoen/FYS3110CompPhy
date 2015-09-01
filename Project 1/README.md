# Project 1

This project aims to solve Poissons equation for a given distribution of charges. This is achieved through a matrix solution and Thomas algorithm.

## Usage
The executable "project1" should be run with the number of meshpoints N and a filename to save the data to. A python script can be called with the file returned from the program as a command line input to plot the data.

## Installation
A makefile is included for easy installation. Open a terminal and type "make" to make the executables.

## Requirements
This program requires the following libraries
* armadillo
* openblas
* lapack
* superlu

If superlu is not installed or the armadillo version is less than five, the makfile needs to be configured to not link to superlu, and "project1.cpp" needs to use the non-sparse matrice methods.