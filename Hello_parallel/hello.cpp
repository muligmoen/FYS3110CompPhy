#include <iostream>
#include <mpi.h>

int main( int argc, char *argv[] )
{
  int numprocs, my_rank;
  
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &numprocs) ;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  
  std::cout << "Hello from process " << my_rank << " of " << numprocs << " available" << std::endl;
  
  
  MPI_Finalize ();
  return 0;
}