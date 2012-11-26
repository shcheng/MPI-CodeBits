#include <iostream>
#include <cmath>
#include "mpi.h"

using namespace std;

const double TT = acos(-1.);

double f(double x);
double trapInteg(double local_a, double local_b, int local_n, double global_h);

int main(int argc, char * argv[]) {

  // Computation-related variables
  int global_n = 1000;
  double global_a = -1, global_b = 1.;
  double global_h = (global_b-global_a)/global_n;
  double integral;
  int local_n;
  double local_a, local_b;
  double local_integral;

  // Parallelization-related variables
  int rank;           // this process
  int nproc;          // number of processes
  int dest=0;   // message destination
  int tag=0;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  local_n = global_n / nproc;
  local_a = global_a + (rank * local_n * global_h);
  local_b = local_a  + (local_n * global_h);
  local_integral = trapInteg(local_a, local_b, local_n, global_h);

  if(rank!=0) {
    MPI_Send(&local_integral, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
  } else {
    integral = local_integral;  // This probably needs to go before the following line
    for(int orig=1; orig<nproc; orig++) {
      MPI_Recv(&local_integral, 1, MPI_DOUBLE, orig, tag, MPI_COMM_WORLD, &status);
      integral += local_integral;
    }
  }

  // Print out the result
  if(rank==0) {
    cout << "With (a,b,n)=(" << global_a << "," << global_b << "," << global_n << ")," << endl;
    cout << "S(a,b,n){f(x)} = " << integral << endl;
  }

  MPI_Finalize();
  return 0;
}



double f(double x) {
  return exp(-x*x/2.)/sqrt(2.*TT);
}

double trapInteg(double local_a, double local_b, int local_n, double global_h) {
  double pIntegral = (f(local_a)+f(local_b))/2.;
  for(int i=1; i<=(local_n-1); i++) pIntegral += f(local_a+(double)i*global_h);
  return pIntegral*global_h;
}
