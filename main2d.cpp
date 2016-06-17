#include <iostream>
#include "poisson.h"

using namespace std;
int main(int argc, char *argv[])
{
bool converge=false;
MPI_Init(&argc,&argv);
// Rectangular domain
double xmin=0;
double xmax=2;
double ymin=0;
double ymax=1;
// Mesh size and overlap
int mg=75;
int ng=75;
int olap=2;
//Tolerance
double tolerance =.01;
poisson P(mg,ng,olap,tolerance,xmin,xmax,ymin,ymax);
//Initalize
P.Initiate_conditions();
MPI_Barrier (MPI_COMM_WORLD);
double t1=MPI_Wtime();

//Iteration
while (!converge)
{
P.Iterate();
converge=P.check_convergence();
if (converge) P.write_result(t1);
}

MPI_Finalize();
}
