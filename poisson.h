#include "mpi.h"
#include <iostream.h>
using namespace std;

class poisson
{
private:
// Number of grid points in my processor
  int m,  n;

// Number of grid points globally
   int mg, ng;

// Number of  processors, and my processor ID.   
   int nproc, myproc;

// Mapping from global to local coordinates iglobal = mmap + ilocal   
   int mmap, nmap,mremain,nremain;

// IDs of neighbour processors, in the two dimensional processor enumeration.
   int mp, mm, np, nm;
   
// Olaped regions
   int olap;
// Tolerance and error ,iterator
   double tol;
   double error_global;
   int it;
// The arrays
   double *x;
   double *y;
   double *a;
   double *b;
   double *f;
   double *old;
   double *nnew;
   double *tempv;
// Store solutions in file
   char* file;
   int message;
// Hide standard constructor
poisson() {};

// Topology

MPI_Comm comm_cart;
int My_ID;
int np1,np2;
int p1,p2;
int dims[2];
int mycoords[2];
int periods[2];

public:
// Constructor, giving global number of points and overlap.
poisson( int i, int j, int olp,double tolerance,double xmin,double xmax,double ymin,double ymax );
void Initiate_conditions();
double * comm_boundary(double * v);
void Iterate();
void write_result(int t1);
bool check_convergence();
void write_file( char* file,double *v );
void write__a_part( char* file,double *v );
double * Jacobi(double * old);
// Destructor, ( gives back memory )
~poisson();
};