#include "poisson.h"
#include <math.h>
#include <iostream>

poisson::poisson(int i,int j,int olp,double tolerance,double xmin,double xmax,double ymin,double ymax)
{
  int s1,s2;
  mg=i;
  ng=j;
  olap=olp;
  tol=tolerance;
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);   
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
  // Define Cartesian Toplogy
  for (int z=0;z<2;z++) {dims[z]=0;periods[z]=0;}
  MPI_Dims_create(nproc,2,dims);np1=dims[0];np2=dims[1];
  MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,0,&comm_cart); 
  MPI_Cart_coords (comm_cart,myproc,2,mycoords);
  p1=mycoords[0];p2=mycoords[1];
  MPI_Cart_rank(comm_cart,mycoords,&My_ID);
  MPI_Cart_shift(comm_cart, 1, 1, &nm, &np);  MPI_Cart_shift(comm_cart, 0, 1, &mm, &mp);
  
  //Domain decomposition
  s1=0; s2=0;
  s1=  (mg+(np1-1)*olap)/(np1);     m=s1;
  s2=  (ng+(np2-1)*olap)/(np2);     n=s2;
  mremain=(mg+(np1-1)*olap)  % np1;
  nremain=(ng+(np2-1)*olap)  % np2;
  mmap=(p1)*(s1-olap)+min(p1,mremain);
  nmap=(p2)*(s2-olap)+min(p2,nremain);
  if (p1+1<=mremain) {m= m+1;}
  if (p2+1<=nremain) {n= n+1;}
  // Allocate memory for variables
  a=new double[n*m];
  b=new double[n*m];
  f=new double[n*m];
  nnew=new double[n*m];
  old=new double[n*m];
  x=new double[ng*mg];
  y=new double[ng*mg];
  // Mesh geberation
   for( j= 1 ; j<=ng ; j++ )
   for( i= 1 ; i<=mg ; i++ )
   {
   x[i-1 + mg* (j-1)]=xmin+(i-1)*((xmax-xmin)/(mg-1));
   y[i-1 + ng* (j-1)]=ymin+(j-1)*((ymax-ymin)/(mg-1));
   }
}

void poisson::Initiate_conditions()
{

// Initialtion
int i,j;
 for( j= 1 ; j<=n ; j++ )

   for( i= 1 ; i<=m ; i++ )
   {
   f[i-1 + m* (j-1)]=x[i-1+mmap + mg* (j-1+nmap)]*exp(y[i-1+mmap + ng* (j-1+nmap)]);
   a[i-1 + m* (j-1)]=0;
   b[i-1 + m* (j-1)]=0;
   }

 // Handling boundary conditions
 if (mp==-1)
 {
for( j= 1 ; j<=n ; j++ )
   {
     a[m-1+m*(j-1)]=2*exp(y[mg-1 + ng* (j-1+nmap)]);
     b[m-1+m*(j-1)]=2*exp(y[mg-1 + ng* (j-1+nmap)]);
   }
}

if (nm==-1)
 {
for( i= 1 ; i<=m ; i++ )
   {
     a[i-1+m*(0)]=x[i-1+mmap + mg* (0)];
     b[i-1+m*(0)]=x[i-1+mmap + mg* (0)];
   }
}

if (np==-1)
 {
for( i= 1 ; i<=m ; i++ )
   {
     a[i-1+m*(n-1)]=exp(1)*x[i-1+mmap + mg* (ng-1)];
     b[i-1+m*(n-1)]=exp(1)*x[i-1+mmap + mg* (ng-1)];
   }
}
}

void poisson::Iterate()
{
a=comm_boundary(a);MPI_Barrier (comm_cart);b= Jacobi(a);
b=comm_boundary(b);MPI_Barrier (comm_cart);a= Jacobi(b);
}

double * poisson::Jacobi(double * old)
{
// Jacobi Iteration
int i,j;
double h=1/(m);
it=it+1;
 // Handling boundary conditions
for( j= 1 ; j<=n ; j++ ) for( i= 1 ; i<=m ; i++ ){ nnew[i-1 + m* (j-1)]=0;}
if (mm==-1){for( j= 1 ; j<=n ; j++ ){nnew[0+m*(j-1)]=1;  nnew[0+m*(j-1)]=1;}}
if (nm==-1){for( i= 1 ; i<=m ; i++ ) {  nnew[i-1+m*(0)]=1; nnew[i-1+m*(0)]=1;}}
for( j= 0 ; j<=n-1 ; j++ ){for( i= 0 ; i<=m-1 ; i++ ){nnew[i+(m)*(j)]=old[i+(m)*(j)];}}

// Update solution
for( j= 1 ; j<=n-2 ; j++ ){for( i= 1 ; i<=m-2 ; i++ )
  {nnew[i+(m)*(j)]=.25*( old[i-1+(m)*(j)]+old[i+(m)*(j+1)]
                        +old[i+(m)*(j-1)]+old[i+1+(m)*(j)])
                        -h*h*f[i+(m)*(j)];}
}
return nnew;
}

double * poisson::comm_boundary(double * v)
{
// // Boundary communication for non-overlapping decomposition
MPI_Status status;
int i,j;
double *vv = new double[m*n];
tempv = new double[m*n+2];
// -------------------   copy internal data  --------------------------------
for( int j= 1 ; j<=n ; j++ ){for( int i= 1 ; i<=m ; i++ )
{tempv[i-1 + (m)*(j-1)+2]=v[i-1 + m*(j-1)];}
}

// -------------------   exchange ghost points in     x direction  ----------
for( j= 1 ; j<=n ; j++ ) {
MPI_Sendrecv(v+1+m*(j-1), 1,MPI_DOUBLE,mm,100,
tempv+m-1+(m)*(j-1)+2,1,MPI_DOUBLE,mp,100,comm_cart,&status);}

for( j= 1 ; j<=n ; j++ ) {
MPI_Sendrecv(v+m-2+m*(j-1), 1,MPI_DOUBLE,mp,100
,tempv+0+(m)*(j-1)+2,1,MPI_DOUBLE,mm,100,comm_cart,&status);}

// -------------------   exchange ghost points in     y direction  ------------------
for( i= 1 ; i<=m ; i++ ){
MPI_Sendrecv(v+i-1+m*(1), 1,MPI_DOUBLE,nm,100,
tempv+i-1+(m)*(n-1)+2,1,MPI_DOUBLE,np,100,comm_cart,&status);}

for( i= 1 ; i<=m ; i++ )
{MPI_Sendrecv(v+i-1 + m* (n-2),1,MPI_DOUBLE,np,100,
tempv+i-1+(m)*(0)+2,1,MPI_DOUBLE,nm,100,comm_cart,&status );};
for( int j= 1 ; j<=n ; j++ ){for( int i= 1 ; i<=m ; i++ )
{vv[i-1 + (m)*(j-1)]=tempv[i-1 + m*(j-1)+2];
}}
return vv;
}
bool poisson::check_convergence()
{
  bool converge;
  int i,j ;
  double S=0;
  double temp=0;
  for( j= 1 ; j<=n-1 ; j++ ){
  for( i= 1 ; i<=m-1 ; i++ ){
  temp=a[i-1+(m)*(j-1)] -  b[i-1+(m)*(j-1)];
  S=S+temp*temp;
  } }
  MPI_Allreduce(&S, &error_global, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  converge=(error_global<tol);
  return converge;
}

void poisson::write_file( char* file,double *v )
{
if( myproc == 0 )
  {
    //     Open file and write dimensions, mg, ng 
    int nr;
    FILE *fp;
    fp = fopen (file , "w" );
    nr  = fwrite(&mg, sizeof(int ), 1, fp );
    nr  = fwrite(&ng, sizeof(int ), 1, fp );
    fclose(fp);
  }

 if( myproc == 0 )
  { 
    //  write_my_part
    int position,nr,i,j,ind;
    FILE *fp ;
     fp=fopen(file,"r+");
      for( int j=1 ; j<=n; j++ )
	{
	  position=mmap+(j-1+nmap)*mg+1;
	  fseek( fp, position*sizeof(double), SEEK_SET );
	  nr  = fwrite( v + m*(j-1), sizeof(double), m, fp );
	}
    fclose(fp);
  }

else
  {
    //  receive a message from myproc-1
    MPI_Status status;
    MPI_Recv(&message,1,MPI_INT,myproc-1,99,MPI_COMM_WORLD,&status);
    //   write_my_part
     write__a_part(file,v);
  }
if( myproc < nproc-1 )
  {
    message=1;
    //  send message to myproc + 1
    MPI_Send(&message,1,MPI_INT,myproc+1,99,MPI_COMM_WORLD);
  }
}

void poisson::write__a_part( char* file,double *v )
{
// Open file
int position,nr,i,j,ind;
FILE *fp ;
 fp=fopen(file,"r+");
 for( int j=1 ; j<=n; j++ )
  {
    position=(mmap)+(j-1+nmap)*mg+1;
    fseek( fp, position*sizeof(double), SEEK_SET );
    nr  = fwrite( v + m*(j-1), sizeof(double), m, fp );
  } 
fclose(fp);
}
void poisson::write_result(int t1)
{
 int i,j,ctmin,ctmax;
 double rtmin,rtmax,t2;
 write_file("solution.dat",b );
  t2=MPI_Wtime()-t1;
  MPI_Reduce(&t2,&rtmax,1,MPI_DOUBLE,MPI_MAX,0,comm_cart);
  MPI_Reduce(&t2,&rtmin,1,MPI_DOUBLE,MPI_MIN,0,comm_cart);
 if(myproc ==0)
  {
  int i,j ;
  cout << "Solution converged after " << it <<" iterations. "<< endl;
  cout << "Solution is saved in file : solution.dat . "<< endl;
  cout << "Real time (s) : " << rtmin<< "  " << rtmax << endl;
  }

  }

poisson::~poisson()
{
delete[] a;
delete[] b;
delete[] f;
delete[] nnew;
delete[] old;
}

