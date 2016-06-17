
#include"iostream"
#include<cmath>
#include"fstream"
#include <stdio.h>
using namespace std;
int main()
{ 
string Extension=".dat",filnam;
string filnamn="sol";
filnamn.append(Extension);
ofstream skriv(filnamn.c_str());
FILE * fp=fopen("solution.dat","r");
double * buffer;
int ind,nr,position,i,j,t1,t2, mg,ng,m1,n1;
nr =fread (&mg, sizeof(int) ,1, fp);
nr= fread (&ng, sizeof(int) ,1, fp);
buffer=new double[mg*ng];

// writing solution
for(  j=1 ; j<=ng; j++ )
    {
    position=(j-1)*mg+1;
    fseek( fp, position*sizeof(double), SEEK_SET );
    nr  = fread( buffer+mg*(j-1) , sizeof(double), mg, fp );
    }
fclose (fp);

 for ( i=0 ; i<mg ;i++)
 { for (j=1;j<ng+1;j++)
    { ind = i + mg*(j-1);   skriv  <<buffer [ind]<< "\n"; } }
return 0;
}


