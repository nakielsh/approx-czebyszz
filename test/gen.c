#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double fun( double x ) {
				double r= ((double)rand() - RAND_MAX/2.0)/RAND_MAX/5; // +-10%
				return (1+r)*(10*x*x*x*x+5*x*x*x+x*x);
}

int main( int argc, char **argv ) {
				int n= argc > 1 ? atoi( argv[1] ) : 100;
				double a= argc > 2 ? atof( argv[2] ) : -3;
				double b= argc > 3 ? atof( argv[3] ) : 2;
				FILE *out= argc > 4 ? fopen( argv[4], "w" ) : stdout;

				srand( argc > 5 ? atoi(argv[5]) : time(NULL) );

				int i;
				double dx = (b-a)/(n-1);

				for( i =0; i < n; i++ ) {
								fprintf( out, "%g %g\n", a+i*dx, fun(a+i*dx) );
				}
				fclose( out );

				return 0;
}
