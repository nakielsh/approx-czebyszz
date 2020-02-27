#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>

void
make_spl (points_t * pts, spline_t * spl)
{

	int n= pts->n - 1;
	matrix_t *eqs= make_matrix( n*3, n*3+1 );
	double *x = pts->x;
	double *y = pts->y;

	int i;

	for( i= 0; i < n; i++ ) {
		double dx= x[i+1] - x[i];
		int if1= 3*i;
		int if2= if1+1;
		int if3= if2+1;
		put_entry_matrix( eqs, if1, if1, dx );
		put_entry_matrix( eqs, if1, if2, dx*dx/2 );
		put_entry_matrix( eqs, if1, if3, dx*dx*dx/6 );
		put_entry_matrix( eqs, if1, n*3, y[i+1]-y[i] );
		put_entry_matrix( eqs, if2, if1, 1 );
		put_entry_matrix( eqs, if2, if2, dx );
		put_entry_matrix( eqs, if2, if3, dx*dx/2 );
		if( if3+1 < n*3 )
			put_entry_matrix( eqs, if2, if3+1, -1 );
		else
			put_entry_matrix( eqs, if2, if1, 0 );
		put_entry_matrix( eqs, if3, if2, 1 );
		put_entry_matrix( eqs, if3, if3, dx );
		if( if3+2 < n*3 )
			put_entry_matrix( eqs, if3, if3+2, -1 );
	}

#ifdef DEBUG
	write_matrix( eqs, stdout );
#endif

	if( piv_ge_solver( eqs ) ) {
		spl->n = 0;
		return;
	}

#ifdef DEBUG
	write_matrix( eqs, stdout );
#endif

  if ( alloc_spl (spl, pts->n) == 0 ) {
    spl->n = pts->n;
		for( i= 0; i < n; i++ ) {
			spl->x[i]= pts->x[i];
			spl->f[i]= pts->y[i];
			spl->f1[i]= get_entry_matrix( eqs, 3*i,   3*n );
			spl->f2[i]= get_entry_matrix( eqs, 3*i+1, 3*n );
			spl->f3[i]= get_entry_matrix( eqs, 3*i+2, 3*n );
		}
		spl->x[n]= pts->x[n];
		spl->f[n]= pts->y[n];
		spl->f1[n]= spl->f1[n-1];
		spl->f2[n]= 0;
		spl->f3[n]= 0;
  }
}
