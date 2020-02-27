#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe: n - liczba funkcji a,b - granice przedzialu aproksymacji i
 * - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */
double
fi( int i, double x)
{


	return (pow(x+sqrt(x*x-1),i)+pow((x-sqrt(x*x-1)),i))/2;

}

/* Pierwsza pochodna fi */
double
dfi(int n, int i, double x)
{

	return -(i*((pow(x - sqrt(-1 + pow(x,2)),i)) - pow(x + sqrt(-1 + pow(x,2)),i)))/(2*sqrt(-1 + pow(x,2)));
}

/* Druga pochodna fi */
double
d2fi( int i, double x)
{
		return ((-1 + i) * i * pow(1 - x/sqrt(-1 + pow(x,2)),2) * pow((x - sqrt(-1 + pow(x,2))),(-2 + i)) + i * (pow(x,2)/pow((-1 + pow(x,2)),(3/2)) - 1/sqrt(-1 + pow(x,2))) * pow((x - sqrt(-1 + pow(x,2))),(-1 + i)) + (-1 + i) * i * pow((1 + x/sqrt(-1 + pow(x,2))),2) * pow((x + sqrt(-1 + pow(x,2))),(-2 + i)) + i*(-(pow(x,2)/pow((-1 + pow(x,2)),(3/2))) + 1/sqrt(-1 + pow(x,2))) * pow((x + sqrt(-1 + pow(x,2))),(-1 + i)))/2;
}

/* Trzecia pochodna fi */
double
d3fi( int i, double x)
{
		return ((-2 + i) * (-1 + i) * i * pow((1 - x/sqrt(-1 + pow(x,2))),3) *pow((x - sqrt(-1 + pow(x,2))),(-3 + i)) + 3* (-1 + i) *i* (pow(x,2)/pow((-1 + pow(x,2)),(3/2)) - 1/sqrt(-1 + pow(x,2)))* (1 - x/sqrt(-1 + pow(x,2))) *pow((x - sqrt(-1 + pow(x,2))),(-2 + i)) + i *((-3 *pow(x,3))/pow((-1 + pow(x,2)),(5/2)) + (3 *x)/pow((-1 + pow(x,2)),(3/2)))* pow((x - sqrt(-1 + pow(x,2))),(-1 + i)) + (-2 + i)* (-1 + i) *i *pow((1 + x/sqrt(-1 + pow(x,2))),3) *pow((x + sqrt(-1 + pow(x,2))),(-3 + i)) + 3 *(-1 + i)* i* (-(pow(x,2)/pow((-1 + pow(x,2)),(3/2))) + 1/sqrt(-1 + pow(x,2)))* (1 + x/sqrt(-1 + pow(x,2))) *pow((x + sqrt(-1 + pow(x,2))),(-2 + i)) + i *((3* pow(x,3))/pow((-1 + pow(x,2)),(5/2)) - (3* x)/pow((-1 + pow(x,2)),(3/2))) *pow((x + sqrt(-1 + pow(x,2))),(-1 + i)))/2;
}

/* Pomocnicza f. do rysowania bazy */
double
xfi(double a, double b, int n, int i, FILE *out)
{
	double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi         [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx      [5];
	int		j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	fprintf( out, "# nb=%d, i=%d: hi=[", n, i );
	for( j= 0; j < 5; j++ )
		fprintf( out, " %d", hi[j] );
	fprintf( out, "] hx=[" );
	for( j= 0; j < 5; j++ )
		fprintf( out, " %g", hx[j] );
	fprintf( out, "]\n" );
}

void
make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);

#ifdef DEBUG
#define TESTBASE 500
	{
		FILE           *tst = fopen("debug_base_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for( j= 0; j < nb; j++ )
			xfi( a, b, nb, j, tst );
		for (i = 0; i < TESTBASE; i++) {
			fprintf(tst, "%g", a + i * dx);
			for (j = 0; j < nb; j++) {
				fprintf(tst, " %g", fi  (a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", dfi (a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d2fi(a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d3fi(a, b, nb, j, a + i * dx));
			}
			fprintf(tst, "\n");
		}
		fclose(tst);
	}
#endif

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi(a, b, nb, i, x[k]) * fi(a, b, nb, j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(a, b, nb, j, x[k]));
	}

#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}
#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * fi  (a, b, nb, k, xx);
				spl->f1[i] += ck * dfi (a, b, nb, k, xx);
				spl->f2[i] += ck * d2fi(a, b, nb, k, xx);
				spl->f3[i] += ck * d3fi(a, b, nb, k, xx);
			}
		}
	}

#ifdef DEBUG
	{
		FILE           *tst = fopen("debug_spline_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for (i = 0; i < TESTBASE; i++) {
			double yi= 0;
			double dyi= 0;
			double d2yi= 0;
			double d3yi= 0;
			double xi= a + i * dx;
			for( k= 0; k < nb; k++ ) {
							yi += get_entry_matrix(eqs, k, nb) * fi(a, b, nb, k, xi);
							dyi += get_entry_matrix(eqs, k, nb) * dfi(a, b, nb, k, xi);
							d2yi += get_entry_matrix(eqs, k, nb) * d2fi(a, b, nb, k, xi);
							d3yi += get_entry_matrix(eqs, k, nb) * d3fi(a, b, nb, k, xi);
			}
			fprintf(tst, "%g %g %g %g %g\n", xi, yi, dyi, d2yi, d3yi );
		}
		fclose(tst);
	}
#endif

}

