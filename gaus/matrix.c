#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

matrix_t *
make_matrix (int rn, int cn)
{
  matrix_t *new_mat = malloc (sizeof *new_mat);
  if (new_mat == NULL)
    return NULL;
  if ((new_mat->e =
       malloc ((size_t) rn * (size_t) cn * sizeof *new_mat->e)) == NULL) {
    free (new_mat);
    return NULL;
  }
  new_mat->rn = rn;
  new_mat->cn = cn;
  memset (new_mat->e, 0, (size_t) (rn * (size_t) cn * sizeof *new_mat->e));
  return new_mat;
}

void
free_matrix (matrix_t * m)
{
  free (m->e);
  free (m);
}

void
put_entry_matrix (matrix_t * m, int i, int j, double val)
{
  if (i >= 0 && i < m->rn && j >= 0 && j <= m->cn)
    m->e[i * m->cn + j] = val;
}

void
add_to_entry_matrix (matrix_t * m, int i, int j, double val)
{
  if (i >= 0 && i < m->rn && j >= 0 && j <= m->cn)
    m->e[i * m->cn + j] += val;
}

double
get_entry_matrix (matrix_t * m, int i, int j )
{
  if (i >= 0 && i < m->rn && j >= 0 && j <= m->cn)
    return m->e[i * m->cn + j];
	else
		return -999;
}

matrix_t *
read_matrix (FILE * in)
{
  int rn, cn;
  int i, j;
  matrix_t *new_mat;
  if (fscanf (in, "%d %d", &rn, &cn) != 2)
    return NULL;

  if ((new_mat = make_matrix (rn, cn)) == NULL)
    return NULL;
  for (i = 0; i < rn; i++)
    for (j = 0; j < cn; j++)
      if (fscanf (in, "%lf", &new_mat->e[i * cn + j]) != 1) {
        free_matrix (new_mat);
        return NULL;
      }

  return new_mat;
}

void
write_matrix (matrix_t * m, FILE * out)
{
  int i, j;
  if (m == NULL) {
    fprintf (out, "Matrix is NULL\n");
    return;
  }

  fprintf (out, "%d %d\n", m->rn, m->cn);
  for (i = 0; i < m->rn; i++) {
    for (j = 0; j < m->cn - 1; j++)
      fprintf (out, "%8.5f ", m->e[i * m->cn + j]);
    fprintf (out, "%8.5f\n", m->e[i * m->cn + j]);
  }
}

matrix_t *
copy_matrix (matrix_t * s)
{
  matrix_t *d = NULL;
  if (s != NULL)
    d = make_matrix (s->rn, s->cn);
  if (d != NULL) {
    memcpy (d->e, s->e, s->rn * s->cn * sizeof *s->e);
    /* int i;
       for( i= 0; i < s->rn*s->cn; i++ )
       *(d->e+i)= *(s->e+i);
     */
    /* d->rn= s->rn; d->cn= s->cn; done in make_matrix */
  }
  return d;
}

matrix_t *
transpose_matrix (matrix_t * s)
{
  matrix_t *d = NULL;
  if (s != NULL)
    d = make_matrix (s->rn, s->cn);
  if (d != NULL) {
    int i, j;
    for (i = 0; i < s->rn; i++)
      for (j = 0; j < s->cn; j++)
        *(d->e + j * d->cn + i) = *(s->e + i * s->cn + j);
    /* d->rn= s->cn; d->cn= s->rn; done in make_matrix */
  }
  return d;
}

void
xchg_rows (matrix_t * m, int i, int j)
{
  if (m != NULL && i >= 0 && i < m->rn && j >= 0 && j < m->rn) {
    int k;
    double tmp;
    for (k = 0; k < m->cn; k++) {
      tmp = *(m->e + i * m->cn + k);
      *(m->e + i * m->cn + k) = *(m->e + j * m->cn + k);
      *(m->e + j * m->cn + k) = tmp;
    }
  }
}

void
xchg_cols (matrix_t * m, int i, int j)
{
  if (m != NULL && i >= 0 && i < m->cn && j >= 0 && j < m->cn) {
    int k;
    double tmp;
    for (k = 0; k < m->rn; k++) {
      tmp = *(m->e + k * m->cn + i);
      *(m->e + k * m->cn + i) = *(m->e + k * m->cn + j);
      *(m->e + k * m->cn + j) = tmp;
    }
  }
}

matrix_t *
mull_matrix (matrix_t * a, matrix_t * b)
{
  if (a == NULL || b == NULL || a->cn != b->rn)
    return NULL;
  else {
    matrix_t *c = make_matrix (a->rn, b->cn);
    int i, j, k;
    if (c == NULL)
      return NULL;
    for (i = 0; i < c->rn; i++) {
      for (j = 0; j < c->cn; j++) {
        double s = 0;
        for (k = 0; k < a->cn; k++)
          s += *(a->e + i * a->cn + k) * *(b->e + k * b->cn + j);
        *(c->e + i * c->cn + j) = s;
      }
    }
    return c;
  }
}

matrix_t *
ge_matrix (matrix_t * a)
{
  matrix_t *c = copy_matrix (a);
  if (c != NULL) {
    int i, j, k;
    int cn = c->cn;
    int rn = c->rn;
    double *e = c->e;
    for (k = 0; k < rn - 1; k++) {      /* eliminujemy (zerujemy) kolumnę nr k */
      for (i = k + 1; i < rn; i++) {    /* pętla po kolejnych
                                           wierszach poniżej diagonalii k,k */
        double d = *(e + i * cn + k) / *(e + k * cn + k);
        for (j = k; j < cn; j++)
          *(e + i * cn + j) -= d * *(e + k * cn + j);
      }
    }
  }
  return c;
}

int
bs_matrix (matrix_t * a)
{
  if (a != NULL) {
    int r, c, k;
    int cn = a->cn;
    int rn = a->rn;
    double *e = a->e;

    for (k = rn; k < cn; k++) { /* pętla po prawych stronach */
      for (r = rn - 1; r >= 0; r--) {   /* petla po niewiadomych */
        double rhs = *(e + r * cn + k); /* wartość prawej strony */
        for (c = rn - 1; c > r; c--)    /* odejmujemy wartości już wyznaczonych niewiadomych *
                                           pomnożone przed odpowiednie współczynniki */
          rhs -= *(e + r * cn + c) * *(e + c * cn + k);
        *(e + r * cn + k) = rhs / *(e + r * cn + r);    /* nowa wartość to prawa strona / el. diagonalny */
      }
    }
    return 0;
  }
  else
    return 1;
}
