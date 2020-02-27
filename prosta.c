#include "makespl.h"

void
make_spl (points_t * pts, spline_t * spl)
{

  if (pts->n < 2 || alloc_spl (spl, 1)) {
    spl->n = 0;
    return;
  }
  else {
    spl->x[0] = pts->x[0];
    spl->f[0] = pts->y[0];
    spl->f1[0] =
      (pts->y[pts->n - 1] - pts->y[0]) / (pts->x[pts->n - 1] - pts->x[0]);
    spl->f2[0] = spl->f3[0] = 0;
    spl->n = 1;
  }
}
