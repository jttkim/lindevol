      #include "lndvals.h"
#include "lndtypes.h"
#include "lndglobl.h"
#include "lnderror.h"
#include "lndlibin.h"
#include "lndlib.h"


int plant_boundingbox(const PLANT *plant, long *xmin, long *xmax, long *ymax)
{
  long i, dx, ww2 = world_width / 2;

  *xmin = 0;
  *xmax = 0;
  *ymax = 0;
  for (i = 0; i < plant->num_cells; i++)
  {
    dx = (plant->cell[i].x - plant->cell[0].x + world_width) % world_width;
    if (dx >= ww2)
      dx -= world_width;
    if (dx < *xmin)
      *xmin = dx;
    if (dx > *xmax)
      *xmax = dx;
    if (plant->cell[i].y > *ymax)
      *ymax = plant->cell[i].y;
  }
  return (0);
}

