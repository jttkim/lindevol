/* $Id: pbbox.c,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: pbbox.c,v $
 * Revision 1.2  2000/01/30 03:11:00  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

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

