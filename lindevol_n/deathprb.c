#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __atarist__
#  include <gfaport.h>
#else
#  include "gfaport.h"
#endif

#include "lndtypes.h"
#include "lnddispl.h"
#include "lnderror.h"
#include "lndsignl.h"
#include "lndlib.h"
#include "lndglobl.h"

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif


double death_probability(long plant_no, double p, double f_energy, double f_numcells, double leanover_penalty)
{
  long leanover, num_aerial_cells, max_leanover, i, ww2;

  if (p > 0.0)
  {
    if (f_energy != 0.0)
    {
      errno = 0;
      f_energy *= log(plant[plant_no]->cellular_energy + 1);
      if (errno)
      {
        fprintf(stderr, "death_probability: error in log(%ld [cellular_energy]): ", plant[plant_no]->cellular_energy + 1);
        perror("");
      }
    }
    if (f_numcells != 0.0)
    {
      errno = 0;
      f_numcells *= log(plant[plant_no]->num_cells);
      if (errno)
      {
        fprintf(stderr, "death_probability: error in log(%ld [num_cells]): ", plant[plant_no]->num_cells);
        perror("");
      }
    }
  }
  errno = 0;
  p *= exp(f_energy + f_numcells);
  if (errno)
  {
    fprintf(stderr, "death_probability: error in exp(%f): ", f_energy + f_numcells);
    perror("");
  }
  if ((leanover_penalty != 0.0) && (plant[plant_no]->num_cells > 1))
  {
    num_aerial_cells = 0;
    leanover = 0;
    ww2 = world_width / 2;
    for (i = 0; i < plant[plant_no]->num_cells; i++)
    {
      if (plant[plant_no]->cell[i].y >= world_soil)
      {
        num_aerial_cells++;
	leanover += ww2 - (ww2 + world_width + plant[plant_no]->cell[i].x - plant[plant_no]->cell[0].x) % world_width;
      }
    }
    if (num_aerial_cells)
    {
      max_leanover = (num_aerial_cells + 1) * num_aerial_cells / 2;
      p += (leanover_penalty * ((leanover < 0) ? -leanover : leanover)) / max_leanover;
    }
  }
  if (p > 1.0)
    p = 1.0;
  if (p < 0.0)
    p = 0.0;
  return (p);
}

