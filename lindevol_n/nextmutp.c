/* $Id: nextmutp.c,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: nextmutp.c,v $
 * Revision 1.2  2000/01/30 03:11:00  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

#include "genomelib.h"
#include "lndtypes.h"
#include "lndlib.h"
#include "lnderror.h"


#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif


long next_mutpos(double m, long pos)
{
  double p;
  double r, p1, lr, lm;

  if (m >= 1.0 - DBL_EPSILON)
    return (pos);
  if (m < DBL_EPSILON)
    m = DBL_EPSILON;
  r = lnd_rnd();
  errno = 0;
  lr = log(1.0 - r);
  lm = log(1.0 - m);
  p1 = lr / lm;
  if (errno)
  {
    printf("next_mutpos: error with log() %d, pos=%lu, m=%f, p1=%f, r=%f\n", errno, pos, m, p1, r);
    if (errno != ERANGE)
    {
      perror("next_mutpos failed");
      /* printf("errno = %d\n", errno); */
    }
    p = DBL_MAX;
  }
  errno = 0;
  p = floor(pos + p1);

#ifdef dreck

  if (errno)
  {
    printf("next_mutpos: error with floor %d, pos=%lu, p1=%f, p=%f, r=%f\n", errno, pos, p1, p, r);
    if (errno != ERANGE)
    {
      perror("next_mutpos failed");
      /* printf("errno = %d\n", errno); */
    }
    /* p = DBL_MAX; */
  }

#endif

  if (p < 0.0)
  {
    p = 0.0;
  }
  return ((LONG_MAX < p) ? LONG_MAX : (long) p);
}

