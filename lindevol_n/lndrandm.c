/* $Id: lndrandm.c,v 1.3 2000/01/30 03:25:18 kim Exp $ */
/*
 * $Log: lndrandm.c,v $
 * Revision 1.3  2000/01/30 03:25:18  kim
 * got rid of #include "lndrandm.h"
 *
 * Revision 1.2  2000/01/30 03:11:00  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

/* This is the random function that is used for all random values needed
   in the process of the actual simulation. It should not be used for any
   purposes outside the simulation, such as on-line analysis, generation
   of random samples etc. */

#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef __atarist__
#  include <unistd.h>
#endif

#ifndef TEST
#  include "lnderror.h"
#endif

#include "jklib.h"
#include "urandom.h"
#include "lndlib.h"

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif

static int  lndrandom_seed;
static unsigned long num_rndcalls = 0;
#ifdef DUMB_RANDOM
static long dmbrandom_state;
#else
static unsigned long rnd_state[64];
#endif


unsigned long lnd_random(unsigned long range)
{
#ifdef DUMB_RANDOM
  if (dmbrandom_state <= 0)
  { 
    dmbrandom_state <<= 1;
    dmbrandom_state ^= 0x1d872b41L;
  }
  else
    dmbrandom_state <<= 1;
  return (dmbrandom_state % range);
#else
  num_rndcalls++;
  return ((unsigned long) urandom_long(range));
#endif
}


double lnd_rnd(void)
{
#ifdef DUMB_RANDOM
  static const double randmax1 = 0x80000000U;

  if (dmbrandom_state <= 0)
  { 
    dmbrandom_state <<= 1;
    dmbrandom_state ^= 0x1d872b41L;
  }
  else
    dmbrandom_state <<= 1;
  return ((double) dmbrandom_state / randmax1 * 0.5 + 0.5);
#else
  num_rndcalls++;
  return (urandom_double());
#endif
}


#ifndef TEST
int write_rndgenerator_state(FILE *f)
{
  if (fprintf(f, "%d\n", lndrandom_seed) == EOF)
  {
    do_error("write_rndgenerator_state: error writing to file");
    return (-1);
  }
#ifdef DUMB_RANDOM
  if (fprintf(f, "%ld\n", dmbrandom_state) == EOF)
  {
    do_error("write_rndgenerator_state: error writing to file");
    return (-1);
  }
  else
  {
    /* printf("saved state %ld\n", dmbrandom_state); */
    return (0);
  }
#else
  fprintf(f, "%ld\n", num_rndcalls);
  return (write_urandom_state(f));
#endif
}


int read_rndgenerator_state(FILE *f)
{
  char buf[514];

  if (fgets(buf, 514, f) != buf)
  {
    do_error("read_rndgenerator_state: error reading state information");
    return (-1);
  }
  lndrandom_seed = strtol(buf, (char **) NULL, 10);
#ifdef DUMB_RANDOM
  if (fgets(buf, 514, f) != buf)
  {
    do_error("read_rndgenerator_state: error reading state information");
    return (-1);
  }
  dmbrandom_state = strtol(buf, (char **) NULL, 10);
  /* printf("restored state %ld\n", dmbrandom_state); */
  return (0);
#else
  if (fgets(buf, 514, f) != buf)
  {
    do_error("read_rndgenerator_state: error reading state information");
    return (-1);
  }
  num_rndcalls = strtol(buf, (char **) NULL, 10);
  return (read_urandom_state(f));
#endif
}
#endif /* TEST */


void seed_lnd_random(int seed)
{
  lndrandom_seed = seed;
#ifdef DUMB_RANDOM
  dmbrandom_state = seed;
#else
  ulong_initstate(lndrandom_seed, rnd_state);
#endif
}


/*
 * Randomly shuffles the num elements of array a
 */

void random_shuffle(long num, long *a)
{
  urandom_shuffle(num, sizeof(long), a);
}

