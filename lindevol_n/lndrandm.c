/* This is the random function that is used for all random values needed
   in the process of the actual simulation. It should not be used for any
   purposes outside the simulation, such as on-line analysis, generation
   of random samples etc. */

#include <errno.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef __atarist__
#  include <unistd.h>
#endif

#include "lnderror.h"

#include "jklib.h"
#include "lndlib.h"

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif

static int lndrandom_seed;
static unsigned long num_rndcalls = 0;
static unsigned long rnd_state[64];



/*** taken from the random() sources by Jan T. Kim ***/

/*
 * Copyright (c) 1983 Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that the above copyright notice and this paragraph are
 * duplicated in all such forms and that any documentation,
 * advertising materials, and other materials related to such
 * distribution and use acknowledge that the software was developed
 * by the University of California, Berkeley.  The name of the
 * University may not be used to endorse or promote products derived
 * from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */

/*
 * This is derived from the Berkeley source:
 *	@(#)random.c	5.5 (Berkeley) 7/6/88
 * It was reworked for the GNU C Library by Roland McGrath.
 */


/* An improved random number generation package.  In addition to the standard
   rand()/srand() like interface, this package also has a special state info
   interface.  The initstate() routine is called with a seed, an array of
   bytes, and a count of how many bytes are being passed in; this array is
   then initialized to contain information for random number generation with
   that much state information.  Good sizes for the amount of state
   information are 32, 64, 128, and 256 bytes.  The state can be switched by
   calling the setstate() function with the same array as was initiallized
   with initstate().  By default, the package runs with 128 bytes of state
   information and generates far better random numbers than a linear
   congruential generator.  If the amount of state information is less than
   32 bytes, a simple linear congruential R.N.G. is used.  Internally, the
   state information is treated as an array of longs; the zeroeth element of
   the array is the type of R.N.G. being used (small integer); the remainder
   of the array is the state information for the R.N.G.  Thus, 32 bytes of
   state information will give 7 longs worth of state information, which will
   allow a degree seven polynomial.  (Note: The zeroeth word of state
   information also has some other information stored in it; see setstate
   for details).  The random number generation technique is a linear feedback
   shift register approach, employing trinomials (since there are fewer terms
   to sum up that way).  In this approach, the least significant bit of all
   the numbers in the state table will act as a linear feedback shift register,
   and will have period 2^deg - 1 (where deg is the degree of the polynomial
   being used, assuming that the polynomial is irreducible and primitive).
   The higher order bits will have longer periods, since their values are
   also influenced by pseudo-random carries out of the lower bits.  The
   total period of the generator is approximately deg*(2**deg - 1); thus
   doubling the amount of state information has a vast influence on the
   period of the generator.  Note: The deg*(2**deg - 1) is an approximation
   only good for large deg, when the period of the shift register is the
   dominant factor.  With deg equal to seven, the period is actually much
   longer than the 7*(2**7 - 1) predicted by this formula.  */



/* For each of the currently supported random number generators, we have a
   break value on the amount of state information (you need at least thi
   bytes of state info to support this random number generator), a degree for
   the polynomial (actually a trinomial) that the R.N.G. is based on, and
   separation between the two lower order coefficients of the trinomial.  */

/* x**63 + x + 1.  */
#define	TYPE_4		4
#define	BREAK_4		256
#define	DEG_4		63
#define	SEP_4		1
#define MAX_TYPES	5

long int xxx_random(void);


/* Array versions of the above information to make code run faster.
   Relies on fact that TYPE_i == i.  */


static unsigned long int randtbl[DEG_4 + 1] = {
	0x00000004, 0x07a6c921, 0x75247e98, 0x958685e2,
	0xeeb14548, 0xd436aed1, 0xeeadcc28, 0xd5a6d6b0, 
	0x882e20e0, 0xc1165f36, 0x2fed1390, 0x44ca0be4,
	0x86da3f50, 0x8d76cf46, 0xcd46b0b8, 0x2d717398, 
	0x422d0ad8, 0xde7747c7, 0x5097ce88, 0xded84046,
	0xdde7cd78, 0xe47c24af, 0x194a39e8, 0xf4d833a8, 
	0x5ffe8570, 0x02787a04, 0x070b6cc0, 0x09f6a498,
	0x9513abc0, 0x2ea921bc, 0xe21aa4f8, 0x30df9bf0, 
	0x7dbc80a8, 0x33786e9d, 0x121d5b78, 0x8893446a,
	0xad28f228, 0x10db401d, 0x80b8ef28, 0x7c226980, 
	0x4f7bb080, 0xcb0bb342, 0x2a73d7f0, 0x2db9414c,
	0x1e9c84b0, 0x94232c82, 0x55f835b8, 0x44543d68, 
	0x7dd87cf8, 0xd78cf6a3, 0xffdd0f68, 0xb5011ece,
	0xcc390358, 0xe52c329b, 0xb47b15e8, 0xb2ef08b8, 
	0x52dc7210, 0xeeab5470, 0xf7a53f20, 0x7d7ace80,
	0x1e483a20, 0x53d0e918, 0x1dc50cf8, 0xaf493880, 
};

/* Fvoid * and Rvoid * are two pointers into the state info, a front and a rear
   pointer.  These two pointers are always rand_sep places aparts, as they
   cycle through the state information.  (Yes, this does mean we could get
   away with just one pointer, but the code for random is more efficient
   this way).  The pointers are left positioned as they would be from the call:
	initstate(1, randtbl, 128);
   (The position of the rear pointer, rptr, is really 0 (as explained above
   in the initialization of randtbl) because the state table pointer is set
   to point to randtbl[1] (as explained below).)  */

static unsigned long int *fptr = &randtbl[SEP_4 + 1];
static unsigned long int *rptr = &randtbl[1];



/* The following things are the pointer to the state information table,
   the type of the current generator, the degree of the current polynomial
   being used, and the separation between the two pointers.
   Note that for efficiency of random, we remember the first location of
   the state information, not the zeroeth.  Hence it is valid to access
   state[-1], which is used to store the type of the R.N.G.
   Also, we remember the last location, since this is more efficient than
   indexing every time to find the address of the last element to see if
   the front and rear pointers have wrapped.  */

static unsigned long int *state = randtbl + 1;

const static int rand_type = TYPE_4;
const static int rand_deg = DEG_4;
const static int rand_sep = SEP_4;

static unsigned long int *end_ptr = &randtbl[sizeof(randtbl) / sizeof(randtbl[0])];

/* Initialize the random number generator based on the given seed.  If the
   type is the trivial no-state-information type, just remember the seed.
   Otherwise, initializes state[] based on the given "seed" via a linear
   congruential generator.  Then, the pointers are set to known locations
   that are exactly rand_sep places apart.  Lastly, it cycles the state
   information a given number of times to get rid of any initial dependencies
   introduced by the L.C.R.N.G.  Note that the initialization of randtbl[]
   for default usage relies on values produced by this routine.  */

void xxx_srandom( unsigned int x)
{
  register long int i;
  state[0] = x;
  /* fprintf(stderr, "srandom: state[%3d] = %12ld (%10lx)\n", 0, state[0], state[0]); */
  for (i = 1; i < rand_deg; ++i)
  {
#ifdef LINUX_RANDOM
    state[i] = 1103515145L * state[i - 1] + 12345;
#else
    state[i] = 1103515245L * state[i - 1] + 12345;
#endif
    /* fprintf(stderr, "srandom: state[%3ld] = %12lu (%10lx)\n", i, state[i], state[i]); */
  }
  fptr = &state[rand_sep];
  rptr = &state[0];
  for (i = 0; i < 10 * rand_deg; ++i)
    (void) xxx_random();
}



/* Initialize the state information in the given array of N bytes for
   future random number generation.  Based on the number of bytes we
   are given, and the break values for the different R.N.G.'s, we choose
   the best (largest) one we can and set things up for it.  srandom is
   then called to initialize the state information.  Note that on return
   from srandom, we set state[-1] to be the type multiplexed with the current
   value of the rear pointer; this is so successive calls to initstate won't
   lose this information and will be able to restart with setstate.
   Note: The first thing we do is save the current state, if any, just like
   setstate so that it doesn't matter when initstate is called.
   Returns a pointer to the old state.  */

unsigned long *xxx_initstate(unsigned int seed, unsigned long *arg_state)
{
  unsigned long *ostate = &state[-1];

  state[-1] = (MAX_TYPES * (rptr - state)) + rand_type;
  state = arg_state + 1;	/* First location.  */
  /* Must set END_void * before srandom.  */
  end_ptr = state + rand_deg;
  xxx_srandom(seed);
  state[-1] = (MAX_TYPES * (rptr - state)) + rand_type;
  return ostate;
}



/* Restore the state from the given state array.
   Note: It is important that we also remember the locations of the pointers
   in the current state information, and restore the locations of the pointers
   from the old state information.  This is done by multiplexing the pointer
   location into the zeroeth word of the state information. Note that due
   to the order in which things are done, it is OK to call setstate with the
   same state as the current state
   Returns a pointer to the old state information.  */

unsigned long *xxx_setstate(unsigned long *arg_state)
{
  register unsigned long int *new_state = arg_state;
  register int type = new_state[0] % MAX_TYPES;
  register int rear = new_state[0] / MAX_TYPES;
  void *ostate = (void *) &state[-1];

  state[-1] = (MAX_TYPES * (rptr - state)) + rand_type;
  switch (type)
  {
  case TYPE_4:
    break;
  default:
    /* State info munged.  */
    fprintf(stderr, "xxx_setstate(): state info munged -- not changed\n");
    errno = EINVAL;
    return NULL;
  }
  state = new_state + 1;
  rptr = state + rear;
  fptr = state + (rear + rand_sep) % rand_deg;
  /* Set end_ptr too.  */
  end_ptr = state + rand_deg;
  return ostate;
}



/* If we are using the trivial TYPE_0 R.N.G., just do the old linear
   congruential bit.  Otherwise, we do our fancy trinomial stuff, which is the
   same in all ther other cases due to all the global variables that have been
   set up.  The basic operation is to add the number at the rear pointer into
   the one at the front pointer.  Then both pointers are advanced to the next
   location cyclically in the table.  The value returned is the sum generated,
   reduced to 31 bits by throwing away the "least random" low bit.
   Note: The code takes advantage of the fact that both the front and
   rear pointers can't wrap on the same call by not testing the rear
   pointer if the front one has wrapped.  Returns a 31-bit random number.  */

long int xxx_random(void)
{
  long int i;

  *fptr += *rptr;
  /* Chucking least random bit.  */
  i = (*fptr >> 1) & 0x7fffffffUL;
  ++fptr;
  if (fptr >= end_ptr)
  {
    fptr = state;
    ++rptr;
  }
  else
  {
    ++rptr;
    if (rptr >= end_ptr)
      rptr = state;
  }
  return i;
}

/*** end of code taken from random() sources ***/


unsigned long lnd_random(unsigned long range)
{
  num_rndcalls++;
  return (((unsigned long) xxx_random()) % range);
}


double lnd_rnd(void)
{
  num_rndcalls++;
  return ((double) (xxx_random() & 0x7fffff) / 0x800000);
}


int write_rndgenerator_state(FILE *f)
{
  int i;
  unsigned long junk[64];

  fprintf(stderr, "write_rndgenerator_state: starting\n");
  if (fprintf(f, "%d\n", lndrandom_seed) == EOF)
  {
    do_error("write_rndgenerator_state: error writing to file");
    return (-1);
  }
  fprintf(f, "%ld\n", num_rndcalls);
  /*** Note: it is necessary to call xxx_initstate() before saving rnd_state ***/
  /*** This call modifies rnd_state, without this modification resuming      ***/
  /*** does not work!                                                        ***/
  xxx_initstate(4711, junk);
  for (i = 0; i < 64; i++)
    fprintf(f, "%lu\n", rnd_state[i]);
  xxx_setstate(rnd_state);
  return (0);
}


int read_rndgenerator_state(FILE *f)
{
  int i;
  char buf[514];
  unsigned long junk[64];

  fprintf(stderr, "read_rndgenerator_state: starting\n");
  xxx_initstate(4711, junk);
  if (fgets(buf, 514, f) != buf)
  {
    do_error("read_rndgenerator_state: error reading state information (seed)");
    return (-1);
  }
  lndrandom_seed = strtol(buf, (char **) NULL, 10);
  if (fgets(buf, 514, f) != buf)
  {
    do_error("read_rndgenerator_state: error reading state information (num_rndcalls)");
    return (-1);
  }
  num_rndcalls = strtol(buf, (char **) NULL, 10);
  for (i = 0; i < 64; i++)
  {
    if (fgets(buf, 514, f) != buf)
    {
      do_error("read_rndgenerator_state: error reading state information (state array)");
      return (-1);
    }
    rnd_state[i] = strtoul(buf, NULL, 10);
  }
  xxx_setstate(rnd_state);
/*
  xxx_srandom(lndrandom_seed);
  for (i = 0; i < num_rndcalls; i++)
  {
    (void) xxx_random();
  }
*/
  /* printf("restored to state string:\n%s\n", buf); */
  /* printf("next random value: %d\n", xxx_random()); */
  return (0);
}


void seed_lnd_random(int seed)
{
  lndrandom_seed = seed;
  xxx_initstate(lndrandom_seed, rnd_state);
}


/*
 * Randomly shuffles the num elements of array a
 */

void random_shuffle(long num, long *a)
{
  long i, r, a_r;

  for (i = 0; i < num; i++)
  {
    r = lnd_random(num);
    a_r = a[r];
    a[r] = a[i];
    a[i] = a_r;
  }
}


#ifdef TEST

int main(int argc, char *argv[])
{
  char rstate[256];
  long i, lr, r, l1[10], l2[10];
  unsigned long virginstate[64], junk[64];
  FILE *f;

  
  xxx_initstate(1, virginstate);
  xxx_initstate(1, junk);
  f = fopen("vstate.dat", "w");
  for (i = 0; i < 64; i += 8)
  {
    for (r = i; r < i + 8; r++)
      fprintf(f, "0x%08lx, ", virginstate[r]);
    fprintf(f, "\n");
  }
  fclose(f);
  initstate(1, rstate, 256);
  seed_lnd_random(1);
  for (i = 0; i < 10; i++)
  {
    lr = xxx_random();
    r = random();
    printf("%3ld: %12ld %12ld\n", i, lr, r);
  }
  printf("\n");
  seed_lnd_random(12345);
  for (i = 0; i < 10; i++)
  {
    lr = lnd_random(100);
    printf("%3ld: %10ld\n", i, lr);
  }
  printf("\n");
  f = fopen("l.rnd", "w");
  write_rndgenerator_state(f);
  fclose(f);
  for (i = 0; i < 10; i++)
  {
    l1[i] = lnd_random(100);
    printf("%3ld: %10ld\n", i, l1[i]);
  }
  printf("\n");
  for (i = 0; i < 10; i++)
  {
    lr = lnd_random(100);
    printf("%3ld: %10ld\n", i, lr);
  }
  printf("\n");
  f = fopen("l.rnd", "r");
  read_rndgenerator_state(f);
  fclose(f);
  f = fopen("lx.rnd", "w");
  write_rndgenerator_state(f);
  fclose(f);
  for (i = 0; i < 10; i++)
  {
    l2[i] = lnd_random(100);
    printf("%3ld: %10ld\n", i, l2[i]);
  }
  printf("\n");
  for (i = 0; i < 10; i++)
  {
    lr = lnd_random(100);
    printf("%3ld: %10ld\n", i, lr);
  }
  printf("\n");
  return (EXIT_SUCCESS);
}


#endif

