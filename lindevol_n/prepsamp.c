#include "lndglobl.h"
#include "lnderror.h"
#include "lndlib.h"
#include "lndlibin.h"


/*
 * Prepare an array of indices of sample_size
 * randomly chosen plants. If population size
 * is below sample_size, the entire population
 * in random order is chosen.
 * Return value: size of sample.
 * Note: prepare_sample uses random_shuffle, which
 *     uses the simulator's random number generator.
 *     To avoid changing the random number genrator's
 *     state, the original state is saved, copied,
 *     the copy is used for shuffling, and the
 *     original state is finally restored.
 */

long prepare_sample(long sample_size, long *sample_index)
{
  long s = 0, i;
  unsigned long junk[64], rstate[64], *saved_state;

  for(i = 0; i < world_width; i++)
    tmp_index[i] = i;
  saved_state = xxx_initstate(4711, junk);
  for (i = 0; i < 64; i++)
    rstate[i] = saved_state[i];
  xxx_setstate(rstate);
  random_shuffle(world_width, tmp_index);
  xxx_setstate(saved_state);
  for(i = 0; i < world_width; i++)
  {
    if (plant[tmp_index[i]])
      sample_index[s++] = tmp_index[i];
  }
  if (s < world_width)
    sample_index[s] = -1;
  return ((s < sample_size) ? s : sample_size);
}

