#include <stdio.h>

long random_seed;
long random_number = 0;

long random(void)
{
  random_number++;
  if (random_seed <= 0)
  { 
    random_seed <<= 1;
    random_seed ^= 0x1d872b41L;
  }
  else
    random_seed <<= 1;
  return random_seed;
}


unsigned long lnd_random(unsigned long range)
{
  return ((random() & 0x7fffffff) % range);
}


double lnd_rnd(void)
{
  long r;

  r = random() & 0x7fffff;
  return (((double) r) / ((double) 0x800000));
}


