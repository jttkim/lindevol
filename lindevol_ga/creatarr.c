#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lndtypes.h"
#include "lnderror.h"

#include "genomelib.h"
#include "gnlib.h"

#include "lndglobl.h"


int create_arrays()
{
  unsigned long i;

  lnd_genome = (LND_GENOME *) malloc(psize * sizeof(LND_GENOME));
  if (lnd_genome == NULL)
    return (-1);
  gi = (unsigned long *) malloc(psize * sizeof(unsigned long));
  if (gi == NULL)
    return (-1);
  gi_tmp = (unsigned long *) malloc(psize * sizeof(unsigned long));
  if (gi_tmp == NULL)
    return (-1);
  species_d = (SPECIES_D *) malloc(psize * sizeof(SPECIES_D));
  if (species_d == NULL)
    return (-1);
  world = (LATTICE_SITE **) malloc(world_width * sizeof (LATTICE_SITE *));
  if (world == NULL)
    return (-1);
  world[0] = (LATTICE_SITE *) malloc(world_width * world_height * sizeof (LATTICE_SITE));
  if (world[0] == NULL)
    return (-1);
  for (i = 1; i < world_width; i++)
    world[i] = world[0] + i * world_height;
  plant = (PLANT *) malloc(psize * sizeof(PLANT));
  if (plant == NULL)
    return (-1);
  pl_index = (unsigned long *) malloc(psize * sizeof(unsigned long));
  if (pl_index == NULL)
    return (-1);
  for (i = 0; i < psize; i++)
  {
    plant[i].num_cells = 0;
    gi[i] = i;
    gi_tmp[i] = i;
    pl_index[i] = i;
  }
  return (0);
}


/* free0() frees a block of memory pointed to by *ptr if *ptr != NULL.
   After freeing, *ptr is set to NULL */

void free0(void **ptr)
{
  if (*ptr != NULL)
  {
    free(*ptr);
    *ptr = NULL;
  }
}


void clear_arrays(void)
{
  long i;

  for (i = 0; i < psize; i++)
    free_genome(&(lnd_genome[i].genome));
  free0((void *) (&lnd_genome));
  free0((void *) (&gi));
  free0((void *) (&gi_tmp));
  free0((void *) (&species_d));
  free0((void *) (&world[0]));
  free0((void *) (&world));
  free0((void *) (&plant));
  free0((void *) (&pl_index));
}

