#include <stdio.h>
#include <stdlib.h>

#include "genomelib.h"

#include "lndvals.h"
#include "lndtypes.h"
#include "lndglobl.h"
#include "lndlibin.h"
#include "lndlib.h"

int create_arrays(void)
{
  long i, j;

  plant = (PLANT **) malloc(world_width * sizeof(PLANT *));
  if (plant == NULL)
  {
    return (-1);
  }
  for (i = 0; i < world_width; i++)
  {
    plant[i] = NULL;  /* init pointers to NULL, otherwise create_plant gets confused with the junk */
  }
  species_d = (SPECIES_D *) malloc(world_width * sizeof(SPECIES_D));
  if (species_d == NULL)
  {
    return (-1);
  }
  world = (LATTICE_SITE **) malloc(world_width * sizeof (LATTICE_SITE *));
  if (world == NULL)
  {
    return (-1);
  }
  world[0] = (LATTICE_SITE *) malloc(world_width * world_height * sizeof (LATTICE_SITE));
  if (world[0] == NULL)
  {
    return (-1);
  }
  for (i = 0; i < world_width; i++)
  {
    world[i] = world[0] + i * world_height;
    for (j = 0; j < world_height; j++)
    {
      world[i][j].plant_no = -1; /* init array to plant no = -1, indicating free positions */
      world[i][j].organic_nutrient = 0;
      world[i][j].nutrient = 0;
    }
  }
  r_index = (long *) malloc(world_width * sizeof(long));
  if (r_index == NULL)
  {
    return (-1);
  }
  pl_index = (long *) malloc(world_width * sizeof(long));
  if (pl_index == NULL)
  {
    return (-1);
  }
  for (i = 0; i < world_width; i++)
  {
    r_index[i] = i;
    pl_index[i] = i;
  }
  tmp_index = (long *) malloc(world_width * world_height * sizeof(long));
  if (tmp_index == NULL)
  {
    return (-1);
  }
  sample_index = (long *) malloc(world_width * sizeof(long));
  if (sample_index == NULL)
  {
    return (-1);
  }
  soil_index = (long *) malloc(world_width * (world_soil + 1) * sizeof(long));
  if (soil_index == NULL)
  {
    return (-1);
  }
  for (i = 0; i < world_width * (world_soil + 1); i++)
    soil_index[i] = i;
  soil_profile_h = (long *) malloc(world_width * sizeof(long));
  if (soil_profile_h == NULL)
    return (-1);
  soil_profile_v = (long *) malloc((world_soil + 1) * sizeof(long));
  if (soil_profile_v == NULL)
    return (-1);
  organic_profile_h = (long *) malloc(world_width * sizeof(long));
  if (organic_profile_h == NULL)
    return (-1);
  organic_profile_v = (long *) malloc((world_soil + 1) * sizeof(long));
  if (organic_profile_v == NULL)
    return (-1);
  return (0);
}


void clear_arrays(void)
{
  long i;

  for (i = 0; i < world_width; i++)
  {
    if (plant[i] != NULL)
    {
      free_genome(&(plant[i]->genome));
      free(plant[i]->cell);
      free0(plant[i]);
    }
  }
  free0(plant);
  free0(species_d);
  free0(world[0]);
  free0(world);
  free0(r_index);
  free0(pl_index);
  free0(tmp_index);
  free0(sample_index);
  free0(soil_index);
  free0(soil_profile_h);
  free0(soil_profile_v);
  free0(organic_profile_h);
  free0(organic_profile_v);
}

