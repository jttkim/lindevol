#include <ctype.h>
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

#include "jklib.h"
#include "genomelib.h"
#include "lndtypes.h"
#include "lnddispl.h"
#include "lnderror.h"
#include "lndsignl.h"
#include "lndglobl.h"

#include "gntypes.h"
#include "gnlib.h"

#include "lndlibin.h"

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif


long num_energyrich_cells(long plant_no)
{
  long cell_no;
  long f = 0;

  for (cell_no = 0; cell_no < plant[plant_no]->num_cells; cell_no++)
  {
    if (plant[plant_no]->cell[cell_no].energy)
    {
      f++;
    }
  }
  /* printf("num_energyrich_cells: plant %lu has %lu cells, f=%lu\n", plant_no, plant[plant_no].num_cells, f); */
  return (f);
}


void add_nutrient_random(void)
{
  long x, y;

  for ( ; ; nutrient_i++)
  {
    if (nutrient_i == world_width * (world_soil + 1))
      nutrient_i = 0;
    if (nutrient_i == 0)
      random_shuffle(world_width * (world_soil + 1), soil_index);
    x = soil_index[nutrient_i] % world_width;
    y = soil_index[nutrient_i] / world_width;
    if (!world[x][y].nutrient)
    {
      world[x][y].nutrient = 1;
      break;
    }
  }
}


void add_organic_nutrient(long x, long y, long n)
{
  if (y > world_soil)
    y = world_soil;
  world[x][y].organic_nutrient += n;
}


/*
 * Returns a randomly chosen free position on the bottom of the world.
 * If none exists, -1 is returned.
 */

long random_free_position(void)
{
  long pi = 0;

  random_shuffle(world_width, r_index);
  while (world[r_index[pi]][world_soil].plant_no > -1)
  {
    pi++;
    if (pi == world_width)
    {
      return (-1);
    }
  }
  return (r_index[pi]);
}


/*
 * Returns the X-coordinate of a free position on the bottom of the
 * world that is as close as possible to pos. If free positions next
 * to pos exist equally distant to the left or to the right, it is
 * chosen randomly which one is returned. If no free site exists at
 * the bottom of the world, -1 is returned.
 */

long next_free_position(long pos)
{
  long offset;

  if (world[pos][world_soil].plant_no == -1)
  {
    return (pos);
  }
  if (lnd_random(2))
  {
    for (offset = 1; offset <= world_width / 2; offset++)
    {
      if (world[(pos + offset) % world_width][world_soil].plant_no == -1)
        return ((pos + offset) % world_width);
      if (world[(world_width + pos - offset) % world_width][world_soil].plant_no == -1)
        return ((world_width + pos - offset) % world_width);
    }
  }
  else
  {
    for (offset = 1; offset <= world_width / 2; offset++)
    {
      if (world[(world_width + pos - offset) % world_width][world_soil].plant_no == -1)
        return ((world_width + pos - offset) % world_width);
      if (world[(pos + offset) % world_width][world_soil].plant_no == -1)
        return ((pos + offset) % world_width);
    }
  }
  return (-1);
}


/*
 * create a plant descending from genome *g. The created plant consists
 * of one germ cell at (pos, 0). generation specifies the generation in
 * which the plant is created, this is stored in the node information.
 * return codes: -1 indicates failure, 0 or positive numbers are indices
 * in the plant array.
 */

long create_plant(long pos, long generation, const GENOME *parent)
{
  long plant_no;

  for (plant_no = 0; (plant_no < world_width) && plant[plant_no]; plant_no++)
    ;
  if (plant_no == world_width)
  {
    do_error("create_plant: trying to create too many plants");
    return (-1);
  }
  plant[plant_no] = (PLANT *) malloc(sizeof(PLANT));
  if (plant[plant_no] == NULL)
  {
    do_error("create_plant: not enough memory for plant structure");
    return (-1);
  }
#ifdef LENGTHCHECK
  if ((parent->length > LENGTHCHECK) || (parent->num_genes > LENGTHCHECK / 2))
    fprintf(stderr, "generation #%5ld: duplication of l = %6ld, n = %6ld\n", generation, parent->length, parent->num_genes);
#endif
  if (duplicate_genome(&(plant[plant_no]->genome), parent, 0xffffffff) < 0)
  {
    do_error("create_plant: failed to duplicate genome (out of memory?)");
    free0(plant[plant_no]);
    return (-1);
  }
  plant[plant_no]->genome.mut_flag = 0;
  plant[plant_no]->genome.num_mutations = 0;
  plant[plant_no]->cellular_energy = 0;
  plant[plant_no]->cellular_nutrient = 0;
  plant[plant_no]->energy_pool = 0;
  plant[plant_no]->nutrient_pool = 0;
  plant[plant_no]->num_cells = 0;
  plant[plant_no]->cell = NULL;
  plant[plant_no]->age = 0;
  create_cell(plant_no, pos, world_soil);
  psize++;
  return (plant_no);
}


/*
 * Remove plant specified by plant_no.
 */

void remove_plant(long plant_no)
{
  long cell_no, x, y;

  if (phyltest_savefreq > 0)
  {
    if (gn_node_death(&gntree, &(plant[plant_no]->node), generation) < 0)
      do_error("remove_plant: error updating genealogy node");
  }
  x = plant[plant_no]->cell[0].x;
  y = plant[plant_no]->cell[0].y;
  add_organic_nutrient(x, y, plant[plant_no]->nutrient_pool);
  for (cell_no = 0; cell_no < plant[plant_no]->num_cells; cell_no++)
  {
    x = plant[plant_no]->cell[cell_no].x;
    y = plant[plant_no]->cell[cell_no].y;
    if (plant[plant_no]->cell[cell_no].nutrient)
      add_organic_nutrient(x, y, 2);
    else
      add_organic_nutrient(x, y, 1);
    world[x][y].plant_no = -1;
    display_cell(x, y);
  }
  display_plant_killed(plant_no);
  free_genome(&(plant[plant_no]->genome));
  free(plant[plant_no]->cell);
  free0(plant[plant_no]);
  psize--;
}


/*
 * Note: seed production costs one nutrient unit. If no new
 *     plant is created, the nutrient is released into the
 *     soil.
 */

int create_seed(long plant_no, long cell_no, long pos)
{
  long newplant_no;

  if ((newplant_no = create_plant(pos, generation, &(plant[plant_no]->genome))) > -1)
  {
    num_new_plants++;
    if (phyltest_savefreq > 0)
    {
      if (gn_new_treenode(&gntree, &(plant[plant_no]->node), generation, NULL, NULL, &(plant[newplant_no]->node)) < 0)
      {
        do_error("create_seed: failed to create new genealogy node");
        return (-1);
      }
    }
  }
  return (0);
}


int create_flying_seed(long plant_no, long cell_no)
{
  long pos;

  num_seeds++;
  num_flying_seeds++;
  pos = random_free_position();
  if (pos < 0)
  {
    add_organic_nutrient(lnd_random(world_width), world_soil, 1);
    return (0);
  }
  return (create_seed(plant_no, cell_no, pos));
}


int create_local_seed(long plant_no, long cell_no)
{
  long pos;

  num_seeds++;
  num_local_seeds++;
  pos = next_free_position(plant[plant_no]->cell[cell_no].x);
  if (pos < 0)
  {
    add_organic_nutrient(plant[plant_no]->cell[cell_no].x, world_soil, 1);
    return (0);
  }
  return (create_seed(plant_no, cell_no, pos));
}


void sunshine(void)
{
  long x, y;

  for (x = 0; x < world_width; x++)
  {
    for (y = world_height - 1; y >= world_soil; y--)
    {
      if ((world[x][y].plant_no > -1) && (lnd_random(2) == 1))
      {
        if (!(plant[world[x][y].plant_no]->cell[world[x][y].cell_no].energy))
        {
          plant[world[x][y].plant_no]->cell[world[x][y].cell_no].energy = 1;
          plant[world[x][y].plant_no]->cellular_energy++;
          display_cell(x, y);
        }
        break;
      }
    }
  }
}


void nutrient_diffusion(void)
{
  long x, y, xnew, ynew, r, i;

  for (i = 0; i < world_width * (world_soil + 1); i++)
    tmp_index[i] = i;
  random_shuffle(world_width * (world_soil + 1), tmp_index);
  for (i = 0; i < world_width * (world_soil + 1); i++)
  {
    x = tmp_index[i] % world_width;
    y = tmp_index[i] / world_width;
    if (!world[x][y].nutrient && (world[x][y].organic_nutrient > 0) && (lnd_rnd() < decomposition_rate))
    {
      world[x][y].organic_nutrient--;
      world[x][y].nutrient = 1;
    }
    if (world[x][y].nutrient && (lnd_rnd() < diffusion_rate))
    {
      r = lnd_random(8);
      xnew = (world_width + x + x_offset[r]) % world_width;
      ynew = y + y_offset[r];
      if ((ynew >= 0) && (ynew <= world_soil) && !world[xnew][ynew].nutrient)
      {
        world[x][y].nutrient = 0;
        world[xnew][ynew].nutrient = 1;
      }
    }
  }
  for (x = 0; x < world_width; x++)
  {
    for (y = 0; y <= world_soil; y++)
    {
      if ((world[x][y].nutrient) && (world[x][y].plant_no > -1))
      {
        if (!plant[world[x][y].plant_no]->cell[world[x][y].cell_no].nutrient)
        {
          plant[world[x][y].plant_no]->cell[world[x][y].cell_no].nutrient = 1;
          plant[world[x][y].plant_no]->cellular_nutrient++;
          world[x][y].nutrient = 0;
        }
      }
    }
  }
}

