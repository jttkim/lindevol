/* $Id: pgrow6.c,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: pgrow6.c,v $
 * Revision 1.2  2000/01/30 03:11:00  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

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

#include "lndtypes.h"
#include "lnddispl.h"
#include "lnderror.h"
#include "lndsignl.h"
#include "lndlib.h"
#include "lnd6.h"
#include "lndglobl.h"

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif


long create_cell(long p, long x, long y)
{
  long c;
  PL_CELL *new_cell;

  if(world[x][y].plant_no > -1)
  {
    do_error("create_cell: ERROR: cell production on occupied site\n");
    fprintf(stderr, "    site (%ld, %ld) is already occupied by plant %ld\n", x, y, world[x][y].plant_no);
    return (-1);
  }
  if (plant[p]->cell)
  {
    /* new_cell = realloc(plant[p]->cell, (plant[p]->num_cells + 1) * sizeof(PL_CELL)); */
    new_cell = malloc((plant[p]->num_cells + 1) * sizeof(PL_CELL));
    if (new_cell)
    {
      memcpy(new_cell, plant[p]->cell, plant[p]->num_cells * sizeof(PL_CELL));
      free(plant[p]->cell);
    }
  }
  else
  {
    new_cell = malloc((plant[p]->num_cells + 1) * sizeof(PL_CELL));
  }
  if (new_cell == 0)
  {
    do_error("create_cell: failed because out of memory\n");
    return (-1);
  }
  plant[p]->cell = new_cell;
  c = plant[p]->num_cells;
  world[x][y].plant_no = p;
  world[x][y].cell_no = c;
  plant[p]->cell[c].x = x;
  plant[p]->cell[c].y = y;
  plant[p]->cell[c].energy = 0;
  plant[p]->cell[c].nutrient = 0;
  plant[p]->cell[c].state = 0;
  plant[p]->cell[c].next_state = 0;
  plant[p]->num_cells++;
  display_cell(x, y);
  return (c);
}


/*
 * divide() returns -1 if the plant has killed itself in the process of
 * trying to perform a cell division.
 * Return: < 0 : error (failed malloc() or position out of world)
 *         >= 0: success, value is array index of newly created cell
 * New in lnd6: new cell's state for next time step is copied from
 *     dividing cell's (internal) state for next time step.
 */

long divide(long plant_no, long cell_no, int direction)
{
  long c = -1;
  long x, y, otherplant_no;

  /* printf("divide: dividing plant %lu, cell %lu (%lu/%lu) to direction %d\n", plant_no, cell_no, plant[plant_no]->cell[cell_no].x, plant[plant_no]->cell[cell_no].y, direction); */
/*
  if (cell_no >= plant[plant_no]->num_cells)
  {
    printf("  error: dividing non-existent cell #%lu of %lu\n", cell_no, plant[plant_no]->num_cells);
  }
*/
  num_divisions++;
  x = (plant[plant_no]->cell[cell_no].x + x_offset[direction] + world_width) % world_width;
  y = plant[plant_no]->cell[cell_no].y + y_offset[direction];
  /* printf("new location: (%ld/%ld) = (%lu + %d/%lu + %d)\n", x, y, plant[plant_no]->cell[cell_no].x, x_offset[direction], plant[plant_no]->cell[cell_no].y, y_offset[direction]); */
  if ((y >= 0) && (y < world_height))
  {
    if ((otherplant_no = world[x][y].plant_no) > -1)
    {
      num_attacks++;
      if (lnd_random(plant[otherplant_no]->cellular_energy + 1) == 0)
      {
        remove_plant(otherplant_no);
        num_deaths++;
        if (otherplant_no == plant_no)
        {
          return (-1);
        }
      }
    }
    if (world[x][y].plant_no == -1)
    {
      /* printf("divide: dividing to (%ld/%ld)\n", x, y); */
      num_new_cells++;
      c = create_cell(plant_no, x, y);
      plant[plant_no]->cell[c].next_state = plant[plant_no]->cell[cell_no].next_state;
    }
    /* else printf("divide: site (%ld/%ld) is occupied by plant %ld\n", x, y, world[x][y].plant_no); */
  }
  return (c);
}


/*
 * interpret the action code determined for a cell.
 */

void cell_activity(long plant_no, long cell_no, long action, unsigned char gene_output)
{
  if (action == LND_DIVIDE)
  {
    if (plant[plant_no]->cell[cell_no].energy)
    {
      plant[plant_no]->cell[cell_no].energy = 0;
      plant[plant_no]->cellular_energy--;
      display_cell(plant[plant_no]->cell[cell_no].x, plant[plant_no]->cell[cell_no].y);
      if (plant[plant_no]->cell[cell_no].nutrient)
      {
	if (divide(plant_no, cell_no, gene_output & 0x07) > -1)
	{
	  plant[plant_no]->cell[cell_no].nutrient = 0;
	  plant[plant_no]->cellular_nutrient--;
	}
      }
    }
  }
  else if ((action == LND_FLYINGSEED) && (plant[plant_no]->num_cells >= seedprod_threshold))
  {
    if (plant[plant_no]->cell[cell_no].energy)
    {
      plant[plant_no]->cell[cell_no].energy = 0;
      plant[plant_no]->cellular_energy--;
      display_cell(plant[plant_no]->cell[cell_no].x, plant[plant_no]->cell[cell_no].y);
      if (plant[plant_no]->cell[cell_no].nutrient)
      {
	create_flying_seed(plant_no, cell_no);
	plant[plant_no]->cell[cell_no].nutrient = 0;
	plant[plant_no]->cellular_nutrient--;
      }
    }
  }
  else if ((action == LND_LOCALSEED) && (plant[plant_no]->num_cells >= seedprod_threshold))
  {
    if (plant[plant_no]->cell[cell_no].energy)
    {
      plant[plant_no]->cell[cell_no].energy = 0;
      plant[plant_no]->cellular_energy--;
      display_cell(plant[plant_no]->cell[cell_no].x, plant[plant_no]->cell[cell_no].y);
      if (plant[plant_no]->cell[cell_no].nutrient)
      {
	create_local_seed(plant_no, cell_no);
	plant[plant_no]->cell[cell_no].nutrient = 0;
	plant[plant_no]->cellular_nutrient--;
      }
    }
  }
  else if (action == LND_MUTMINUS)
  {
    if (plant[plant_no]->cell[cell_no].energy)
    {
      plant[plant_no]->cell[cell_no].energy = 0;
      plant[plant_no]->cellular_energy--;
      display_cell(plant[plant_no]->cell[cell_no].x, plant[plant_no]->cell[cell_no].y);
      num_mutminus++;
      plant[plant_no]->genome.mut_flag--;
    }
  }
  else if (action == LND_MUTPLUS)
  {
    if (plant[plant_no]->cell[cell_no].energy)
    {
      plant[plant_no]->cell[cell_no].energy = 0;
      plant[plant_no]->cellular_energy--;
      display_cell(plant[plant_no]->cell[cell_no].x, plant[plant_no]->cell[cell_no].y);
      num_mutplus++;
      plant[plant_no]->genome.mut_flag++;
    }
  }
  else
  {
    do_error("Trying to interpret illegal action code");
  }
}


static void translate_genome(const GENOME *genome, GENE_SPEC *gene_spec)
{
  long gpos, gene_no;

  gpos = 0;
  gene_no = 0;
  while (gpos < genome->length)
  {
    while ((gpos < genome->length) && !l4_promoter(genome->g[gpos]))
      gpos++;
    gpos++;
    if (gene_no < genome->num_genes)
      gene_spec[gene_no].output = -1;
    if (gpos >= genome->length)
      break;
    gpos = compute_statespec(genome, gpos, &(gene_spec[gene_no].statespec), num_bitchecks);
    if (gpos >= genome->length)
      break;
    if (l4_terminator(genome->g[gpos]))
      gene_spec[gene_no].output = genome->g[gpos];
    gene_no++;
  }
}


void plant_growth(long plant_no)
{
  long num_existing_cells, cell_no, action;
  unsigned char gene_output;
  double death_p;
  GENE_SPEC *gene_spec = NULL;

  if (plant[plant_no]->age > 0)
  {
    num_grownplants++;
    death_p = death_probability(plant_no, p_random_death, rdeath_f_energy, rdeath_f_numcells, leanover_penalty);
    min_deathprob = (min_deathprob < death_p) ? min_deathprob : death_p;
    max_deathprob = (max_deathprob > death_p) ? max_deathprob : death_p;
    sum_deathprob += death_p;
    if (lnd_rnd() < death_p)
    {
      remove_plant(plant_no);
      num_deaths++;
    }
    else
    {
      plant[plant_no]->genome.mut_flag = 0;
#ifndef SUPPRESS_PRETRANSLATION
      if (plant[plant_no]->genome.num_genes)
      {
        if ((gene_spec = (GENE_SPEC *) malloc(plant[plant_no]->genome.num_genes * sizeof(GENE_SPEC))) != NULL)
          translate_genome(&(plant[plant_no]->genome), gene_spec);
        else
          do_error("plant_growth: no pretranslation of genome because malloc() failed");
      }
#endif
      num_existing_cells = plant[plant_no]->num_cells;
      for (cell_no = 0 ; cell_no < num_existing_cells; cell_no++)
      {
        action = process_cell(plant_no, cell_no, &gene_output, gene_spec);
	if (action > -1)
	{
	  cell_activity(plant_no, cell_no, action, gene_output);
	  if (plant[plant_no] == NULL)  /* if plant has killed itself */
	    break;
	}
      }
      if (gene_spec)
        free0(gene_spec);
    }
    if (plant[plant_no] != NULL)
    {
      for (cell_no = 0; cell_no < num_existing_cells; cell_no++)
      {
        plant[plant_no]->cell[cell_no].state = plant[plant_no]->cell[cell_no].next_state;
      }
      mutate_genome(&(plant[plant_no]->genome), m_replacement, m_insertion, m_deletion, m_duplication, m_factor);
      display_plant_mutated(plant_no);
    }
  }
  if (plant[plant_no] != NULL)
  {
    plant[plant_no]->age++;
  }
}

