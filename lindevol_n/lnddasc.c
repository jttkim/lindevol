/* $Id: lnddasc.c,v 1.2 2000/01/30 03:10:59 kim Exp $ */
/*
 * $Log: lnddasc.c,v $
 * Revision 1.2  2000/01/30 03:10:59  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

/* This module contains the functions for displaying the LindEvol world.
   At present, there is just a primitive ASCII displaying system. */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif

#include "lnderror.h"
#include "lndvals.h"
#include "lndtypes.h"

#include "lndglobl.h"
#include "lnddispl.h"


static long curr_generation;

FILE  *world_file = NULL;
char   world_file_name[MAX_SLEN];


int open_world_file(const char *simname, const char *mode)
{
  sprintf(world_file_name, "%s%swld", simname, ".");
  world_file = fopen(world_file_name, mode);
  if (world_file == NULL)
  {
    perror("open_world_file failed");
    return (-1);
  }
  return (0);
}


void close_world_file(void)
{
  if (world_file != NULL)
    fclose(world_file);
}


void show_world(long generation, FILE *f)
{
  long x, y, l, p;

  fprintf(f, "world at generation %ld, sized %ld * %ld\n", curr_generation, world_width, world_height);
  for (y = world_height - 1;y >= 0 ; y--)
  {
    for (x = 0; x < world_width; x++)
    {
      if ((p = world[x][y].plant_no) > -1)
      {
        if (plant[p]->cell[world[x][y].cell_no].energy)
        {
          fprintf(f, "%c", (char) ('A' + p % 26));
        }
        else
        {
          fprintf(f, "%c", (char) ('a' + p % 26));
        }
      }
      else
      {
        fprintf(f, " ");
      }
    }
    fprintf(f, "\n");
  }
  fprintf(f, "\n");
  for (l = 0; l < 6; l++)
  {
    for (x = 0; x < l; x++) fprintf(f, " ");
    for (x = l; x < world_width; x+= 6)
    {
      if ((p = world[x][0].plant_no) > -1)
      {
        if (world[x][0].cell_no == 0)
        {
          fprintf(f, "%-5ld ", plant[p]->cellular_energy);
        }
        else
        {
          fprintf(f, "      ");
        }
      }
      else
      {
        fprintf(f, "      ");
      }
    }
    fprintf(f, "\n");
  }
  /* fprintf(f, "\nmax. fitness: %ld, min. fitness: %ld\n", genome[gi[0]].fitness, genome[gi[psize - 1]].fitness); */
}


void display_start(long generation)
{
  curr_generation = generation;
}


void display_done(void)
{
  if (!quietmode)
    show_world(curr_generation, stdout);
  if ((worldmode != 0) && (world_file != NULL))
    show_world(curr_generation, world_file);
}


void display_plant_killed(long plant_no)
{
}


void display_plant_mutated(long plant_no)
{
}


void display_cell(long x, long y)
{
}


void poll_user_interface(void)
{
}


void init_world_display(const char *simname)
{
  if (worldmode)
  {
    open_world_file(simname, "a");
  }
}


void close_world_display(void)
{
  if (worldmode)
  {
    close_world_file();
  }
}

