/* This module contains the functions for displaying the LindEvol world.
   At present, there is just a primitive ASCII displaying system. */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "lnderror.h"
#include "lndvals.h"
#include "lndtypes.h"


extern char quietmode;
extern char worldmode;

extern unsigned long world_width, world_height; /* these are defined and */
extern unsigned long psize;

extern PLANT *plant;                            /* created in the main module */
extern LATTICE_SITE **world;
extern LND_GENOME *lnd_genome;
extern unsigned long *gi;


static unsigned long curr_generation;

static char **w_asc = (char **) NULL;

FILE  *world_file;
char   world_file_name[MAX_SLEN];


int open_world_file(const char *simname)
{
  sprintf(world_file_name, "%s%swld", simname, ".");
  world_file = fopen(world_file_name, "w");
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


void show_world(unsigned long generation, FILE *f)
{
  long x, y, i, s, p;

  fprintf(f, "world at generation %lu, sized %lu * %lu\n", generation, world_width, world_height);
  for (y = world_height - 1; ; y--)
  {
    for (x = 0; x < world_width; x++)
    {
      fprintf(f, "%c", (int) w_asc[x][y]);
    }
    fprintf(f, "\n");
    if (y == 0) break;
  }
  fprintf(f, "\n");
  s = 6L / ((long) world_width / (long) psize);
  s = (s < 1) ? 1 : s;
  for (i = 0; i < s; i++)
  {
    p = 0;
    for (x = i; x < psize; x += s)
    {
      for ( ;p < plant[x].cell[0].x; p++)
        fprintf(f, " ");
      fprintf(f, "%-5lu ", lnd_genome[gi[x]].fitness);
      p += 6;
    }
    fprintf(f, "\n");
  }
  /* fprintf(f, "\nmax. fitness: %lu, min. fitness: %lu\n", genome[gi[0]].fitness, genome[gi[psize - 1]].fitness); */
}


void display_on(unsigned long generation)
{
  unsigned long x, y;

  curr_generation = generation;
  for (x = 0; x < world_width; x++)
  {
    for (y = 0; y < world_height; y++)
    {
      w_asc[x][y] = ' ';
    }
  }
}


void display_off(void)
{
  if (w_asc != NULL)
  {
    if (!quietmode)
      show_world(curr_generation, stdout);
    if ((worldmode != 0) && (world_file != NULL))
      show_world(curr_generation, world_file);
  }
}


void draw_cell(unsigned long x, unsigned long y)
{
  if (w_asc != NULL)
  {
    if (world[x][y].plant_no == -1)
    {
      w_asc[x][y] = ' ';
      /* printf("cell erased at (%lu/%lu)\n", x, y); */
    }
    else
    {
      if (plant[world[x][y].plant_no].cell[world[x][y].cell_no].energy == 0)
      {
        /* w_asc[x][y] = 'o'; */
        w_asc[x][y] = 'a' + world[x][y].plant_no % 26;
        /* printf("energyless cell drawn at (%lu/%lu)\n", x, y); */
      }
      else
      {
        /* w_asc[x][y] = '*'; */
        w_asc[x][y] = 'A' + world[x][y].plant_no % 26;
        /* printf("energyrich cell drawn at (%lu/%lu)\n", x, y); */
      }
    }
  }
}


void poll_user_interface(void)
{
}


void clear_world_display(void)
{
}


void show_world_display(void)
{
}


void init_world_display(void)
{
  unsigned long x;

  w_asc = (char **) malloc(world_width * sizeof (char *));
  if (w_asc == NULL)
  {
    do_error("init_world_display failed: Out of memory");
  }
  else
  {
    w_asc[0] = (char *) malloc(world_height * world_width * sizeof(char));
    if (w_asc[0] == NULL)
    {
      free((void *) w_asc);
      w_asc = (char **) NULL;
      do_error("init_world_display failed: Out of memory");
    }
    else
    {
      for (x = 1; x < world_width; x++)
        w_asc[x] = w_asc[0] + x * world_height;
    }
  }
}


void close_world_display(void)
{
  if (w_asc != NULL)
  {
    free(w_asc[0]);
    free(w_asc);
    w_asc = (char **) NULL;
  }
}

