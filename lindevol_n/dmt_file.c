/* $Id: dmt_file.c,v 1.2 2000/01/30 03:10:59 kim Exp $ */
/*
 * $Log: dmt_file.c,v $
 * Revision 1.2  2000/01/30 03:10:59  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

#include <stdio.h>

#include "jklib.h"
#include "lndvals.h"
#include "lndtypes.h"
#include "lndglobl.h"
#include "lndlibin.h"
#include "lndlib.h"


FILE *dmt_file = NULL;


int open_dmt_file(const char *mode)
{
  sprintf(buf, "%sb", mode);
  dmt_file = fopen(dmt_file_name, buf);
  if (dmt_file == NULL)
  {
    perror("open_dmt_file failed");
    return (-1);
  }
  if (*mode == 'w')
  {
    fprintf(dmt_file, "BINARY 16BIT\r\n");
    fprintf(dmt_file, "500\r\n");
  }
  return (0);
}


/*
PROCEDURE distanzmatrix_speichern
  ' ********************************************************************
  ' * Saves the matrix of edit distances between all different genomes *
  ' * in the population.                                               *
  ' ********************************************************************
*/

void write_dmt(void)
{
  long i, j;
  short d;

  fprintf(dmt_file, "%ld\r\n", psize);
  fprintf(dmt_file, "%ld\r\n", generation);
  j = 0;
  for (i = 0; i < world_width; i++)
  {
    if (plant[i] != NULL)
    {
      tmp_index[j++] = i;
    }
  }
  for (i = 0; i < psize - 1; i++)
  {
    for (j = i + 1; j < psize; j++)
    {
      d = (short) edit_distance(plant[tmp_index[i]]->genome.length, (char *) plant[tmp_index[i]]->genome.g,
                                plant[tmp_index[j]]->genome.length, (char *) plant[tmp_index[j]]->genome.g);
      fwrite_int16array(&d, 1, dmt_file);
    }
  }
  fflush(dmt_file);
}

