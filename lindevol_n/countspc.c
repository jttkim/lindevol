/* $Id: countspc.c,v 1.2 2000/01/30 03:10:59 kim Exp $ */
/*
 * $Log: countspc.c,v $
 * Revision 1.2  2000/01/30 03:10:59  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

#include <math.h>
#include <stdlib.h>

#include "jklib.h"
#include "genomelib.h"
#include "lndvals.h"
#include "lndtypes.h"
#include "lndglobl.h"
#include "lnderror.h"
#include "lndlib.h"


/*
PROCEDURE arten_zaehlen
  '
  ' ***************************************************************************
  ' * This procedure determines the number of different genomes and the most  *
  ' * frequent genome. This is done by sorting a copy of the g$() array       *
  ' * alphabetically. Then the number of times where g1$(i%) <> g1$(i%)       *
  ' * is the number of different species. The number of copies of each        *
  ' * species is determined in the same loop by incrementing the number if    *
  ' * two adjacent genomes are identical and initializing a new number        *
  ' * otherwise. To determine the most frequent genome, the array of numbers  *
  ' * of copies of each species is sorted along with an array of pointers     *
  ' * pointing to the index of a genome of that species. Then, the first      *
  ' * pointer points to the genome for which the number of copies is largest, *
  ' * i.e. the most frequent genome.                                          *
  ' ***************************************************************************
*/

static int cmp_lexical(const void *plant1, const void *plant2)
{
  /* These obfuscated initializers evaluate to the addresses of the genome of the
     plant with the index pointed to by plant1 and plant2, respectively. plant1
     and plant2 have to point to longs in the array to be sorted. */

  const GENOME *g1 = &(plant[*((long *) plant1)]->genome);
  const GENOME *g2 = &(plant[*((long *) plant2)]->genome);
  long i, l;

  /* printf("comparing genomes %lu and %lu\n", g1, g2); */
  l = (g1->length < g2->length) ? g1->length : g2->length;
  for (i = 0; i < l; i++)
  {
    if (g1->g[i] < g2->g[i])
    {
      return (-1);
    }
    else if (g1->g[i] > g2->g[i])
    {
      return (1);
    }
  }
  return ((g1->length < g2->length) ? -1 : ((g1->length > g2->length) ? 1 : 0));
}

/* cmp_species_d() returns -1 if s1 is more frequent than s2 and 1 if s1 is less frequent
   than s2, in order to lead to descending qsorting of the species_d[] array */

static int cmp_species_d(const void *s1, const void *s2)
{
  long n1 = ((SPECIES_D *) s1)->num;
  long n2 = ((SPECIES_D *) s2)->num;

  return (n1 < n2) ? 1 : ((n1 > n2) ? -1 : 0);
}

void count_species(void)
{
  GENOME main_species;
  long i, j;
  double f;

  /* printf("counting species\n"); */
  num_species = 0;
  mainspec_distance = 0;
  genetic_diversity = 0.0;
  rel_genetic_diversity = 0.0;
  if (psize)
  {
    j = 0;
    for (i = 0; i < world_width; i++)
    {
      if (plant[i] != NULL)
      {
        tmp_index[j++] = i;
      }
    }
    qsort((void *) tmp_index, (size_t) psize, sizeof(long), &cmp_lexical);
    /* printf("sorting done\n"); */
    num_species = 0;
    species_d[0].num = 1;
    species_d[0].index = tmp_index[0];
    for (i = 1; i < psize; i++)
    {
      if (cmp_lexical((void *) &(tmp_index[i - 1]), (void *) &(tmp_index[i])) == 0)
      {
        species_d[num_species].num++;
      }
      else
      {
        num_species++;
        species_d[num_species].num = 1;
        species_d[num_species].index = tmp_index[i];
      }
    }
    num_species++;
    if (num_species > 1)
    {
      for (i = 0; i < num_species; i++)
      {
        f = (double) species_d[i].num / psize;
        genetic_diversity -= f * log(f);
      }
      rel_genetic_diversity = -genetic_diversity / log(1.0 / psize);
    }
    qsort((void *) species_d, (size_t) num_species, sizeof(SPECIES_D), &cmp_species_d);
    if (duplicate_genome(&main_species, &(plant[species_d[0].index]->genome), 0) < 0)
    {
      do_error("count_species: failed to allocate memory for main species");
      mainspec_distance = -1;
      free_genome(&old_main_species);
      return;
    }
    if (generation == 0)
    {
      mainspec_distance = 0;
    }
    else
    {
      mainspec_distance = edit_distance(main_species.length, (char *) main_species.g, old_main_species.length, (char *)old_main_species.g);
      free_genome(&old_main_species);
    }
    if (mainspec_distance == -1)
    {
      do_error("count_species: error while computing distance");
    }
    old_main_species = main_species;
  }
}

