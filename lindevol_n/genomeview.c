/* $Id: genomeview.c,v 1.2 2000/01/30 03:10:59 kim Exp $ */
/*
 * $Log: genomeview.c,v $
 * Revision 1.2  2000/01/30 03:10:59  kim
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

#include "genomelib.h"
#include "lndtypes.h"
#include "lnddispl.h"
#include "lnderror.h"
#include "lndsignl.h"
#include "lndlib.h"
#include "lnd6.h"

#include "gnlib.h"

#include "lndglobl.h"

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif

#define DO_NEWSIM         1
#define DO_RESUMESIM         2
#define DO_EXIT         3
#define DO_TESTSIM         4

char buf[MAX_SLEN];

char par_fname[MAX_SLEN];       /* name of parameter file */
char g_fname[MAX_SLEN] = "\0";  /* name of genome file */

char quietmode = 0;
char worldmode = 0;

char finish_flag = 0;
char save_data = 0;

char simname[MAX_SLEN];

long psize_init;               /* initial population size */
long psize;                    /* current population size */
double m_replacement;          /* replacement rate */
double m_insertion;            /* insertion rate */
double m_deletion;             /* deletion rate */

long world_width, world_height;
long glen_init;                /* length of randomly initialized genomes */
double p_random_death;         /* probability for random death */
double rdeath_f_energy;        /* total energy dependent factor for random death */
double rdeath_f_numcells;      /* total number of cells dependent factor for random death */
double leanover_penalty;       /* probability of dying per unit of cell imbalance */
                               /* and time step */
long seedprod_threshold;       /* min. number of cells required for a plant */
                               /* to be able to reproduce */

GSYS_PARAMETERS gsys_parameters;

long start_generation;         /* the generation at which the simulation should start */
long num_generations;          /* number of generations to simulate */
long dmt_savefreq;             /* distance matrix is saved every so many generations */
long ddistr_savefreq;          /* frequency of saving distance distributions */
long phyltest_savefreq;        /* genomes are dumped every so many generations */
long bp_ceiling;               /* ceiling for B&P evolutionary activity analysis */

int random_seed = 1;           /* the random seed */

long generation;               /* the current generation */

long min_num_cells, max_num_cells;
double average_num_cells;
long min_cellular_energy, max_cellular_energy;
double average_cellular_energy;
double average_genome_length;
long min_genome_length, max_genome_length;
double average_num_used;
long min_num_used, max_num_used;
long num_species;
long mainspec_distance;
double distance_entropy = 0.0;
double distance_entropy_rel = 0.0;
long num_seeds, num_local_seeds, num_flying_seeds;
long num_new_plants;
long num_divisions, num_new_cells;
long num_attacks, num_deaths;
long min_age, max_age;
double average_age;
long num_mutminus, num_mutplus;
long min_num_mutations, max_num_mutations;
long num_unmutated_genomes;
long min_num_genes, max_num_genes;
double average_num_genes;
double genetic_diversity, rel_genetic_diversity;

GN_TREE gntree;

PLANT         **plant = (PLANT **) NULL;        /* array of pointers to plants */
GENOME          old_main_species;               /* the main species of previous time step */

SPECIES_D      *species_d = (SPECIES_D *) NULL; /* array of species descriptors (for counting species etc) */

LATTICE_SITE  **world = (LATTICE_SITE **) NULL; /* the array containing the plant cell locations */
long  *r_index = (long *) NULL;                 /* general purpose index array (for sorting, randomizing etc.) */
long  *pl_index = (long *) NULL;                /* plant index array for randomly shuffling the order of growth processing */
long  *tmp_index = (long *) NULL;               /* this index array is not initialized upon creation
                                                   and may be written to (in order to be used only partly) */

const int x_offset[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
const int y_offset[8] = {-1, -1, -1, 0, 0, 1, 1, 1};

long node_number;


int interpret_commandline(int argc, char *argv[])
{
  int i, ret_code = 0;

  par_fname[0] = '\0';
  for (i = 1; i < argc; i++)
  {
    if (!strcmp(argv[i], "-t"))
    {
      ret_code = DO_TESTSIM;
    }
    if (!strcmp(argv[i], "-q"))
      quietmode = 1;
    if (!strcmp(argv[i], "-w"))
      worldmode = 1;
    if (!strcmp(argv[i], "-f"))
    {
      if (i + 1 < argc)
      {
        i++;
        strncpy(par_fname, argv[i], MAX_SLEN);
        ret_code = DO_NEWSIM;
      }
    }
    if (!strcmp(argv[i], "-r"))
    {
      if (i + 1 < argc)
      {
        i++;
        strncpy(simname, argv[i], MAX_SLEN);
        ret_code = DO_RESUMESIM;
      }
    }
    if (!strcmp(argv[i], "-g"))
    {
      if (i + 1 < argc)
      {
        i++;
        strncpy(g_fname, argv[i], MAX_SLEN);
      }
    }
    if (!strcmp(argv[i], "-h"))
    {
      printf("\nLindEvol commandline options:\n\n");
      printf("-q        : start running in quiet mode\n");
      printf("-t        : run test simulation (with builtin test parameter set)\n");
      printf("-w        : create a world file\n");
      printf("-f <file> : read control parameters from <file>\n");
      printf("-r <name> : resume simulation <name>\n");
      printf("-g <file> : start with genomes read from <file>\n");
      printf("\nSignals:\n\n");
      printf("HUP  : write world state into *.dat file after current generation\n");
      printf("TERM : write *.dat file and terminate after current generation\n");
      printf("USR1 : toggle quiet mode\n");
      exit (0);
    }
  }
  return (ret_code);
}


void fprint_dnachars(FILE *f, unsigned char d)
{
  int i;
  char dna[5];

  for (i = 0; i < 4; i++)
  {
    switch (d & 0x03)
    {
    case 0:
      dna[i] = 'a';
      break;
    case 1:
      dna[i] = 'c';
      break;
    case 2:
      dna[i] = 'g';
      break;
    case 3:
      dna[i] = 't';
      break;
    }
    d >>= 2;
  }
  dna[4] = '\0';
  fprintf(f, dna);
}


void print_plant(FILE *f, long plant_no)
{
  long j, gpos, gene_no;
  unsigned char gene_output;
  unsigned long state_bits = 0, valid_bits = 0;
  unsigned long regulator_bit;
  long activity;

  fprintf(f, "** plant #%1ld, %ld cells, age = %ld ***n\n", plant_no, plant[plant_no]->num_cells, plant[plant_no]->age);
  fprintf(f, "raw listing:\n");
  gpos = 0;
  gene_no = 0;
  while (gpos < plant[plant_no]->genome.length)
  {
    if (l4_promoter(plant[plant_no]->genome.g[gpos]))
    {
      activity = -1;
      gpos++;
      if (gpos == plant[plant_no]->genome.length)
        break;
      while ((gpos < plant[plant_no]->genome.length) && !l4_promoter(plant[plant_no]->genome.g[gpos]))
      {
        if (l4_terminator(plant[plant_no]->genome.g[gpos]))
          break;
        regulator_bit = 1 << (plant[plant_no]->genome.g[gpos] & 0x7);
        valid_bits |= regulator_bit;
        if (plant[plant_no]->genome.g[gpos] & 0x20)
          state_bits |= regulator_bit;
        gpos++;
      }
      if (gpos == plant[plant_no]->genome.length)
        break;
      if (l4_terminator(plant[plant_no]->genome.g[gpos]))
      {
        gene_output = plant[plant_no]->genome.g[gpos];
        activity = gene_activity(gene_output);
        gpos++;
      }
      if (plant[plant_no]->genome.usg_count[gene_no] > 0)
        fprintf(f, "* ");
      else
        fprintf(f, "  ");
      fprintf(f, "%4ld (%4ld): ", gene_no, gpos);
      for (j = 31; j >= 0; j--)
      {
        if ((valid_bits & (1 << j)))
        {
          if ((state_bits & (1 << j)))
            fprintf(f, "1");
          else
            fprintf(f, "0");
        }
        else
        {
          fprintf(f, "-");
        }
      }
      fprintf(f, " -> ");
      if (activity == -1)
        fprintf(f, "no action    ");
      else if (activity == LND_DIVIDE)
        fprintf(f, "divide %1d     ", gene_output & 0x07);
      else if (activity == LND_FLYINGSEED)
        fprintf(f, "flying seed  ");
      else if (activity == LND_LOCALSEED)
        fprintf(f, "local seed   ");
      else if (activity == LND_MUTMINUS)
        fprintf(f, "mut-         ");
      else if (activity == LND_MUTPLUS)
        fprintf(f, "mut+         ");
      if (plant[plant_no]->genome.flags & GNM_BP)
        fprintf(f, "[%4d / %6d]\n", plant[plant_no]->genome.usg_count[gene_no], plant[plant_no]->genome.bp_count[gene_no]);
      else
        fprintf(f, "[%4d / ------]\n", plant[plant_no]->genome.usg_count[gene_no]);
      gene_no++;
    }
    gpos++;
  }
  fprintf(f, "\n");
}


int main(int argc, char *argv[])
{
  FILE *outfile;
  long plant_no;

  if (argc > 1)
  {
    strncpy(simname, argv[1], MAX_SLEN);
  }
  else
  {
    input_str("name of lindevol-4 simulation run: ", simname, MAX_SLEN);
  }
  if (argc > 2)
  {
    if ((outfile = fopen(argv[2], "w")) == NULL)
    {
      fprintf(stderr, "failed to open \"%s\" for output\n", argv[2]);
      return (-1);
    }
  }
  else
  {
    outfile = stdout;
  }
  if (load_savefile() != 0)
  {
    fprintf(stderr, "failed to read data for \"%s\n", simname);
    return (-1);
  }
  calculate_activity_codes(&gsys_parameters);
  for (plant_no = 0; plant_no < world_width; plant_no++)
  {
    if (plant[plant_no] != NULL)
    {
      print_plant(outfile, plant_no);
    }
  }
  return (0);
}

