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

#include "jklib.h"
#include "gnlib.h"

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif

#define DO_NEWSIM     1
#define DO_RESUMESIM  2
#define DO_EXIT       3
#define DO_TESTSIM    4

#define LND_MAIN
#include "lndglobl.h"


long node_number;

int main_menu(void)
{
  int ret_code;

  printf("Please choose\n");
#ifdef DO_NEWSIM
  printf("%1d - start new simulation\n", DO_NEWSIM);
#endif
#ifdef DO_RESUMESIM
  printf("%1d - resume a simulation\n", DO_RESUMESIM);
#endif
#ifdef DO_EXIT
  printf("%1d - exit program\n", DO_EXIT);
#endif
#ifdef DO_TESTSIM
  printf("%1d - test\n", DO_TESTSIM);
#endif
  fgets(buf, MAX_SLEN, stdin);
  if ((ret_code = strtol(buf, (char **) NULL, 10)) == DO_RESUMESIM)
  input_str("name of simulation to be resumed: ", simname, MAX_SLEN - 4);
  return (ret_code);
}


int get_control_parameters(void)
{
  (void) input_str("name of simulation run: ", simname, MAX_SLEN - 4);
  psize_init = input_long("initial population size: ");
  m_replacement =   input_dbl("replacement rate: ");
  m_insertion = input_dbl("insertion rate  : ");
  m_deletion = input_dbl("deletion rate   : ");
  m_deletion = input_dbl("gene duplication rate   : ");
  m_factor = input_dbl("mutation factor : ");

  world_width = input_long("world width: ");
  world_height = input_long("world height: ");
  world_soil = input_long("soil level: ");
  glen_init = input_long("initial genome length (in genes): ");
  p_random_death = input_dbl("probability for random death: ");
  rdeath_f_energy = input_dbl("random death probability factor per energy unit: ");
  rdeath_f_numcells = input_dbl("random death probability factor per energy unit: ");
  leanover_penalty = input_dbl("random death probability factor per energy unit: ");
  seedprod_threshold = input_long("threshold # of cells for seed production: ");
  nutrient_init = input_long("initial amount of nutrient in world: ");
  diffusion_rate = input_dbl("diffusion rate: ");
  decomposition_rate = input_dbl("decomposition rate: ");

  gsys_parameters.num_divide = input_long("number of code bytes for divide: ");
  gsys_parameters.num_flyingseed = input_long("number of code bytes for flying seed: ");
  gsys_parameters.num_localseed = input_long("number of code bytes for local seed: ");
  gsys_parameters.num_mutminus = input_long("number of code bytes for mut-: ");
  gsys_parameters.num_mutplus = input_long("number of code bytes for mut+: ");
  gsys_parameters.num_to_epool = input_long("number of code bytes for to_epool: ");
  gsys_parameters.num_to_npool = input_long("number of code bytes for to_npool: ");
  gsys_parameters.num_from_epool = input_long("number of code bytes for from_epool: ");
  gsys_parameters.num_from_npool = input_long("number of code bytes for from_npool: ");
  gsys_parameters.num_statebit = input_long("number of code bytes for set statebit: ");

  random_seed = input_long("random seed: ");

  num_generations = input_long("number of generations to simulate: ");
  dmt_savefreq = input_long("distance matrix every (0 for none): ");
  phyltest_savefreq = input_long("phylogenetic test sets every (0 for none): ");
  phyltest_samplesize = input_long("phylogenetic sample size: ");
  ddistr_savefreq = input_long("distance distribution saved every (0 for none): ");
  ddistr_samplesize = input_long("sample size for distance distribution analysis: ");
  bp_ceiling = input_long("ceiling for B&P analysis (0 for none): ");
  soil_savefreq = input_long("soil statistics saved every (0 for none): ");

  return (1);
}


int get_test_parameters(void)
{
  strcpy(simname, "test");
  psize_init = 50;
  m_replacement = 0.03;
  m_insertion = 0.01;
  m_deletion = 0.01;
  m_duplication = 0.0;
  m_factor = 1.0;

  world_width = 50;
  world_height = 20;
  world_soil = 10;
  glen_init = 50;
  p_random_death = 0.01;
  rdeath_f_energy = 0.95;
  rdeath_f_numcells = 1.0;
  leanover_penalty = 0.0;
  seedprod_threshold = 1;
  nutrient_init = 400;
  diffusion_rate = 1.0;
  decomposition_rate = 1.0;

  gsys_parameters.num_divide = 32;
  gsys_parameters.num_flyingseed = 2;
  gsys_parameters.num_localseed = 2;
  gsys_parameters.num_mutminus = 2;
  gsys_parameters.num_mutplus = 2;
  gsys_parameters.num_to_epool = 2;
  gsys_parameters.num_to_npool = 2;
  gsys_parameters.num_from_epool = 2;
  gsys_parameters.num_from_npool = 2;
  gsys_parameters.num_statebit = 16;

  random_seed = 1;

  dmt_savefreq = 0;
  phyltest_savefreq = 0;
  phyltest_samplesize = LONG_MAX;
  ddistr_savefreq = 0;
  ddistr_samplesize = LONG_MAX;
  bp_ceiling = 0;
  ddistr_savefreq = 20;
  num_generations = 1000;

  return (1);
}


char *identify_parameter(const char *line, const char *pname)
{
  size_t i;
  char *p;

  if ((p = strstr(line, pname)) == NULL)
    return (NULL);
  p += strlen(pname);
  i = strspn(p, " :=");
  if (i == 0)
    return (NULL);
  else
    return (p + i);
}


#define par_id_long(pname, pvar) \
  if ((s = identify_parameter(line, pname )) != NULL) \
  { \
    pvar = strtol(s, (char **) NULL, 10); \
    return (0); \
  }

#define par_id_double(pname, pvar) \
  if ((s = identify_parameter(line, pname )) != NULL) \
  { \
    pvar = strtod(s, (char **) NULL); \
    return (0); \
  }


int interpret_parameter(const char *line)
{
  char *s;
  int i;

  if ((s = identify_parameter(line, "simname")) != NULL)
  {
    i = 0;
    while (!iscntrl(s[i]))
    {
      simname[i] = s[i];
      i++;
    }
    simname[i] = '\0';
    return (0);
  }
  if ((s = identify_parameter(line, "g_fname")) != NULL)
  {
    i = 0;
    while (!iscntrl(s[i]))
    {
      g_fname[i] = s[i];
      i++;
    }
    g_fname[i] = '\0';
    return (0);
  }
  par_id_long("psize_init", psize_init);
  par_id_double("m_replacement", m_replacement);
  par_id_double("m_insertion", m_insertion);
  par_id_double("m_deletion", m_deletion);
  par_id_double("m_duplication", m_duplication);
  par_id_double("m_factor", m_factor);
  par_id_long("world_width", world_width);
  par_id_long("world_height", world_height);
  par_id_long("world_soil", world_soil);
  par_id_long("glen_init", glen_init);
  par_id_double("p_random_death", p_random_death);
  par_id_double("rdeath_f_energy", rdeath_f_energy);
  par_id_double("rdeath_f_numcells", rdeath_f_numcells);
  par_id_double("leanover_penalty", leanover_penalty);
  par_id_long("seedprod_threshold", seedprod_threshold);
  par_id_long("nutrient_init", nutrient_init);
  par_id_double("diffusion_rate", diffusion_rate);
  par_id_double("decomposition_rate", decomposition_rate);

  par_id_long("num_divide", gsys_parameters.num_divide);
  par_id_long("num_flyingseed", gsys_parameters.num_flyingseed);
  par_id_long("num_localseed", gsys_parameters.num_localseed);
  par_id_long("num_mutminus", gsys_parameters.num_mutminus);
  par_id_long("num_mutplus", gsys_parameters.num_mutplus);
  par_id_long("num_to_epool", gsys_parameters.num_to_epool);
  par_id_long("num_to_npool", gsys_parameters.num_to_npool);
  par_id_long("num_from_epool", gsys_parameters.num_from_epool);
  par_id_long("num_from_npool", gsys_parameters.num_from_npool);
  par_id_long("num_statebit", gsys_parameters.num_statebit);

  par_id_long("random_seed", random_seed);

  par_id_long("num_generations", num_generations);
  par_id_long("dmt_savefreq", dmt_savefreq);
  par_id_long("phyltest_savefreq", phyltest_savefreq);
  par_id_long("phyltest_samplesize", phyltest_samplesize);
  par_id_long("ddistr_savefreq", ddistr_savefreq);
  par_id_long("ddistr_samplesize", ddistr_samplesize);
  par_id_long("bp_ceiling", bp_ceiling);
  par_id_long("soil_savefreq", soil_savefreq);

  return (-1);
}


int load_control_parameters(const char *par_fname)
{
  FILE *f;
  char  errmsg[MAX_SLEN + 15];
  int   lineno = 0;

  if ((f = fopen(par_fname, "r")) == NULL)
  {
    do_error("Cannot open parameter file");
    return (-1);
  }
  get_test_parameters();
  num_generations = 0;
  while (!feof(f))
  {
    fgets(buf, MAX_SLEN, f);
    lineno++;
    if (buf[0] != '#')
    {
      if (strlen(buf) > 1)
      {
        if (interpret_parameter(buf) != 0)
        {
          fclose(f);
          sprintf(errmsg, "cannot understand line %d: %s\n", lineno, buf);
          do_error(errmsg);
          return (lineno);
        }
      }
    }
  }
  fclose(f);
  if (num_generations == 0)
  {
    num_generations = input_long("number of generations to simulate: ");
  }
  return (0);
}


int init_soil(long num_nutrient)
{
  long i;

  fprintf(stderr, "init_soil: introducing %ld nutrient units\n", num_nutrient);
  for (i = 0; i < num_nutrient; i++)
    add_nutrient_random();
  return (0);
}


/*
 * init_plants() creates the initial population of plants. They
 * are created randomly scattered on the world, each one with a
 * randomized genome of length glen_init and consisting of one
 * energyless cell.
 */

int init_plants(const char *g_fname, long num_plants, long glen)
{
  long   i, j, pos, plant_no;
  GENOME genome;
  FILE  *f;
  int ch, cl;
  unsigned long flags;
  long num_nutrient;

  num_nutrient = nutrient_init;
  psize = 0;
  generation = -1;
  flags = GNM_USG;
  if (bp_ceiling)
    flags |= GNM_BP;
  gn_init_tree(&gntree);
  if ((g_fname == NULL) || !strlen(g_fname))
  {
    if (alloc_genome(&genome, glen, 0, flags) < 0)
    {
      do_error("failed to allocate dummy ancestor genome");
      free_genome(&genome);
      return (-1);
    }
    for (i = 0; i < num_plants; i++)
    {
      for (j = 0; j < glen; j++)
        genome.g[j] = (unsigned char) lnd_random(256);
      if (resize_genome(&genome, glen, l4_num_genes(&genome)) < 0)
      {
        free_genome(&genome);
        return (-1);
      }
      if (flags & GNM_USG)
      {
        for (j = 0; j < genome.num_genes; j++)
          genome.usg_count[j] = 0;
      }
      if (flags & GNM_BP)
      {
        for (j = 0; j < genome.num_genes; j++)
          genome.bp_count[j] = 0;
      }
      pos = random_free_position();
      plant_no = create_plant(pos, -1, &genome);
      if (plant_no > -1)
        num_nutrient--;
      else
        do_error("init_plants: failed to create plant");
      if (phyltest_savefreq > 0)
      {
        if (gn_new_treenode(&gntree, NULL, -1, NULL, NULL, &(plant[plant_no]->node)) < 0)
        {
          do_error("init_plants: failed to create new genealogy node");
          free_genome(&genome);
          return (-1);
        }
      }
    }
  }
  else
  {
    if ((f = fopen(g_fname, "r")) == NULL)
    {
      do_error("init_plants: failed to open genome file");
      return (-1);
    }
    while (!feof(f) && !ferror(f))
    {
      if (fgets(buf, MAX_SLEN, f) == NULL)
      {
        break;
      }
      glen = strtol(buf, (char **) NULL, 10);
      if (alloc_genome(&genome, glen, 0, flags) < 0)
      {
        free_genome(&genome);
        return (-1);
      }
      for (j = 0; !iscntrl(ch = fgetc(f)); j++)
      {
        if (j >= genome.length)
        {
          do_error("init_plants: skipping characters while reading genome");
          fgets(buf, MAX_SLEN, f);
          break;
        }
        cl = fgetc(f);
        if ((ch == EOF) || (cl == EOF) || (iscntrl(cl)))
        {
          do_error("init_plants: unexpected EOF or control character in control character");
          return (-1);
        }
        genome.g[j] = (hexval(ch) << 4) | hexval(cl);
      }
      if (resize_genome(&genome, glen, l4_num_genes(&genome)) < 0)
      {
        free_genome(&genome);
        return (-1);
      }
      if (genome.flags & GNM_USG)
      {
        for (j = 0; j < genome.num_genes; j++)
          genome.usg_count[j] = 0;
      }
      if (genome.flags & GNM_BP)
      {
        for (j = 0; j < genome.num_genes; j++)
          genome.bp_count[j] = 0;
      }
      pos = random_free_position();
      plant_no = create_plant(pos, -1, &genome);
      if (plant_no > -1)
        num_nutrient--;
      else
        do_error("init_plants: failed to create plant");
      if (phyltest_savefreq > 0)
      {
        if (gn_new_treenode(&gntree, NULL, -1, NULL, NULL, &(plant[plant_no]->node)) < 0)
        {
          do_error("create_plant: failed to create new genealogy node");
          return (-1);
        }
      }
      free_genome(&genome);
    }
  }
  free_genome(&genome);
  old_main_species.flags = 0;
  old_main_species.length = 0;
  old_main_species.num_genes = 0;
  old_main_species.g = NULL;
  old_main_species.usg_count = NULL;
  old_main_species.bp_count = NULL;
  old_main_species.mut_flag = 0;
  old_main_species.num_mutations = 0;
  init_soil(num_nutrient);
  return (0);
}


/*
 * Reset the counters for some events to zero. Called before each generation
 */

void zero_counters(void)
{
  long i;

  num_seeds = 0;
  num_local_seeds = 0;
  num_flying_seeds = 0;
  num_new_plants = 0;
  num_divisions = 0;
  num_new_cells = 0;
  num_attacks = 0;
  num_deaths = 0;
  num_mutminus = 0;
  num_mutplus = 0;
  num_unmutated_genomes = 0;
  min_deathprob = 1.0;
  max_deathprob = 0.0;
  sum_deathprob = 0.0;
  num_grownplants = 0;
  num_to_epool = 0;
  num_to_npool = 0;
  num_from_epool = 0;
  num_from_npool = 0;
  for (i = 0; i < NUM_STATEBITS; i++)
  {
    num_bitchecks[i] = 0;
    num_bitsets[i] = 0;
  }
}


/*
PROCEDURE auswerten
  ' *****************************************************************
  ' * This procedure calculates some global values of a generation: *
  ' * The minimal, maximal and average fitness,                     *
  ' * the minimal, maximal and average genome length in genes,      *
  ' * the minimal, maximal and average number of used genes.        *
  ' *****************************************************************
*/

void statistics(void)
{
  long num_cells_sum = 0;
  long cellular_energy_sum = 0;
  long cellular_nutrient_sum = 0;
  long glen_sum = 0;
  long num_used_sum = 0;
  long age_sum = 0;
  long num_genes_sum = 0;
  long num_used;
  long mutcounter_sum = 0;
  long energy_pool_sum = 0, nutrient_pool_sum = 0;
  long specificity_sum = 0, total_num_genes = 0, tng_check = 0, specificity;
  STATE_SPEC statespec;
  long i, j;

  min_num_cells = LONG_MAX;
  max_num_cells = 0;
  min_cellular_energy = LONG_MAX;
  max_cellular_energy = 0;
  min_cellular_nutrient = LONG_MAX;
  max_cellular_nutrient = 0;
  min_genome_length = LONG_MAX;
  max_genome_length = 0;
  min_num_used = LONG_MAX;
  max_num_used = 0;
  min_age = LONG_MAX;
  max_age = 0;
  min_num_genes = LONG_MAX;
  max_num_genes = 0;
  min_mutcounter = LONG_MAX;
  max_mutcounter = LONG_MIN;
  min_specificity = LONG_MAX;
  max_specificity = 0;
  min_energy_pool = LONG_MAX;
  max_energy_pool = 0;
  min_nutrient_pool = LONG_MAX;
  max_nutrient_pool = 0;
  num_free_nutrient = 0;
  num_organic_nutrient = 0;
  num_biomass_nutrient = 0;
  for (i = 0; i < world_width; i++)
  {
    if (plant[i] != NULL)
    {
      num_used = 0;
      for (j = 0; j < plant[i]->genome.num_genes; j++)
      {
        if (plant[i]->genome.usg_count[j] > 0)
          num_used++;
      }
      for (j = 0; j < plant[i]->genome.length; j++)
      {
        if (l4_promoter(plant[i]->genome.g[j]))
        {
          compute_statespec(&(plant[i]->genome), j + 1, &statespec, NULL);
          specificity = num_setbits(statespec.valid_bits);
          min_specificity = (min_specificity < specificity) ? min_specificity : specificity;
          max_specificity = (max_specificity > specificity) ? max_specificity : specificity;
          specificity_sum += specificity;
          tng_check++;
        }
      }
      total_num_genes += plant[i]->genome.num_genes;
      if (total_num_genes != tng_check)
        fprintf(stderr, "statistics: detected numgenes inconsistency: plant #%ld, n=%ld, check=%ld\n", i, total_num_genes, tng_check);
      num_cells_sum += plant[i]->num_cells;
      min_num_cells = (min_num_cells < plant[i]->num_cells) ? min_num_cells : plant[i]->num_cells;
      max_num_cells = (max_num_cells > plant[i]->num_cells) ? max_num_cells : plant[i]->num_cells;
      cellular_energy_sum += plant[i]->cellular_energy;
      min_cellular_energy = (min_cellular_energy < plant[i]->cellular_energy) ? min_cellular_energy : plant[i]->cellular_energy;
      max_cellular_energy = (max_cellular_energy > plant[i]->cellular_energy) ? max_cellular_energy : plant[i]->cellular_energy;
      cellular_nutrient_sum += plant[i]->cellular_nutrient;
      min_cellular_nutrient = (min_cellular_nutrient < plant[i]->cellular_nutrient) ? min_cellular_nutrient : plant[i]->cellular_nutrient;
      max_cellular_nutrient = (max_cellular_nutrient > plant[i]->cellular_nutrient) ? max_cellular_nutrient : plant[i]->cellular_nutrient;
      energy_pool_sum += plant[i]->energy_pool;
      min_energy_pool = (min_energy_pool < plant[i]->energy_pool) ? min_energy_pool : plant[i]->energy_pool;
      max_energy_pool = (max_energy_pool > plant[i]->energy_pool) ? max_energy_pool : plant[i]->energy_pool;
      nutrient_pool_sum += plant[i]->nutrient_pool;
      min_nutrient_pool = (min_nutrient_pool < plant[i]->nutrient_pool) ? min_nutrient_pool : plant[i]->nutrient_pool;
      max_nutrient_pool = (max_nutrient_pool > plant[i]->nutrient_pool) ? max_nutrient_pool : plant[i]->nutrient_pool;
      glen_sum += plant[i]->genome.length;
      min_genome_length = (min_genome_length < plant[i]->genome.length) ? min_genome_length : plant[i]->genome.length;
      max_genome_length = (max_genome_length > plant[i]->genome.length) ? max_genome_length : plant[i]->genome.length;
      num_used_sum += num_used;
      min_num_used = (min_num_used < num_used) ? min_num_used : num_used;
      max_num_used = (max_num_used > num_used) ? max_num_used : num_used;
      age_sum += plant[i]->age;
      min_age = (min_age < plant[i]->age) ? min_age : plant[i]->age;
      max_age = (max_age > plant[i]->age) ? max_age : plant[i]->age;
      min_num_genes = (min_num_genes < plant[i]->genome.num_genes)
              ? min_num_genes : plant[i]->genome.num_genes;
      max_num_genes = (max_num_genes > plant[i]->genome.num_genes)
              ? max_num_genes : plant[i]->genome.num_genes;
      num_genes_sum += plant[i]->genome.num_genes;
      mutcounter_sum += plant[i]->genome.mut_flag;
      min_mutcounter = (min_mutcounter < plant[i]->genome.mut_flag)
                          ? min_mutcounter : plant[i]->genome.mut_flag;
      max_mutcounter = (max_mutcounter > plant[i]->genome.mut_flag)
                          ? max_mutcounter : plant[i]->genome.mut_flag;
      num_biomass_nutrient += plant[i]->num_cells + plant[i]->cellular_nutrient + plant[i]->nutrient_pool;
    }
  }
  if (psize > 0)
  {
    average_num_cells = ((double) num_cells_sum) / psize;
    average_cellular_energy = ((double) cellular_energy_sum) / psize;
    average_cellular_nutrient = ((double) cellular_nutrient_sum) / psize;
    average_energy_pool = ((double) energy_pool_sum) / psize;
    average_nutrient_pool = ((double) nutrient_pool_sum) / psize;
    average_genome_length = ((double) glen_sum) / psize;
    average_num_used = ((double) num_used_sum) / psize;
    average_age = ((double) age_sum) / psize;
    average_mutcounter = (double) mutcounter_sum / psize;
    average_num_genes = ((double) num_genes_sum) / psize;
  }
  else
  {
    min_num_cells = 0;
    min_cellular_energy = 0;
    min_genome_length = 0;
    min_num_used = 0;
    min_num_genes = 0;
    min_age = 0;
    min_num_genes = 0;
    min_mutcounter = 0;
    max_mutcounter = 0;
    min_cellular_energy = 0;
    min_cellular_nutrient = 0;
    min_energy_pool = 0;
    min_nutrient_pool = 0;
    average_num_cells = 0.0;
    average_cellular_energy = 0.0;
    average_genome_length = 0.0;
    average_num_used = 0.0;
    average_age = 0.0;
    average_mutcounter = 0.0;
    average_num_genes = 0.0;
    average_cellular_energy = 0.0;
    average_cellular_nutrient = 0.0;
    average_energy_pool = 0.0;
    average_nutrient_pool = 0.0;
  }
  if (total_num_genes > 0)
  {
    average_specificity = ((double) specificity_sum) / ((double) total_num_genes);
  }
  else
  {
    min_specificity = 0;
    average_specificity = 0.0;
  }
  for (i = 0; i < world_width; i++)
  {
    for (j = 0; j <= world_soil; j++)
    {
      num_organic_nutrient += world[i][j].organic_nutrient;
      if (world[i][j].nutrient)
        num_free_nutrient++;
    }
  }
}


void nutrient_check(void)
{
  long i, j, sum, pn;

  sum = 0;
  for (i = 0; i < world_width; i++)
  {
    for (j = 0; j <= world_soil; j++)
    {
      if (world[i][j].nutrient)
        sum++;
    }
  }
  for (i = 0; i < world_width; i++)
  {
    if (plant[i])
    {
      pn = 0;
      for (j = 0; j < plant[i]->num_cells; j++)
      {
        if (plant[i]->cell[j].nutrient)
          pn++;
      }
      if (pn != plant[i]->cellular_nutrient)
        fprintf(stderr, "generation %ld, plant %ld: cellular nutrient is %ld, actually is %ld\n",
                generation, i, plant[i]->cellular_nutrient, pn);
      sum += pn + plant[i]->nutrient_pool + plant[i]->num_cells;
    }
  }
  fprintf(stderr, "generation %5ld: total nutrient is %ld\n", generation, sum);
}


void soil_statistics(void)
{
  long i, j;

  for (i = 0; i < world_width; i++)
  {
    soil_profile_h[i] = 0;
    organic_profile_h[i] = 0;
    for (j = 0; j <= world_soil; j++)
    {
      organic_profile_h[i] += world[i][j].organic_nutrient;
      if (world[i][j].nutrient)
        soil_profile_h[i]++;
    }
  }
  for (j = 0; j <= world_soil; j++)
  {
    organic_profile_v[j] = 0;
    soil_profile_v[j] = 0;
    for (i = 0; i < world_width; i++)
    {
      organic_profile_v[j] += world[i][j].organic_nutrient;
      if (world[i][j].nutrient)
        soil_profile_v[j]++;
    }
  }
}


/*
PROCEDURE evolution
  ' **********************************************************************
  ' * This is the master simulation procedure. It simulates generation   *
  ' * until either the run is halted by the user or the generation limit *
  ' * is reached.                                                        *
  ' **********************************************************************
*/

void evolution(void)
{
  long finish_generation, save_generation = -1;
  long plant_no, s;

  finish_generation = start_generation + num_generations;
  init_world_display(simname);
  if (savetime_fname[0])
  {
    open_savetime_file(savetime_fname);
    save_generation = savetime_next(start_generation);
  }
  for (generation = start_generation; (generation < finish_generation) && (psize > 0) && !finish_flag; generation++)
  {
    /* fprintf(stderr, "starting generation %ld\n", generation); */
    count_species();           /* determine number of different species */
    if (bp_ceiling > 0)
    {
      write_bpe();             /* write B & P histogram for current generation */
    }
    node_number = 0;
    zero_counters();
    random_shuffle(world_width, pl_index);
    display_start(generation);    /* entertain the users ... ;-) */
    sunshine();
    nutrient_diffusion();
    for (plant_no = 0; plant_no < world_width; plant_no++)
    {
      if (plant[pl_index[plant_no]] != NULL)
      {
        plant_growth(pl_index[plant_no]);
      }
    }
    display_done();
    /* nutrient_check(); */
    statistics();              /* ... of fitness, genome length, gene usage */
    if ((ddistr_savefreq > 0) && ((generation % ddistr_savefreq) == 0))
    {
      s = prepare_sample(ddistr_samplesize, sample_index);
      write_dst(1, s, sample_index);       /* write the distance distribution */
    }
    if (dmt_savefreq > 0)
    {
      if ((generation % dmt_savefreq) == 0)
      {
        write_dmt();           /* write distance matrix (if desired) */
      }
    }
    if (phyltest_savefreq > 0)
    {
      if ((generation % phyltest_savefreq) == 0)
      {
        s = prepare_sample(phyltest_samplesize, sample_index);
        write_genomes(s, sample_index);
        write_jf(s, sample_index);             /* write phylogenetic trees */
      }
    }
    if ((soil_savefreq > 0) && ((generation % soil_savefreq) == 0))
      soil_statistics();
    /* fprintf(stderr, "writing protocol\n"); */
    write_pro();               /* write out the protocol for current generation */
    if ((generation == save_generation) || (save_data))
    {
      sprintf(buf, "%s-%07ld.dat", simname, generation);
      write_named_savefile(buf);
      if (generation == save_generation)
        save_generation = savetime_next(generation);
      save_data = 0;
    }
    if (pixmode)
      write_pixfile();
    /* printf("generation %ld: %ld plants, %ld genealogy nodes\n", generation, psize, gntree.num_nodes); */
#if defined(MEMDEBUG) && defined(FULL_MEMDEBUG)
    check_memdebugError();
    set_MemdebugOptions(c_Yes, /* general statistics */
                        c_Yes, /* alphabetical list */
                      c_No, /* not-free list */
                      c_No,  /* call sequence list */
                      c_Yes, /* spurious free list */
                      c_No,  /* print contents */
                      c_Yes, /* destroy contents */
                      0,     /* GenerateErrorCount */
                      0,     /* max. memory available */
                      "lndmem.dat", "lndmem.err");
    print_MemdebugStatistics();
#endif
  }
  close_savetime_file();
  write_savefile();
  free_genome(&old_main_species);
  close_world_display();
}


int interpret_commandline(int argc, char *argv[])
{
  int i, ret_code = 0;

  par_fname[0] = '\0';
  savetime_fname[0] = '\0';
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
    if (!strcmp(argv[i], "-p"))
      pixmode = 1;
    if (!strcmp(argv[i], "-f"))
    {
      if (i + 1 < argc)
      {
        i++;
        strncpy(par_fname, argv[i], MAX_SLEN);
        ret_code = DO_NEWSIM;
      }
    }
    if (!strcmp(argv[i], "-s"))
    {
      if (i + 1 < argc)
      {
        i++;
        strncpy(savetime_fname, argv[i], MAX_SLEN);
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
      printf("-p        : create a pixel file\n");
      printf("-f <file> : read control parameters from <file>\n");
      printf("-r <name> : resume simulation <name>\n");
      printf("-g <file> : start with genomes read from <file>\n");
      printf("-s <file> : read generations at which to save state from <file>\n");
      printf("\nSignals:\n\n");
      printf("HUP  : write world state into *.dat file after current generation\n");
      printf("TERM : write *.dat file and terminate after current generation\n");
      printf("USR1 : toggle quiet mode\n");
      exit (0);
    }
  }
  return (ret_code);
}


int main(int argc, char *argv[])
{
  int what_now;

  printf("%s, written by Jan T. Kim\n", SIMPRGNAME);
  printf("compiled %s %s\n", __DATE__, __TIME__);
  what_now = interpret_commandline(argc, argv);
  if (what_now == 0)
    what_now = main_menu();
  if (what_now == DO_NEWSIM)
  {
    if (strlen(par_fname))
    {
      if (load_control_parameters(par_fname))
      {
        return (EXIT_FAILURE);
      }
    }
    else
    {
      get_control_parameters();
    }
  }
  else if (what_now == DO_TESTSIM)
  {
    get_test_parameters();
    what_now = DO_NEWSIM;
  }
  if (what_now == DO_NEWSIM)
  {
    if (create_arrays() < 0)
    {
      clear_arrays();
      do_error("Out of memory error");
      return (EXIT_FAILURE);
    }
    if (open_data_files("w") < 0)
    {
      return (EXIT_FAILURE);
    }
    start_generation = 0;
    seed_lnd_random(random_seed);
    /* fprintf(stderr, "initializing\n"); */
    init_plants(g_fname, psize_init, glen_init);
  }
  else if (what_now == DO_RESUMESIM)
  {
    if (load_savefile() < 0)
    {
      return (EXIT_FAILURE);
    }
    if (open_data_files("a") < 0)
    {
      return (EXIT_FAILURE);
    }
    if (num_generations == 0)
    {
      num_generations = input_long("number of generations to simulate: ");
    }
  }
  if (what_now != DO_EXIT)
  {
    if (num_generations == 0)
    {
      num_generations = input_long("number of generations to simulate: ");
    }
    calculate_activity_codes(&gsys_parameters);
    init_signal_handling();
    evolution();
    close_data_files();
    if (worldmode)
    {
      close_world_file();
    }
    clear_arrays();
  }
  gn_free_tree(&gntree);

#ifdef MEMDEBUG
#  ifdef FULL_MEMDEBUG
  set_MemdebugOptions(c_Yes, /* general statistics */
                      c_Yes, /* alphabetical list */
                      c_Yes, /* not-free list */
                      c_No,  /* call sequence list */
                      c_Yes, /* spurious free list */
                      c_No,  /* print contents */
                      c_Yes, /* destroy contents */
                      0,     /* GenerateErrorCount */
                      0,     /* max. memory available */
                      "lndmem.dat", "lndmem.err");
#  endif
#endif

  return (0);
}

