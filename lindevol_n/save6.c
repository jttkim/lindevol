/* $Id: save6.c,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: save6.c,v $
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lndvals.h"
#include "lndtypes.h"
#include "lndglobl.h"
#include "lnderror.h"
#include "lndlibin.h"
#include "lndlib.h"
#include "lnd6.h"

#include "genomelib.h"
#include "gnlib.h"


static void write_plant_cell(FILE *f, const PL_CELL *c)
{
  fprintf(f, "%ld\n", c->x);
  fprintf(f, "%ld\n", c->y);
  fprintf(f, "%d\n", c->energy);
  fprintf(f, "%d\n", c->nutrient);
  fprintf(f, "%lu\n", c->state);
  fprintf(f, "%lu\n", c->next_state);
}


static void write_gsys_parameters(FILE *f, const GSYS_PARAMETERS *gp)
{
  fprintf(f, "%ld\n", gp->num_divide);
  fprintf(f, "%ld\n", gp->num_flyingseed);
  fprintf(f, "%ld\n", gp->num_localseed);
  fprintf(f, "%ld\n", gp->num_mutminus);
  fprintf(f, "%ld\n", gp->num_mutplus);
  fprintf(f, "%ld\n", gp->num_to_epool);
  fprintf(f, "%ld\n", gp->num_to_npool);
  fprintf(f, "%ld\n", gp->num_from_epool);
  fprintf(f, "%ld\n", gp->num_from_npool);
  fprintf(f, "%ld\n", gp->num_statebit);
}


static void read_plant_cell(FILE *f, PL_CELL *c)
{
  char buf[MAX_SLEN];

  fgets(buf, MAX_SLEN, f);
  c->x = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  c->y = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  c->energy = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  c->nutrient = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  c->state = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  c->next_state = strtol(buf, (char **) NULL, 10);
}


static void read_gsys_parameters(FILE *f, GSYS_PARAMETERS *gp)
{
  fgets(buf, MAX_SLEN, f);
  gp->num_divide = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  gp->num_flyingseed = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  gp->num_localseed = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  gp->num_mutminus = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  gp->num_mutplus = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  gp->num_to_epool = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  gp->num_to_npool = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  gp->num_from_epool = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  gp->num_from_npool = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  gp->num_statebit = strtol(buf, (char **) NULL, 10);
}


int write_named_savefile(const char *save_file_name)
{
  long i, j;
  FILE *f;
  
  if ((f = fopen(save_file_name, "w")) == NULL)
  {
    return (-1);
  }
  fprintf(f, "Data dump of %s simulation, run %s, generation: %ld\n", SIMPRGNAME, simname, generation);
  fprintf(f, "%s\n", simname);
  fprintf(f, "%ld\n", generation);
  fprintf(f, "%ld\n", psize_init);
  fprintf(f, "%ld\n", psize);
  fprintf(f, "%1.12g\n", m_replacement);
  fprintf(f, "%1.12g\n", m_insertion);
  fprintf(f, "%1.12g\n", m_deletion);
  fprintf(f, "%1.12g\n", m_duplication);
  fprintf(f, "%1.12g\n", m_factor);
  fprintf(f, "%ld\n", world_width);
  fprintf(f, "%ld\n", world_height);
  fprintf(f, "%ld\n", world_soil);
  fprintf(f, "%ld\n", nutrient_init);
  fprintf(f, "%ld\n", nutrient_per_timestep);
  fprintf(f, "%1.12g\n", diffusion_rate);
  fprintf(f, "%1.12g\n", organic_diffusion_rate);
  fprintf(f, "%1.12g\n", decomposition_rate);
  fprintf(f, "%1.12g\n", p_random_death);
  fprintf(f, "%1.12g\n", rdeath_f_energy);
  fprintf(f, "%1.12g\n", rdeath_f_numcells);
  fprintf(f, "%1.12g\n", leanover_penalty);
  fprintf(f, "%ld\n", seedprod_threshold);
  fprintf(f, "%ld\n", dmt_savefreq);
  fprintf(f, "%ld\n", phyltest_savefreq);
  fprintf(f, "%ld\n", phyltest_samplesize);
  fprintf(f, "%ld\n", ddistr_savefreq);
  fprintf(f, "%ld\n", ddistr_samplesize);
  fprintf(f, "%ld\n", soil_savefreq);
  fprintf(f, "%ld\n", bp_ceiling);
  write_gsys_parameters(f, &gsys_parameters);
  fprintf(f, "%ld\n", random_seed);
  write_rndgenerator_state(f);
  if (phyltest_savefreq > 0)
  {
    /* gn_print_tree(&gntree); */
    gn_save_tree(&gntree, f);
  }
  write_genome(f, &old_main_species, 0);
  for (i = 0; i < world_width; i++)
  {
    fprintf(f, "%ld\n", r_index[i]);
    fprintf(f, "%ld\n", pl_index[i]);
  }
  fprintf(f, "%ld\n", nutrient_i);
  for (i = 0; i < (world_soil + 1) * world_width; i++)
    fprintf(f, "%ld\n", soil_index[i]);
  for (i = 0; i < world_width; i++)
  {
    if (plant[i] == NULL)
    {
      fprintf(f, "-1\n");
    }
    else
    {
      fprintf(f, "%ld\n", plant[i]->num_cells);
      fprintf(f, "%ld\n", plant[i]->age);
      fprintf(f, "%ld\n", plant[i]->energy_pool);
      fprintf(f, "%ld\n", plant[i]->nutrient_pool);
      if (phyltest_savefreq > 0)
        gn_save_id(&(plant[i]->node), f);
      for (j = 0; j < plant[i]->num_cells; j++)
      {
        write_plant_cell(f, &(plant[i]->cell[j]));
      }
      write_genome(f, &(plant[i]->genome), 0xffffffff);
    }
  }
  for (i = 0; i < world_width; i++)
  {
    for (j = 0; j <= world_soil; j++)
    {
      fprintf(f, "%u\n", world[i][j].nutrient);
      fprintf(f, "%ld\n", world[i][j].organic_nutrient);
    }
  }
  fclose(f);
  return (0);
}


int write_savefile(void)
{
  return (write_named_savefile(save_file_name));
}


int load_named_savefile(const char *save_file_name)
{
  FILE *f;
  char buf[MAX_SLEN + 1];
  long i, j, expected_psize;

  prepare_filenames();
  if ((f = fopen(save_file_name, "r")) == NULL)
  {
    do_error("load_savefile: failed to open save file");
    return (-1);
  }
  fgets(buf, MAX_SLEN, f);
  fgets(simname, MAX_SLEN, f);
  i = strlen(simname) - 1;
  if ((i >= 0) && (simname[i] == '\n'))
    simname[i] = '\0';
  prepare_filenames();
  fgets(buf, MAX_SLEN, f);
  start_generation = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  psize_init = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  expected_psize = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  m_replacement = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  m_insertion = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  m_deletion = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  m_duplication = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  m_factor = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  world_width = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  world_height = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  world_soil = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  nutrient_init = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  nutrient_per_timestep = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  diffusion_rate = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  organic_diffusion_rate = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  decomposition_rate = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  p_random_death = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  rdeath_f_energy = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  rdeath_f_numcells = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  leanover_penalty = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  seedprod_threshold = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  dmt_savefreq = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  phyltest_savefreq = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  phyltest_samplesize = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  ddistr_savefreq = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  ddistr_samplesize = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  soil_savefreq = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  bp_ceiling = strtol(buf, (char **) NULL, 10);
  read_gsys_parameters(f, &gsys_parameters);
  fgets(buf, MAX_SLEN, f);
  random_seed = strtol(buf, (char **) NULL, 10);
  read_rndgenerator_state(f);
  if (phyltest_savefreq > 0)
  {
    gn_read_tree(&gntree, f);
    /* gn_print_tree(&gntree); */
  }
  if (create_arrays() < 0)
  {
    do_error("load_savefile: out of memory");
    clear_arrays();
    fclose(f);
    return (-1);
  }
  if (read_genome(f, &old_main_species, 0) < 0)
  {
    do_error("load_savefile: error reading old main species");
    clear_arrays();
    fclose(f);
    return (-1);
  }
  for (i = 0; i < world_width; i++)
  {
    fgets(buf, MAX_SLEN, f);
    r_index[i] = strtol(buf, (char **) NULL, 10);
    fgets(buf, MAX_SLEN, f);
    pl_index[i] = strtol(buf, (char **) NULL, 10);
  }
  fgets(buf, MAX_SLEN, f);
  nutrient_i = strtol(buf, (char **) NULL, 10);
  for (i = 0; i < (world_soil + 1) * world_width; i++)
  {
    fgets(buf, MAX_SLEN, f);
    soil_index[i] = strtol(buf, (char **) NULL, 10);
  }
  psize = 0;
  for (i = 0; i < world_width; i++)
  {
    fgets(buf, MAX_SLEN, f);
    j = strtol(buf, (char **) NULL, 10);
    if (j > -1)
    {
      plant[i] = malloc(sizeof(PLANT));
      if (plant[i] == NULL)
      {
        do_error("load_savefile: failed to allocate memory for plant");
        clear_arrays();
	fclose(f);
        return (-1);
      }
      plant[i]->num_cells = j;
      fgets(buf, MAX_SLEN, f);
      plant[i]->age = strtol(buf, (char **) NULL, 10);
      fgets(buf, MAX_SLEN, f);
      plant[i]->energy_pool = strtol(buf, (char **) NULL, 10);
      fgets(buf, MAX_SLEN, f);
      plant[i]->nutrient_pool = strtol(buf, (char **) NULL, 10);
      if (phyltest_savefreq > 0)
        gn_read_id(&(plant[i]->node), f);
      plant[i]->cellular_energy = 0;
      plant[i]->cellular_nutrient = 0;
      plant[i]->cell = malloc(plant[i]->num_cells * sizeof(PL_CELL));
      if (plant[i]->cell == NULL)
      {
        do_error("load_savefile: failed to allocate memory for plant cells");
        clear_arrays();
	fclose(f);
        return (-1);
      }
      for (j = 0; j < plant[i]->num_cells; j++)
      {
        read_plant_cell(f, &(plant[i]->cell[j]));
        plant[i]->cellular_energy += plant[i]->cell[j].energy;
        plant[i]->cellular_nutrient += plant[i]->cell[j].nutrient;
        world[plant[i]->cell[j].x][plant[i]->cell[j].y].plant_no = i;
        world[plant[i]->cell[j].x][plant[i]->cell[j].y].cell_no = j;
      }
      if (read_genome(f, &(plant[i]->genome), 0xffffffff) < 0)
      {
        clear_arrays();
	fclose(f);
        return (-1);
      }
      psize++;
    }
  }
  for (i = 0; i < world_width; i++)
  {
    for (j = 0; j <= world_soil; j++)
    {
      fgets(buf, MAX_SLEN, f);
      world[i][j].nutrient = strtol(buf, (char **) NULL, 10);
      fgets(buf, MAX_SLEN, f);
      world[i][j].organic_nutrient = strtol(buf, (char **) NULL, 10);
    }
  }
  fclose(f);
  if (expected_psize != psize)
  {
    clear_arrays();
    printf("expected %ld plants, found %ld plants\n", expected_psize, psize);
    do_error("load_savefile: expected and actual population size inconsistent");
    return (-1);
  }
  return (0);
}


int load_savefile(void)
{
  prepare_filenames();
  return (load_named_savefile(save_file_name));
}

