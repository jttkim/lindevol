#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "jklib.h"

#include "lndvals.h"
#include "lndtypes.h"
#include "lndglobl.h"
#include "lndlibin.h"
#include "lndlib.h"
#include "lnd6.h"
#include "lnderror.h"


FILE *pro_file = NULL;


int open_pro_file(const char *mode)
{
  time_t simtime;

  pro_file = fopen(pro_file_name, mode);
  if (pro_file == NULL)
  {
    perror("open_pro_file failed");
    return (-1);
  }
  if (*mode == 'w')
  {
    fprintf(pro_file, "%s\n", SIMPRGNAME);
    simtime = time(NULL);
    fprintf(pro_file, "%s", ctime(&simtime));
    fprintf(pro_file, "%s", ctime(&simtime));
    fprintf(pro_file, "# run name                               : %s\n", simname);
    fprintf(pro_file, "# initial genome length               (L): %ld\n", glen_init);
    fprintf(pro_file, "# initial population size             (P): %ld\n", psize_init);
    fprintf(pro_file, "# mutation rate for replacement       (M): %1.12g\n", m_replacement);
    fprintf(pro_file, "# mutation rate for insertion        (Mi): %1.12g\n", m_insertion);
    fprintf(pro_file, "# mutation rate for deletion         (Md): %1.12g\n", m_deletion);
    fprintf(pro_file, "# mutation rate for gene duplication (Mg): %1.12g\n", m_duplication);
    fprintf(pro_file, "# mutation factor                    (Mf): %1.12g\n", m_factor);
    fprintf(pro_file, "# world width                            : %ld\n", world_width);
    fprintf(pro_file, "# world height                           : %ld\n", world_height);
    fprintf(pro_file, "# probability for random death           : %1.12g\n", p_random_death);
    fprintf(pro_file, "# energy dependent rdeath factor         : %1.12g\n", rdeath_f_energy);
    fprintf(pro_file, "# number of cells dependent rdeath factor: %1.12g\n", rdeath_f_numcells);
    fprintf(pro_file, "# leanover penalty                       : %1.12g\n", leanover_penalty);
    fprintf(pro_file, "# seed production threshold              : %ld\n", seedprod_threshold);
    fprintf(pro_file, "# number of divide actions               : %ld\n", gsys_parameters.num_divide);
    fprintf(pro_file, "# number of flying seed actions          : %ld\n", gsys_parameters.num_flyingseed);
    fprintf(pro_file, "# number of local seed actions           : %ld\n", gsys_parameters.num_localseed);
    fprintf(pro_file, "# number of mut- actions                 : %ld\n", gsys_parameters.num_mutminus);
    fprintf(pro_file, "# number of mut+ actions                 : %ld\n", gsys_parameters.num_mutplus);
    fprintf(pro_file, "# number of set statebit actions         : %ld\n", gsys_parameters.num_statebit);
    fprintf(pro_file, "# random seed                            : %ld\n", random_seed);

    fprintf(pro_file, "@\n");
    fprintf(pro_file, "SIMULATOR=%s\n", SIMPRGNAME);
    fprintf(pro_file, "STARTTIME=%s", ctime(&simtime));
    fprintf(pro_file, "SIMNAME=%s\n", simname);
    fprintf(pro_file, "PSIZE_INIT=%ld\n", psize_init);
    fprintf(pro_file, "MUTRATES=%f %f %f %f, f=%f\n", m_replacement, m_insertion, m_deletion, m_duplication, m_factor);
    fprintf(pro_file, "MR=%f\n", m_replacement);
    fprintf(pro_file, "MI=%f\n", m_insertion);
    fprintf(pro_file, "MD=%f\n", m_deletion);
    fprintf(pro_file, "MG=%f\n", m_duplication);
    fprintf(pro_file, "MF=%f\n", m_factor);
    fprintf(pro_file, "GLEN_INIT=%ld\n", glen_init);
    fprintf(pro_file, "DEATH_P=%f\n", p_random_death);
    fprintf(pro_file, "DEATH_FC=%f\n", rdeath_f_numcells);
    fprintf(pro_file, "DEATH_FE=%f\n", rdeath_f_energy);
    fprintf(pro_file, "DEATH_LEANOVER=%f\n", leanover_penalty);
    fprintf(pro_file, "SEED_THRESHOLD=%ld\n", seedprod_threshold);
    fprintf(pro_file, "NUTRIENT_INIT=%ld\n", nutrient_init);
    fprintf(pro_file, "DIFFUSION_RATE=%f\n", diffusion_rate);
    fprintf(pro_file, "DECOMPOSITION_RATE=%f\n", decomposition_rate);
    fprintf(pro_file, "RANDOM_SEED=%ld\n", random_seed);
    fprintf(pro_file, "GSYS_CONF=d%ld f%ld l%ld -%ld +%ld et%ld nt%ld ef%ld nf%ld s%ld\n",
            gsys_parameters.num_divide,
            gsys_parameters.num_flyingseed, gsys_parameters.num_localseed,
            gsys_parameters.num_mutminus, gsys_parameters.num_mutplus,
	    gsys_parameters.num_to_epool, gsys_parameters.num_to_npool,
	    gsys_parameters.num_from_epool, gsys_parameters.num_from_npool,
            gsys_parameters.num_statebit);
    fprintf(pro_file, "WORLD_WIDTH=%ld\n", world_width);
    fprintf(pro_file, "WORLD_HEIGHT=%ld\n", world_height);
    fprintf(pro_file, "WORLD_SOIL=%ld\n", world_soil);
    fprintf(pro_file, "PHYLTEST_SAVEFREQ=%ld\n", phyltest_savefreq);
    fprintf(pro_file, "PHYLTEST_SAMPLESIZE=%ld\n", phyltest_samplesize);
    fprintf(pro_file, "DDISTR_SAVEFREQ=%ld\n", ddistr_savefreq);
    fprintf(pro_file, "DDISTR_SAMPLESIZE=%ld\n", ddistr_samplesize);
    fprintf(pro_file, "@\n");

    fprintf(pro_file, "n\n");
    fprintf(pro_file, "generation #\n");                             /* 000 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "population size\n");                          /* 001 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. number of cells\n");                     /* 002 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. number of cells\n");                     /* 003 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average number of cells\n");                  /* 004 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. cellular energy\n");                     /* 005 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. cellular energy\n");                     /* 006 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average cellular energy\n");                  /* 007 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of different genomes\n");              /* 008 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "copies of most frequent genome\n");           /* 009 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "copies of second most frequent genome\n");    /* 010 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "copies of third most frequent genome\n");     /* 011 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "edit distance to previous main species\n");   /* 012 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. genome length in chars\n");              /* 013 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. genome length in chars\n");              /* 014 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average genome length in chars\n");           /* 015 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. number of genes\n");                     /* 016 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. number of genes\n");                     /* 017 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average number of genes\n");                  /* 018 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. number of used genes\n");                /* 019 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. number of used genes\n");                /* 020 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average number of used genes\n");             /* 021 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. specificity\n");                         /* 022 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. specificity\n");                         /* 023 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average specificity\n");                      /* 024 */
    if (ddistr_savefreq > 0)
    {
      fprintf(pro_file, "n i%1ld\n", ddistr_savefreq);
      fprintf(pro_file, "distance distribution entropy\n");            /* 022 */
      fprintf(pro_file, "n i%1ld\n", ddistr_savefreq);
      fprintf(pro_file, "relative distance distribution entropy\n");   /* 023 */
    }
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of seeds\n");                          /* 024 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of local seeds\n");                    /* 025 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of flying seeds\n");                   /* 026 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of new plants\n");                     /* 027 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of cell divisions\n");                 /* 028 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of new cells\n");                      /* 029 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of attacks\n");                        /* 030 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of deaths\n");                         /* 031 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. age\n");                                 /* 032 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. age\n");                                 /* 033 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average age\n");                              /* 034 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of mut- actions\n");                   /* 035 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of mut+ actions\n");                   /* 036 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. mutation counter value\n");              /* 037 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. mutation counter value\n");              /* 038 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average mutation counter value\n");           /* 039 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "genetic diversity\n");                        /* 040 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "relative genetic diversity\n");               /* 041 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "minimal death probability\n");                /* 042 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "maximal death probability\n");                /* 043 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average death probability\n");                /* 044 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. cellular nutrient\n");                   /* 045 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. cellular nutrient\n");                   /* 046 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average cellular nutrient\n");                /* 047 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. energy pool\n");                         /* 048 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. energy pool\n");                         /* 049 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average energy pool\n");                      /* 050 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. nutrient pool\n");                       /* 051 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. nutrient pool\n");                       /* 052 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average nutrient pool\n");                    /* 053 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "free nutrient\n");                            /* 054 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "organic nutrient\n");                         /* 055 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "biomass nutrient\n");                         /* 056 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of to_epool\n");                       /* 057 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of to_npool\n");                       /* 058 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of from_epool\n");                     /* 059 */
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of from_npool\n");                     /* 060 */
    fprintf(pro_file, "d\n");
    fprintf(pro_file, "bitcheck distribution\n");                    /* 061 */
    fprintf(pro_file, "d\n");
    fprintf(pro_file, "bitset distribution\n");                      /* 062 */
    if (soil_savefreq > 0)
    {
      fprintf(pro_file, "d i%1ld\n", soil_savefreq);
      fprintf(pro_file, "horizontal soil profile\n");                /* 063 */
      fprintf(pro_file, "d i%1ld\n", soil_savefreq);
      fprintf(pro_file, "vertical soil profile\n");                  /* 064 */
      fprintf(pro_file, "d i%1ld\n", soil_savefreq);
      fprintf(pro_file, "horizontal organic profile\n");             /* 065 */
      fprintf(pro_file, "d i%1ld\n", soil_savefreq);
      fprintf(pro_file, "vertical organic profile\n");               /* 066 */
    }
    if (ddistr_savefreq > 0)
    {
      fprintf(pro_file, "d f%s i%1ld\n", dst_file_name, ddistr_savefreq);
      fprintf(pro_file, "distance distribution\n");
    }
    if (bp_ceiling > 0)
    {
      fprintf(pro_file, "d f%s\n", bpe_file_name);
      fprintf(pro_file, "B&P evolutionary activity\n");
    }
    if (phyltest_savefreq > 0)
    {
      fprintf(pro_file, "s %s i%1ld\n", jf_file_name, phyltest_savefreq);
      fprintf(pro_file, "phylogenetic trees\n");
    }
    fprintf(pro_file, "@\n");
  }
  return (0);
}


void write_pro(void)
{
  long i, l;
  long *c_hgram, ch_length;

  fprintf(pro_file, "%ld\n", generation);
  fprintf(pro_file, "%ld\n", psize);
  fprintf(pro_file, "%ld\n", min_num_cells);
  fprintf(pro_file, "%ld\n", max_num_cells);
  fprintf(pro_file, "%1.12g\n", average_num_cells);
  fprintf(pro_file, "%ld\n", min_cellular_energy);
  fprintf(pro_file, "%ld\n", max_cellular_energy);
  fprintf(pro_file, "%1.12g\n", average_cellular_energy);
  fprintf(pro_file, "%ld\n", num_species);
  fprintf(pro_file, "%ld\n", species_d[0].num);
  fprintf(pro_file, "%ld\n", species_d[1].num);
  fprintf(pro_file, "%ld\n", species_d[2].num);
  fprintf(pro_file, "%ld\n", mainspec_distance);
  fprintf(pro_file, "%ld\n", min_genome_length);
  fprintf(pro_file, "%ld\n", max_genome_length);
  fprintf(pro_file, "%1.12g\n", average_genome_length);
  fprintf(pro_file, "%ld\n", min_num_genes);
  fprintf(pro_file, "%ld\n", max_num_genes);
  fprintf(pro_file, "%1.12g\n", average_num_genes);
  fprintf(pro_file, "%ld\n", min_num_used);
  fprintf(pro_file, "%ld\n", max_num_used);
  fprintf(pro_file, "%1.12g\n", average_num_used);
  fprintf(pro_file, "%ld\n", min_specificity);
  fprintf(pro_file, "%ld\n", max_specificity);
  fprintf(pro_file, "%1.12g\n", average_specificity);
  if ((ddistr_savefreq > 0) && ((generation % ddistr_savefreq) == 0))
  {
    fprintf(pro_file, "%1.12g\n", distance_entropy);
    fprintf(pro_file, "%1.12g\n", distance_entropy_rel);
  }
  fprintf(pro_file, "%ld\n", num_seeds);
  fprintf(pro_file, "%ld\n", num_local_seeds);
  fprintf(pro_file, "%ld\n", num_flying_seeds);
  fprintf(pro_file, "%ld\n", num_new_plants);
  fprintf(pro_file, "%ld\n", num_divisions);
  fprintf(pro_file, "%ld\n", num_new_cells);
  fprintf(pro_file, "%ld\n", num_attacks);
  fprintf(pro_file, "%ld\n", num_deaths);
  fprintf(pro_file, "%ld\n", min_age);
  fprintf(pro_file, "%ld\n", max_age);
  fprintf(pro_file, "%1.12g\n", average_age);
  fprintf(pro_file, "%ld\n", num_mutminus);
  fprintf(pro_file, "%ld\n", num_mutplus);
  fprintf(pro_file, "%ld\n", min_mutcounter);
  fprintf(pro_file, "%ld\n", max_mutcounter);
  fprintf(pro_file, "%1.12g\n", average_mutcounter);
  fprintf(pro_file, "%1.12g\n", genetic_diversity);
  fprintf(pro_file, "%1.12g\n", rel_genetic_diversity);
  fprintf(pro_file, "%1.12g\n", min_deathprob);
  fprintf(pro_file, "%1.12g\n", max_deathprob);
  if (num_grownplants)
    fprintf(pro_file, "%1.12g\n", sum_deathprob / num_grownplants);
  else
    fprintf(pro_file, "0.0\n");
  fprintf(pro_file, "%ld\n", min_cellular_nutrient);
  fprintf(pro_file, "%ld\n", max_cellular_nutrient);
  fprintf(pro_file, "%1.12g\n", average_cellular_nutrient);
  fprintf(pro_file, "%ld\n", min_energy_pool);
  fprintf(pro_file, "%ld\n", max_energy_pool);
  fprintf(pro_file, "%1.12g\n", average_energy_pool);
  fprintf(pro_file, "%ld\n", min_nutrient_pool);
  fprintf(pro_file, "%ld\n", max_nutrient_pool);
  fprintf(pro_file, "%1.12g\n", average_nutrient_pool);
  fprintf(pro_file, "%ld\n", num_free_nutrient);
  fprintf(pro_file, "%ld\n", num_organic_nutrient);
  fprintf(pro_file, "%ld\n", num_biomass_nutrient);
  fprintf(pro_file, "%ld\n", num_to_epool);
  fprintf(pro_file, "%ld\n", num_to_npool);
  fprintf(pro_file, "%ld\n", num_from_epool);
  fprintf(pro_file, "%ld\n", num_from_npool);
  ch_length = world_width;
  ch_length = (ch_length > (world_soil + 1)) ? ch_length : (world_soil + 1);
  ch_length = (ch_length > NUM_STATEBITS) ? ch_length : NUM_STATEBITS;
  /* fprintf(stderr, "ch_length = %ld\n", ch_length); */
  c_hgram = (long *) malloc(ch_length * sizeof(long));
  if (c_hgram == NULL)
  {
    do_error("write_pro: failed to allocate c_hgram, saving histograms uncompressed");
    fprintf(pro_file, "%ld\n", (long) NUM_STATEBITS);
    for (i = 0; i < NUM_STATEBITS; i++)
      fprintf(pro_file, "%ld ", num_bitchecks[i]);
    fprintf(pro_file, "\n");
    fprintf(pro_file, "%ld\n", (long) NUM_STATEBITS);
    for (i = 0; i < NUM_STATEBITS; i++)
      fprintf(pro_file, "%ld ", num_bitsets[i]);
    fprintf(pro_file, "\n");
    if ((soil_savefreq > 0) && ((generation % soil_savefreq) == 0))
    {
      fprintf(pro_file, "%ld\n", world_width);
      for (i = 0; i < world_width; i++)
	fprintf(pro_file, "%ld ", soil_profile_h[i]);
      fprintf(pro_file, "\n");
      fprintf(pro_file, "%ld\n", world_soil + 1);
      for (i = 0; i < world_soil + 1; i++)
	fprintf(pro_file, "%ld ", soil_profile_v[i]);
      fprintf(pro_file, "\n");
      fprintf(pro_file, "%ld\n", world_width);
      for (i = 0; i < world_width; i++)
	fprintf(pro_file, "%ld ", organic_profile_h[i]);
      fprintf(pro_file, "\n");
      fprintf(pro_file, "%ld\n", world_soil + 1);
      for (i = 0; i < world_soil + 1; i++)
	fprintf(pro_file, "%ld ", organic_profile_v[i]);
      fprintf(pro_file, "\n");
    }
  }
  else
  {
    l = compress_histogram_long(NUM_STATEBITS, num_bitchecks, c_hgram);
    /* fprintf(stderr, "compressed length of num_bitchecks: %ld\n", l); */
    fprintf(pro_file, "%ld\n", l);
    for (i = 0; i < l; i++)
      fprintf(pro_file, "%ld ", c_hgram[i]);
    fprintf(pro_file, "\n");
    l = compress_histogram_long(NUM_STATEBITS, num_bitsets, c_hgram);
    /* fprintf(stderr, "compressed length of num_bitsets: %ld\n", l); */
    fprintf(pro_file, "%ld\n", l);
    for (i = 0; i < l; i++)
      fprintf(pro_file, "%ld ", c_hgram[i]);
    fprintf(pro_file, "\n");
    if ((soil_savefreq > 0) && ((generation % soil_savefreq) == 0))
    {
      l = compress_histogram_long(world_width, soil_profile_h, c_hgram);
      fprintf(pro_file, "%ld\n", l);
      for (i = 0; i < l; i++)
	fprintf(pro_file, "%ld ", c_hgram[i]);
      fprintf(pro_file, "\n");
      l = compress_histogram_long(world_soil + 1, soil_profile_v, c_hgram);
      fprintf(pro_file, "%ld\n", l);
      for (i = 0; i < l; i++)
	fprintf(pro_file, "%ld ", c_hgram[i]);
      fprintf(pro_file, "\n");
      l = compress_histogram_long(world_width, organic_profile_h, c_hgram);
      fprintf(pro_file, "%ld\n", l);
      for (i = 0; i < l; i++)
	fprintf(pro_file, "%ld ", c_hgram[i]);
      fprintf(pro_file, "\n");
      l = compress_histogram_long(world_soil + 1, organic_profile_v, c_hgram);
      fprintf(pro_file, "%ld\n", l);
      for (i = 0; i < l; i++)
	fprintf(pro_file, "%ld ", c_hgram[i]);
      fprintf(pro_file, "\n");
    }
    free(c_hgram);
  }
  fflush(pro_file);
}

