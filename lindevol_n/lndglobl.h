#ifndef H_LNDGLOBL
#define H_LNDGLOBL

#include "lndtypes.h"
#include "gntypes.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef LND_MAIN

char buf[MAX_SLEN];

char par_fname[MAX_SLEN];       /* name of parameter file */
char g_fname[MAX_SLEN] = "\0";  /* name of genome file */

char quietmode = 0;
char worldmode = 0;
char pixmode   = 0;

char finish_flag = 0;
char save_data = 0;

char simname[MAX_SLEN];
char savetime_fname[MAX_SLEN];

long psize_init;               /* initial population size */
long psize;                    /* current population size */
double m_replacement;          /* replacement rate */
double m_insertion;            /* insertion rate */
double m_deletion;             /* deletion rate */
double m_duplication;          /* gene duplication rate */
double m_factor;               /* multiplication factor per mut- / mut+ */

long world_width, world_height;
long world_soil;
long glen_init;                /* length of randomly initialized genomes */
double p_random_death;         /* probability for random death */
double rdeath_f_energy;        /* total energy dependent factor for random death */
double rdeath_f_numcells;      /* total number of cells dependent factor for random death */
double leanover_penalty;       /* probability of dying per unit of cell imbalance */
                               /* and time step */
long seedprod_threshold;       /* min. number of cells required for a plant */
                               /* to be able to reproduce */
long nutrient_init;            /* initial amount of nutrient in the world */
double diffusion_rate;         /* probability of nutrient unit to diffuse per time step */
double decomposition_rate;     /* probability of nutrient unit to move from organic to free state per time step */

GSYS_PARAMETERS gsys_parameters;

long start_generation;         /* the generation at which the simulation should start */
long num_generations;          /* number of generations to simulate */
long dmt_savefreq;             /* distance matrix is saved every so many generations */
long ddistr_savefreq;          /* frequency of saving distance distributions */
long ddistr_samplesize;        /* sample size for distance distribution analysis */
long phyltest_savefreq;        /* genomes are dumped every so many generations */
long phyltest_samplesize;      /* generate random samples of that size for phylogenetic analysis */
long bp_ceiling;               /* ceiling for B&P evolutionary activity analysis */
long soil_savefreq;            /* frequency of saving horizontal and vertical soil profile */

long random_seed = 1;           /* the random seed */

long generation;               /* the current generation */

long num_grownplants;          /* number of plants processed in time step */

long min_num_cells, max_num_cells;
double average_num_cells;
long min_cellular_energy, max_cellular_energy;
double average_cellular_energy;
long min_cellular_nutrient, max_cellular_nutrient;
double average_cellular_nutrient;
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
long min_mutcounter, max_mutcounter;
double average_mutcounter;
long num_unmutated_genomes;
long min_num_genes, max_num_genes;
double average_num_genes;
double genetic_diversity, rel_genetic_diversity;
double min_deathprob, max_deathprob, sum_deathprob;
long min_specificity, max_specificity;
double average_specificity;
long num_bitchecks[NUM_STATEBITS], num_bitsets[NUM_STATEBITS];
long min_energy_pool, max_energy_pool;
double average_energy_pool;
long min_nutrient_pool, max_nutrient_pool;
double average_nutrient_pool;
long num_free_nutrient, num_organic_nutrient, num_biomass_nutrient;
long num_to_epool, num_to_npool, num_from_epool, num_from_npool;
long *soil_profile_h = NULL, *soil_profile_v = NULL;
long *organic_profile_h = NULL, *organic_profile_v = NULL;

char pro_file_name[MAX_SLEN];
char dmt_file_name[MAX_SLEN];
char dst_file_name[MAX_SLEN];
char jf_file_name[MAX_SLEN];
char genome_file_name[MAX_SLEN];
char bpe_file_name[MAX_SLEN];
char save_file_name[MAX_SLEN];
char pixel_file_name[MAX_SLEN];

GN_TREE gntree;

PLANT         **plant = (PLANT **) NULL;        /* array of pointers to plants */
GENOME          old_main_species;               /* the main species of previous time step */

SPECIES_D      *species_d = (SPECIES_D *) NULL; /* array of species descriptors (for counting species etc) */

LATTICE_SITE  **world = (LATTICE_SITE **) NULL; /* the array containing the plant cell locations */
long  *r_index = (long *) NULL;                 /* general purpose index array (for sorting, randomizing etc.) */
long  *pl_index = (long *) NULL;                /* plant index array for randomly shuffling the order of growth processing */
long  *tmp_index = (long *) NULL;               /* this index array is not initialized upon creation
                                                   and may be written to (in order to be used only partly) */
long *sample_index = (long *) NULL;             /* index for preparing random samples of population */
long nutrient_i = 0;                            /* index for soil_index array */
long *soil_index = (long *) NULL;               /* index for randomly choosing soil lattice sites */

const int x_offset[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
const int y_offset[8] = {-1, -1, -1, 0, 0, 1, 1, 1};

#endif /* LND_MAIN */

extern char buf[MAX_SLEN];

extern char quietmode;
extern char worldmode;
extern char pixmode;

extern char finish_flag;
extern char save_data;

extern char simname[MAX_SLEN];
extern char savetime_fname[MAX_SLEN];

extern long psize_init;               /* initial population size */
extern long psize;                    /* current population size */
extern double m_replacement;          /* replacement rate */
extern double m_insertion;            /* insertion rate */
extern double m_deletion;             /* deletion rate */
extern double m_duplication;          /* gene duplication rate */
extern double m_factor;               /* multiplication factor per mut- / mut+ */

extern long world_width, world_height;
extern long world_soil;
extern long glen_init;                /* length of randomly initialized genomes */
extern double p_random_death;         /* probability for random death */
extern double rdeath_f_energy;        /* total energy dependent factor for random death */
extern double rdeath_f_numcells;
extern double leanover_penalty;
extern long seedprod_threshold;       /* min. number of cells required for a plant */
                                      /* to be able to reproduce */
extern long nutrient_init;            /* initial amount of nutrient in the world */
extern double diffusion_rate;         /* probability of nutrient unit to diffuse per time step */
extern double decomposition_rate;     /* probability of nutrient unit to move from organic to free state per time step */

extern GSYS_PARAMETERS gsys_parameters;

extern long random_seed;               /* the random seed */

extern long start_generation;         /* the generation at which the simulation should start */
extern long num_generations;          /* number of generations to simulate */
extern long dmt_savefreq;             /* distance matrix is saved every so many generations */
extern long ddistr_savefreq;          /* distance distributions are saved every so many generations */
extern long ddistr_samplesize;        /* sample size for distance distribution analysis */
extern long phyltest_savefreq;        /* genomes are dumped every so many generations */
extern long phyltest_samplesize;      /* generate random samples of that size for phylogenetic analysis */
extern long bp_ceiling;               /* ceiling for B&P evolutionary activity analysis */
extern long soil_savefreq;            /* frequency of saving horizontal and vertical soil profile */

extern long generation;               /* the current generation */

extern long num_grownplants;          /* number of plants processed in time step */

extern long min_num_cells, max_num_cells;
extern double average_num_cells;
extern long min_cellular_energy, max_cellular_energy;
extern double average_cellular_energy;
extern long min_cellular_nutrient, max_cellular_nutrient;
extern double average_cellular_nutrient;
extern double average_genome_length;
extern long min_genome_length, max_genome_length;
extern double average_num_used;
extern long min_num_used, max_num_used;
extern long num_species;
extern long mainspec_distance;
extern double distance_entropy;
extern double distance_entropy_rel;
extern long num_seeds, num_local_seeds, num_flying_seeds;
extern long num_new_plants;
extern long num_divisions, num_new_cells;
extern long num_attacks, num_deaths;
extern long min_age, max_age;
extern double average_age;
extern long num_mutminus, num_mutplus;
extern long min_mutcounter, max_mutcounter;
extern double average_mutcounter;
extern long min_num_genes, max_num_genes;
extern double average_num_genes;
extern long num_unmutated_genomes;
extern double genetic_diversity, rel_genetic_diversity;
extern double min_deathprob, max_deathprob, sum_deathprob;
extern long min_specificity, max_specificity;
extern double average_specificity;
extern long num_bitchecks[NUM_STATEBITS], num_bitsets[NUM_STATEBITS];
extern long min_energy_pool, max_energy_pool;
extern double average_energy_pool;
extern long min_nutrient_pool, max_nutrient_pool;
extern double average_nutrient_pool;
extern long num_free_nutrient, num_organic_nutrient, num_biomass_nutrient;
extern long num_to_epool, num_to_npool, num_from_epool, num_from_npool;
extern long *soil_profile_h, *soil_profile_v, *organic_profile_h, *organic_profile_v;

extern char pro_file_name[MAX_SLEN];
extern char dmt_file_name[MAX_SLEN];
extern char dst_file_name[MAX_SLEN];
extern char jf_file_name[MAX_SLEN];
extern char genome_file_name[MAX_SLEN];
extern char bpe_file_name[MAX_SLEN];
extern char save_file_name[MAX_SLEN];
extern char pixel_file_name[MAX_SLEN];

extern GN_TREE gntree;

extern PLANT         **plant;             /* array of pointers to plants */
extern GENOME          main_species;      /* the main species */
extern GENOME          old_main_species;

extern SPECIES_D      *species_d;         /* array of species descriptors (for counting species etc) */

extern LATTICE_SITE **world;              /* the array containing the plant cell locations */
extern long  *r_index;                    /* general purpose index array (for sorting, randomizing etc.) */
extern long  *pl_index;                   /* plant index array for randomly shuffling the order of growth processing */
extern long  *tmp_index;                  /* this index array is not initialized upon creation
                                             and may be written to (in order to be used only partly) */
extern long *sample_index;                /* index for preparing random samples of population */
extern long nutrient_i;                   /* index for soil_index array */
extern long *soil_index;                  /* index for randomly choosing soil lattice sites */

extern const int x_offset[8], y_offset[8];

#ifdef __cplusplus
}
#endif

#endif /* H_LNDGLOBL */

