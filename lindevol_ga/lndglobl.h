#include "lndvals.h"
#include "lndtypes.h"

#include "gnlib.h"


#ifdef LND_MAIN

char testmode  = 0;
char quietmode = 0;
char worldmode = 0;
char pixmode   = 0;
char pixel_file_name[MAX_SLEN];

char finish_flag = 0;

char simname[MAX_SLEN];
char savetime_fname[MAX_SLEN];
char par_fname[MAX_SLEN];

unsigned long psize;          /* population size */
double srate;                 /* selection rate */
unsigned long num_select;     /* number of discarded individuals per generation */
double recombination;         /* recombination rate */
unsigned long num_recombine;  /* number of recombinations per generation */
double m_replacement;         /* replacement rate */
double m_insertion;           /* insertion rate */
double m_deletion;            /* deletion rate */
double m_factor;
GSYS_PARAMETERS gsys_parameters;

long random_seed;

long world_width, world_height;
long plant_distance;          /* number of lattice sites between plants */
unsigned long num_days;       /* number of days in plant growth simulation */
unsigned long glen_init;      /* length of randomly initialized genomes */

unsigned long start_generation;
unsigned long num_generations; /* number of generations to simulate */
unsigned long dmt_savefreq;    /* distance matrix is saved every so many generations */
unsigned long bp_ceiling;      /* ceiling for B&P evolutionary activity analysis */
unsigned long ddistr_freq;     /* frequency of saving of distance distribution (0: none) */
unsigned long phhist_freq;     /* frequency of saving of phylogenetic trees (0: none) */

unsigned long generation;      /* the current generation */

double average_fitness;
double average_genome_length;
unsigned long min_genome_length, max_genome_length;
double average_num_used;
unsigned long min_num_used, max_num_used;
unsigned long num_species;
long mainspec_distance;
double distance_entropy = 0.0;
double distance_entropy_rel = 0.0;
long min_mutflag, max_mutflag;
double average_mutflag;
double genetic_diversity, rel_genetic_diversity;

LND_GENOME         *lnd_genome = (LND_GENOME *) NULL;          /* the array of genomes */
unsigned long  *gi = (unsigned long *) NULL;       /* genome index for standard access. */
unsigned long  *gi_tmp = (unsigned long *) NULL;   /* genome index for counting species, may be used for other purposes */
GENOME          old_main_species;

SPECIES_D      *species_d = (SPECIES_D *) NULL;    /* array of species descriptors (for counting species etc) */

LATTICE_SITE  **world = (LATTICE_SITE **) NULL;    /* the array containing the plant cell locations */
PLANT          *plant = (PLANT *) NULL;            /* the coordinates and the state of the plant cells */
unsigned long  *pl_index = (unsigned long *) NULL; /* plant index for randomly shuffling the order of growing */

GN_TREE  gn_tree;

const int x_offset[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
const int y_offset[8] = {-1, -1, -1, 0, 0, 1, 1, 1};

#endif /* LND_MAIN */

extern char testmode;
extern char quietmode;
extern char worldmode;
extern char pixmode;
extern char pixel_file_name[MAX_SLEN];

extern char simname[MAX_SLEN];
extern char savetime_fname[MAX_SLEN];
extern char par_fname[MAX_SLEN];

extern unsigned long psize;          /* population size */
extern double srate;                 /* selection rate */
extern unsigned long num_select;     /* number of discarded individuals per generation */
extern double recombination;         /* recombination rate */
extern unsigned long num_recombine;  /* number of recombinations per generation */
extern double m_replacement;         /* replacement rate */
extern double m_insertion;           /* insertion rate */
extern double m_deletion;            /* deletion rate */
extern double m_factor;
extern GSYS_PARAMETERS gsys_parameters;

extern long random_seed;

extern long world_width, world_height;
extern long plant_distance;          /* number of lattice sites between plants */
extern unsigned long num_days;       /* number of days in growth simulation */
extern unsigned long glen_init;      /* length of randomly initialized genomes */

extern unsigned long start_generation;
extern unsigned long num_generations; /* number of generations to simulate */
extern unsigned long dmt_savefreq;    /* distance matrix is saved every so many generations */
extern unsigned long bp_ceiling;      /* ceiling for B&P evolutionary activity analysis */
extern unsigned long ddistr_freq;             /* frequency of saving of distance distribution (0: none) */
extern unsigned long phhist_freq;             /* frequency of saving of phylogenetic trees (0: none) */

extern unsigned long generation;      /* the current generation */

extern double average_fitness;
extern double average_genome_length;
extern unsigned long min_genome_length, max_genome_length;
extern double average_num_used;
extern unsigned long min_num_used, max_num_used;
extern unsigned long num_species;
extern long mainspec_distance;
extern double distance_entropy, distance_entropy_rel;
extern long min_mutflag, max_mutflag;
extern double average_mutflag;
extern double genetic_diversity, rel_genetic_diversity;

extern LND_GENOME         *lnd_genome;       /* the array of genomes */
extern unsigned long  *gi;           /* genome index for standard access. */
extern unsigned long  *gi_tmp;       /* genome index for counting species, may be used for other purposes */
extern GENOME          main_species; /* the main species */
extern GENOME          old_main_species;

extern SPECIES_D      *species_d;    /* array of species descriptors (for counting species etc) */

extern LATTICE_SITE  **world;        /* the array containing the plant cell locations */
extern PLANT          *plant;        /* the coordinates and the state of the plant cells */
extern unsigned long  *pl_index;     /* plant index for randomly shuffling the order of growing */

extern GN_TREE gn_tree;

extern const int x_offset[8];
extern const int y_offset[8];

