#ifndef H_LNDLIB
#define H_LNDLIB

#include <stdio.h>  /* so FILE is defined */

#include "lndtypes.h"

#define free0(p) free(p); (p) = NULL

#ifdef __cplusplus
extern "C" {
#endif

extern int create_arrays(void);
extern void clear_arrays(void);
extern void count_species(void);
extern void prepare_filenames(void);
extern int open_pro_file(const char *mode);
extern int open_dmt_file(const char *mode);
extern int open_jf_file(const char *mode);
extern int open_genome_file(const char *mode);
extern int open_dst_file(const char *mode);
extern int open_bpe_file(const char *mode);
extern int open_data_files(const char *mode);
extern void close_data_files(void);
extern void write_pro(void);
extern void write_jf(long num_samples, const long *sample_index);
extern void write_dmt(void);
extern void write_dst(char write_it, long sample_size, long *sample_index);
extern void write_bpe(void);
extern void write_genomes(long num_samples, const long *sample_index);

extern long prepare_sample(long sample_size, long *sample_index);

extern unsigned long *xxx_initstate(unsigned int seed, unsigned long *arg_state);
extern unsigned long *xxx_setstate(unsigned long *arg_state);
extern unsigned long lnd_random(unsigned long range);
extern double lnd_rnd(void);
extern int write_rndgenerator_state(FILE *f);
extern int read_rndgenerator_state(FILE *f);
extern void seed_lnd_random(int seed);
extern void random_shuffle(long num, long *a);

extern long next_mutpos(double m, long pos);

extern double death_probability(long plant_no, double p, double f_energy, double f_numcells, double leanover_penalty);

extern int plant_boundingbox(const PLANT *plant, long *xmin, long *xmax, long *ymax);

extern int open_savetime_file(const char *fname);
extern void close_savetime_file(void);
extern long savetime_next(long generation);

extern int open_pixfile(const char *mode);
extern int write_pixfile(void);

extern long num_energyrich_cells(long plant_no);
extern void add_nutrient_random(void);
extern long random_free_position(void);
extern long next_free_position(long pos);
extern long create_plant(long pos, long generation, const GENOME *parent);
extern void remove_plant(long plant_no);
extern int create_flying_seed(long plant_no, long cell_no);
extern int create_local_seed(long plant_no, long cell_no);
extern void sunshine(void);
extern void nutrient_diffusion(void);

#ifdef __cplusplus
}
#endif

#endif /* H_LNDLIB */

