#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "jklib.h"
#include "genomelib.h"
#include "gfaport.h"
#include "lnderror.h"
#include "lndtypes.h"
#include "lndvals.h"
#include "lndrandm.h"
#include "creatarr.h"
#include "pixfile.h"
#include "lndio.h"

#include "gnlib.h"

#include "lndglobl.h"


static FILE *pro_file; /* the protocol file */
static FILE *dmt_file; /* the distance matrices file */
static FILE *dst_file; /* the distance distribution file */
static FILE *stb_file; /* the phylogenetic tree file */
static FILE *bpe_file; /* the Bedau and Packard evolutionary activity waves file */
static FILE *genome_file; /* the file for saving genomes for phylogenetic reconstruction studies */
static FILE *savetime_file = NULL;

char pro_file_name[MAX_SLEN];
char dmt_file_name[MAX_SLEN];
char dst_file_name[MAX_SLEN];
char stb_file_name[MAX_SLEN];
char genome_file_name[MAX_SLEN];
char bpe_file_name[MAX_SLEN];
char save_file_name[MAX_SLEN];
char savetime_file_name[MAX_SLEN] = "";


void prepare_filenames(void)
{
  sprintf(pro_file_name, "%s%spro", simname, ".");
  sprintf(dmt_file_name, "%s%sdmt", simname, ".");
  sprintf(stb_file_name, "%s%sstb", simname, ".");
  sprintf(genome_file_name, "%s%sgen", simname, ".");
  sprintf(dst_file_name, "%s%sdst", simname, ".");
  sprintf(bpe_file_name, "%s%sbpe", simname, ".");
  sprintf(save_file_name, "%s.dat", simname);
  sprintf(pixel_file_name, "%s.pix", simname);
}


int calculate_activity_codes(GSYS_PARAMETERS *gp)
{
  long sum = gp->num_divide + gp->num_statebit + gp->num_broadcastbit
             + gp->num_flyingseed + gp->num_localseed
             + gp->num_mutminus + gp->num_mutplus;

  if (sum != 256L)
  {
    fprintf(stderr, "sum of numbers of code bytes does not amount to 256\n");
    gp->num_divide = gp->num_divide * 256 / sum;
    gp->num_statebit = gp->num_statebit * 256 / sum;
    gp->num_broadcastbit = gp->num_broadcastbit * 256 / sum;
    gp->num_flyingseed = gp->num_flyingseed * 256 / sum;
    gp->num_localseed = gp->num_localseed * 256 / sum;
    gp->num_mutminus = gp->num_mutminus * 256 / sum;
    gp->num_mutplus = 256 - gp->num_divide - gp->num_statebit - gp->num_broadcastbit
                      -gp->num_flyingseed - gp->num_localseed - gp->num_mutminus;
  }
  gp->divide_code = 0;
  gp->statebit_code = gp->divide_code + gp->num_divide;
  gp->broadcastbit_code = gp->statebit_code + gp->num_statebit;
  gp->flyingseed_code = gp->broadcastbit_code + gp->num_broadcastbit;
  gp->localseed_code = gp->flyingseed_code + gp->num_flyingseed;
  gp->mutminus_code = gp->localseed_code + gp->num_localseed;
  gp->mutplus_code = gp->mutminus_code + gp->num_mutminus;
  return (0);
}


/*
 * return the gene activity corresponding to the output passed
 * in argument.
 */

long gene_activity(unsigned char gene_output)
{
  if (gsys_parameters.num_divide && (gsys_parameters.divide_code <= gene_output) && (gene_output < gsys_parameters.statebit_code))
    return (LND_DIVIDE);
  if (gsys_parameters.num_statebit && (gsys_parameters.statebit_code <= gene_output) && (gene_output < gsys_parameters.broadcastbit_code))
    return (LND_STATEBIT);
  if (gsys_parameters.num_broadcastbit && (gsys_parameters.broadcastbit_code <= gene_output) && (gene_output < gsys_parameters.flyingseed_code))
    return (LND_BROADCASTBIT);
  if (gsys_parameters.num_flyingseed && (gsys_parameters.flyingseed_code <= gene_output) && (gene_output < gsys_parameters.localseed_code))
    return (LND_FLYINGSEED);
  if (gsys_parameters.num_localseed && (gsys_parameters.localseed_code <= gene_output) && (gene_output < gsys_parameters.mutminus_code))
    return (LND_LOCALSEED);
  if (gsys_parameters.num_mutminus && (gsys_parameters.mutminus_code <= gene_output) && (gene_output < gsys_parameters.mutplus_code))
    return (LND_MUTMINUS);
  /* if (gsys_parameters.mutplus_code <= gene_output) */
  return (LND_MUTPLUS);
}


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
    fprintf(pro_file, "%s protocol file\n", SIMPRGNAME);
    simtime = time(NULL);
    fprintf(pro_file, "%s", ctime(&simtime));
    fprintf(pro_file, "%s", ctime(&simtime));
    fprintf(pro_file, "Spielname                              : %s\n", simname);
    fprintf(pro_file, "Stringlnge                         (L): variabel\n");
    fprintf(pro_file, "Populationsgre                    (P): %lu\n", psize);
    fprintf(pro_file, "Selektionsrate                      (S): %1.12g\n", srate);
    fprintf(pro_file, "Anteil rekombinativer Vermehrung    (G): %1.12g\n", recombination);
    fprintf(pro_file, "Anteil vegetativer Vermehrung       (V): %1.12g\n", 1.0 - recombination);
    fprintf(pro_file, "Mutationsrate                       (M): %1.12g\n", m_replacement);
    fprintf(pro_file, "Insertionsrate                     (Mi): %1.12g\n", m_insertion);
    fprintf(pro_file, "Deletionsrate                      (Md): %1.12g\n", m_deletion);
    fprintf(pro_file, "Faktor fuer Mutationsratenmodifikation : %1.12g\n", m_factor);
    fprintf(pro_file, "Random seed                            : %ld\n", random_seed);
    fprintf(pro_file, "Breite der Welt                        : %lu\n", world_width);
    fprintf(pro_file, "Hhe der Welt                          : %lu\n", world_height);
    fprintf(pro_file, "Maximale Anzahl Zellen pro Pflanze     : %lu\n", MAX_NUM_CELLS);
    fprintf(pro_file, "Anzahl Tage pro Wachstumsperiode       : %lu\n", num_days);
    fprintf(pro_file, "Anzahl Gene am Anfang                  : %lu\n", glen_init);
    fprintf(pro_file, "@\n");
    fprintf(pro_file, "SIMULATOR=%s\n", SIMPRGNAME);
    fprintf(pro_file, "SIMNAME=%s\n", simname);
    fprintf(pro_file, "STARTTIME=%s", ctime(&simtime));
    fprintf(pro_file, "GLEN_INIT=%lu\n", glen_init);
    fprintf(pro_file, "POPSIZE=%lu\n", psize);
    fprintf(pro_file, "MUTRATES=%f %f %f, f=%f\n", m_replacement, m_insertion, m_deletion, m_factor);
    fprintf(pro_file, "MD=%f\n", m_replacement);
    fprintf(pro_file, "MI=%f\n", m_insertion);
    fprintf(pro_file, "MD=%f\n", m_deletion);
    fprintf(pro_file, "MF=%f\n", m_factor);
    fprintf(pro_file, "RANDOM_SEED=%ld\n", random_seed);
    fprintf(pro_file, "SRATE=%f\n", srate);
    fprintf(pro_file, "WORLD_WIDTH=%ld\n", world_width);
    fprintf(pro_file, "WORLD_HEIGHT=%ld\n", world_height);
    fprintf(pro_file, "PLANT_DISTANCE=%ld\n", plant_distance);
    fprintf(pro_file, "NUM_DAYS=%lu\n", num_days);
    fprintf(pro_file, "GSYS_CONF=d%ld f%ld l%ld -%ld +%ld\n",
	    gsys_parameters.num_divide,
	    gsys_parameters.num_flyingseed, gsys_parameters.num_localseed,
	    gsys_parameters.num_mutminus, gsys_parameters.num_mutplus);
    fprintf(pro_file, "MUT_FACTOR=%f\n", m_factor);
    fprintf(pro_file, "NUM_DAYS=%lu\n", num_days);
    fprintf(pro_file, "@\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "generation #\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. Fitness\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. Fitness\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average fitness\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "number of different genomes\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "copies of most frequent genome\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "copies of second most frequent genome\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "copies of third most frequent genome\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "edit distance to previous main species\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. genome length in genes\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. genome length in genes\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average genome length in genes\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. number of used genes\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. number of used genes\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average number of used genes\n");
    if (ddistr_freq > 0)
    {
      fprintf(pro_file, "n i%1lu\n", ddistr_freq);
      fprintf(pro_file, "distance distribution entropy\n");
      fprintf(pro_file, "n i%1lu\n", ddistr_freq);
      fprintf(pro_file, "relative distance distribution entropy\n");
    }
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "min. mutation modification exponent\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "max. mutation modification exponent\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "average mutation modification exponent\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "genetic diversity\n");
    fprintf(pro_file, "n\n");
    fprintf(pro_file, "relative genetic diversity\n");
    if (ddistr_freq > 0)
    {
      fprintf(pro_file, "d f%s i%1lu\n", dst_file_name, ddistr_freq);
      fprintf(pro_file, "distance distribution\n");
    }
    if (bp_ceiling > 0)
    {
      fprintf(pro_file, "e f%s i%1lu\n", bpe_file_name, (unsigned long) 1);
      fprintf(pro_file, "B&P evolutionary activity\n");
    }
    if (phhist_freq > 0)
    {
      fprintf(pro_file, "s f%s i%1lu\n", stb_file_name, phhist_freq);
      fprintf(pro_file, "phylogenetic record\n");
    }
    fprintf(pro_file, "@\n");
  }
  return (0);
}


int open_dmt_file(const char *mode)
{
  char m[32];

  sprintf(m, "%sb", mode);
  dmt_file = fopen(dmt_file_name, m);
  if (dmt_file == NULL)
  {
    perror("open_dmt_file failed");
    return (-1);
  }
  if (*mode == 'w')
  {
    fprintf(dmt_file, "INT16\r\n");
    fprintf(dmt_file, "500\r\n");
  }
  return (0);
}


int open_stb_file(const char *mode)
{
  stb_file = fopen(stb_file_name, mode);
  if (stb_file == NULL)
  {
    perror("open_stb_file failed");
    return (-1);
  }
  return (0);
}


int open_genome_file(const char *mode)
{
  genome_file = fopen(genome_file_name, mode);
  if (genome_file == NULL)
  {
    perror("open_genome_file failed");
    return (-1);
  }
  return (0);
}


int open_dst_file(const char *mode)
{
  char m[32];

  sprintf(m, "%sb", mode);
  dst_file = fopen(dst_file_name, m);
  if (dst_file == NULL)
  {
    perror("open_dst_file failed");
    return (-1);
  }
  if (*mode == 'w')
  {
    fprintf(dst_file, "ProPlot-1.0 INT16\r\n");
    fprintf(dst_file, "i%1lu\r\n", ddistr_freq);
  }
  return (0);
}


int open_bpe_file(const char *mode)
{
  char m[32];

  sprintf(m, "%sb", mode);
  bpe_file = fopen(bpe_file_name, m);
  if (bpe_file == NULL)
  {
    perror("open_bpe_file failed");
    return (-1);
  }
  if (*mode == 'w')
  {
    fprintf(bpe_file, "ProPlot-1.0 INT16\r\n");
    fprintf(bpe_file, "i%1lu\r\n", 1UL);
  }
  return (0);
}


void close_data_files(void)
{
  if (pro_file != NULL)
    fclose(pro_file);
  if (dmt_file != NULL)
    fclose(dmt_file);
  if (stb_file != NULL)
    fclose(stb_file);
  if (genome_file != NULL)
    fclose(genome_file);
  if (dst_file != NULL)
    fclose(dst_file);
  if (bpe_file != NULL)
    fclose(bpe_file);
}


void write_pro(void)
{
  fprintf(pro_file, "%lu\n", generation);
  fprintf(pro_file, "%lu\n", lnd_genome[gi[psize - 1]].fitness);
  fprintf(pro_file, "%lu\n", lnd_genome[gi[0]].fitness);
  fprintf(pro_file, "%1.12g\n", average_fitness);
  fprintf(pro_file, "%lu\n", num_species);
  fprintf(pro_file, "%lu\n", species_d[0].num);
  fprintf(pro_file, "%lu\n", species_d[1].num);
  fprintf(pro_file, "%lu\n", species_d[2].num);
  fprintf(pro_file, "%ld\n", mainspec_distance);
  fprintf(pro_file, "%lu\n", min_genome_length);
  fprintf(pro_file, "%lu\n", max_genome_length);
  fprintf(pro_file, "%1.12g\n", average_genome_length);
  fprintf(pro_file, "%lu\n", min_num_used);
  fprintf(pro_file, "%lu\n", max_num_used);
  fprintf(pro_file, "%1.12g\n", average_num_used);
  if ((ddistr_freq > 0) && ((generation % ddistr_freq) == 0))
  {
    fprintf(pro_file, "%1.12g\n", distance_entropy);
    fprintf(pro_file, "%1.12g\n", distance_entropy_rel);
  }
  fprintf(pro_file, "%ld\n", min_mutflag);
  fprintf(pro_file, "%ld\n", max_mutflag);
  fprintf(pro_file, "%1.12g\n", average_mutflag);
  fprintf(pro_file, "%1.12g\n", genetic_diversity);
  fprintf(pro_file, "%1.12g\n", rel_genetic_diversity);
}


/*
PROCEDURE stammbaum
  ' *************************************************************************
  ' * Procedure for saving an exact phylogenetic tree containing a node for *
  ' * each genome in each generation.                                       *
  ' *************************************************************************
*/

void write_stb(void)
{
  int ret_code;

  fprintf(stb_file, "g %lu\n", generation);
  if ((ret_code = gn_print_jftrees(&gn_tree, generation, stb_file)) < 0)
  {
    fprintf(stderr, "*** error %d upon saving tree(s) at generation %lu\n", ret_code, generation);
  }
}


void write_genomes(void)
{
  unsigned long i;

  fprintf(genome_file, "g %lu\n", generation);
  fprintf(genome_file, "%lu\n", psize);
  for (i = 0; i < psize; i++)
  {
    gn_save_id(&(lnd_genome[i].node_id), genome_file);
    write_genome(genome_file, &(lnd_genome[i].genome), 0);
  }
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
  unsigned long i, j;
  short d;

  fprintf(dmt_file, "%lu\r\n", psize);
  fprintf(dmt_file, "%lu\r\n", generation);
  for (i = 0; i < psize - 1; i++)
  {
    for (j = i + 1; j < psize; j++)
    {
      d = (short) edit_distance(lnd_genome[i].genome.length, (char *) lnd_genome[i].genome.g, lnd_genome[j].genome.length, (char *) lnd_genome[j].genome.g);
      fwrite_int16array(&d, 1, dmt_file);
    }
  }
}


/*
PROCEDURE distance_distribution
  ' *******************************************************************
  ' * Saves the distribution of distance values in the matrix of edit *
  ' * distances between all genomes in the population.                *
  ' *******************************************************************
*/

void write_dst(void)
{
  unsigned long i, j;
  unsigned long compress_len;
  long d;
  long max_distance = 0;
  long max_len;
  short dd_rel[100];

  short *dd, *compress;

  dd = (short *) malloc((max_genome_length * 2 + 1) * sizeof(short));
  if (dd == NULL)
  {
    do_error("write_dst failed: couldn't allocate memory");
  }
  else
  {
    compress = (short *) malloc((max_genome_length * 2 + 1) * sizeof(short));
    if (compress == NULL)
    {
      do_error("write_dst failed: couldn't allocate memory");
    }
    else
    {
      /* printf("write_dst: starting\n"); */
      for (i = 0; i <= max_genome_length * 2; i++)
      {
        dd[i] = 0;
      }
      /* printf("write_dst: dd[] initialized to zero\n"); */
      for (i = 0; i < 100; i++)
      {
        dd_rel[i] = 0;
      }
      for (i = 0; i < psize - 1; i++)
      {
        for (j = i + 1; j < psize; j++)
        {
          d = edit_distance(lnd_genome[i].genome.length, (char *) lnd_genome[i].genome.g, lnd_genome[j].genome.length, (char *) lnd_genome[j].genome.g);
          /* printf("d(%lu, %lu) = %ld\n", i, j, d); */
          if (d >= 0)
          {
            dd[d]++;
            max_distance = (d > max_distance) ? d : max_distance;
            max_len = (lnd_genome[i].genome.num_genes > lnd_genome[j].genome.num_genes) ? lnd_genome[i].genome.num_genes : lnd_genome[j].genome.num_genes;
            if (max_len > 0)
              d = (d * 50) / max_len;
            else
              d = 0;
            dd_rel[(d < 100) ? d : 99]++;
          }
        }
      }
      /* printf("write_dst: dd[] filled, max_distance=%ld\n", max_distance); */
      compress_len = compress_histogram(max_genome_length * 2 + 1, dd, compress);
      /* printf("write_dst: dd[] compressed, length now: %lu\n", compress_len); */
      fprintf(dst_file, "%lu\r\n", compress_len);
      fwrite_int16array(compress, compress_len, dst_file);
      distance_entropy = shannon_short(max_distance, dd);
      distance_entropy_rel = shannon_short(100, dd_rel);
      /* printf("s = %lf, rel. s = %lf\n", distance_entropy, distance_entropy_rel); */
      free((void *) compress);
    }
    free((void *) dd);
  }
}

/*
PROCEDURE bedau_and_packard
  '
  ' **************************************************************************
  ' * determines the Bedau and Packard evolutionary activity in a generation *
  ' * and saves it.                                                          *
  ' **************************************************************************
*/

void write_bpe(void)
{
  unsigned long i, j;
  unsigned long compress_len;
  unsigned long bp_max = 0;

  short *bp, *compress;

  bp = (short *) malloc((bp_ceiling + 1) * sizeof(short));
  if (bp == NULL)
  {
    do_error("write_bpe failed: couldn't allocate memory");
  }
  else
  {
    compress = (short *) malloc((bp_ceiling + 1) * sizeof(short));
    if (compress == NULL)
    {
      do_error("write_bpe failed: couldn't allocate memory");
    }
    else
    {
      /* printf("write_bpe: starting\n"); */
      for (i = 0; i <= bp_ceiling; i++)
      {
        bp[i] = 0;
      }
      for (i = 0; i < psize; i++)
      {
        /* printf("write_bpe: processing genome %lu\n", i); */
        /* printf("write_bpe: length of genome: %lu\n", genome[i].len); */
        for (j = 0; j < lnd_genome[i].genome.num_genes; j++)
        {
          /* !!! difference form LNDEVL15.GFA: omitted range checking !!! */
          /* !!! because it is already done in growth(). !!! */
          bp_max = (lnd_genome[i].genome.bp_count[j] > bp_max) ? lnd_genome[i].genome.bp_count[j] : bp_max;
          bp[lnd_genome[i].genome.bp_count[j]]++;
        }
      }
      /* printf("write_bpe: loop done\n"); */
      compress_len = compress_histogram(bp_max + 1, bp, compress);
      /* printf("write_bpe: histogram compressed\n"); */
      fprintf(bpe_file, "%lu\r\n", compress_len);
      fwrite_int16array(compress, (size_t) compress_len, bpe_file);
      /* printf("write_bpe: %lu elements (%lu bytes) written\n", compress_len, compress_len * sizeof(compress[0])); */
      free((void *) compress);
    }
    free((void *) bp);
  }
}


static void write_gsys_parameters(FILE *f, const GSYS_PARAMETERS *gp)
{
  fprintf(f, "%ld\n", gp->num_divide);
  fprintf(f, "%ld\n", gp->num_flyingseed);
  fprintf(f, "%ld\n", gp->num_localseed);
  fprintf(f, "%ld\n", gp->num_mutminus);
  fprintf(f, "%ld\n", gp->num_mutplus);
  fprintf(f, "%ld\n", gp->num_statebit);
}


static void read_gsys_parameters(FILE *f, GSYS_PARAMETERS *gp)
{
  char buf[MAX_SLEN + 1];
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
  gp->num_statebit = strtol(buf, (char **) NULL, 10);
}


static void write_plant(FILE *f, long plant_no)
{
  long cell_no;

  write_genome(f, &(lnd_genome[gi[plant_no]].genome), 0xffffffff);
  gn_save_id(&(lnd_genome[gi[plant_no]].node_id), f);
  fprintf(f, "%lu\n", plant[plant_no].num_cells);
  for (cell_no = 0; cell_no < plant[plant_no].num_cells; cell_no++)
  {
    fprintf(f, "%ld\n", plant[plant_no].cell[cell_no].x);
    fprintf(f, "%ld\n", plant[plant_no].cell[cell_no].y);
    fprintf(f, "%d\n", plant[plant_no].cell[cell_no].energy);
  }
}


static int read_plant(FILE *f, long plant_no)
{
  char buf[MAX_SLEN + 1];
  long cell_no;

  read_genome(f, &(lnd_genome[gi[plant_no]].genome), 0xffffffff);
  gn_read_id(&(lnd_genome[gi[plant_no]].node_id), f);
  fgets(buf, MAX_SLEN, f);
  plant[plant_no].num_cells = strtol(buf, (char **) NULL, 10);
  for (cell_no = 0; cell_no < plant[plant_no].num_cells; cell_no++)
  {
    fgets(buf, MAX_SLEN, f);
    plant[plant_no].cell[cell_no].x = strtol(buf, (char **) NULL, 10);
    fgets(buf, MAX_SLEN, f);
    plant[plant_no].cell[cell_no].y = strtol(buf, (char **) NULL, 10);
    if (ferror(f) || (fgets(buf, MAX_SLEN, f) == NULL))
      return (-1);
    plant[plant_no].cell[cell_no].energy = strtol(buf, (char **) NULL, 10);
    world[plant[plant_no].cell[cell_no].x][plant[plant_no].cell[cell_no].y].plant_no = plant_no;
    world[plant[plant_no].cell[cell_no].x][plant[plant_no].cell[cell_no].y].cell_no = cell_no;
  }
  return (0);
}


int write_named_savefile(const char *savefile_name)
{
  FILE *f;
  long i;

  f = fopen(savefile_name, "w");
  if (f == NULL)
    return (-1);
  fprintf(f, "LindEvol-1.20 savefile\n");
  fprintf(f, "%s\n", simname);
  fprintf(f, "%ld\n", generation);
  fprintf(f, "%ld\n", psize);
  fprintf(f, "%f\n", srate);
  fprintf(f, "%1.12g\n", m_replacement);
  fprintf(f, "%1.12g\n", m_insertion);
  fprintf(f, "%1.12g\n", m_deletion);
  fprintf(f, "%1.12g\n", m_factor);
  fprintf(f, "%ld\n", world_width);
  fprintf(f, "%ld\n", world_height);
  fprintf(f, "%ld\n", num_days);
  fprintf(f, "%ld\n", glen_init);
  fprintf(f, "%ld\n", num_generations);
  fprintf(f, "%ld\n", bp_ceiling);
  fprintf(f, "%ld\n", ddistr_freq);
  fprintf(f, "%ld\n", phhist_freq);
  write_gsys_parameters(f, &gsys_parameters);
  fprintf(f, "%ld\n", random_seed);
  write_rndgenerator_state(f);
  write_genome(f, &old_main_species, 0);
  gn_save_tree(&gn_tree, f);
  fprintf(f, "hack: end saved tree\n");
  for (i = 0; i < psize; i++)
  {
    fprintf(f, "%lu\n", gi[i]);
    fprintf(f, "%lu\n", gi_tmp[i]);
    fprintf(f, "%lu\n", pl_index[i]);
    write_plant(f, i);
  }
  fclose(f);
  return (0);
}


int load_named_savefile(const char *savefile_name)
{
  FILE *f;
  char buf[MAX_SLEN + 1];
  long i, x, y;

  if ((f = fopen(savefile_name, "r")) == NULL)
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
  psize = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  srate = strtod(buf, (char **) NULL);
  num_select = bas_int(srate * psize);
  fgets(buf, MAX_SLEN, f);
  m_replacement = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  m_insertion = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  m_deletion = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  m_factor = strtod(buf, (char **) NULL);
  fgets(buf, MAX_SLEN, f);
  world_width = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  world_height = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  num_days = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  glen_init = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  num_generations = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  bp_ceiling = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  ddistr_freq = strtol(buf, (char **) NULL, 10);
  fgets(buf, MAX_SLEN, f);
  phhist_freq = strtol(buf, (char **) NULL, 10);
  read_gsys_parameters(f, &gsys_parameters);
  fgets(buf, MAX_SLEN, f);
  random_seed = strtol(buf, (char **) NULL, 10);
  fprintf(stderr, "random seed was %ld\n", random_seed);
  read_rndgenerator_state(f);
  if (read_genome(f, &old_main_species, 0) < 0)
  {
    do_error("load_savefile: error reading old main species");
    clear_arrays();
    fclose(f);
    return (-1);
  }
  if (gn_read_tree(&gn_tree, f))
    fprintf(stderr, "Problems reading genealogy tree\n");
  fgets(buf, MAX_SLEN, f);
  fprintf(stderr, "%s", buf);
  if (create_arrays() < 0)
  {
    do_error("load_savefile: out of memory");
    clear_arrays();
    fclose(f);
    return (-1);
  }
  for (y = 0; y < world_height; y++)
  {
    for (x = 0; x < world_width; x++)
    {
      world[x][y].plant_no = -1;
    }
  }
  for (i = 0; i < psize; i++)
  {
    fgets(buf, MAX_SLEN, f);
    gi[i] = strtol(buf, (char **) NULL, 10);
    fgets(buf, MAX_SLEN, f);
    gi_tmp[i] = strtol(buf, (char **) NULL, 10);
    fgets(buf, MAX_SLEN, f);
    pl_index[i] = strtol(buf, (char **) NULL, 10);
    if (read_plant(f, i) < 0)
      break;
  }
  fclose(f);
  if (i == psize)
    return (0);
  else
    return (-1);
}


int write_savefile(void)
{
  return (write_named_savefile(save_file_name));
}


int load_savefile(void)
{
  prepare_filenames();
  return (load_named_savefile(save_file_name));
}


int open_savetime_file(const char *fname)
{
  strncpy(savetime_file_name, fname, MAX_SLEN);
  if ((savetime_file = fopen(savetime_file_name, "r")) == NULL)
    return (-1);
  return (0);
}


void close_savetime_file(void)
{
  if (savetime_file)
    fclose(savetime_file);
  savetime_file = NULL;
}


long savetime_next(long generation)
{
  char buf[MAX_SLEN], e[MAX_SLEN * 2], *s;
  size_t i;
  long n;

  if (savetime_file == NULL)
    return (-1);
  for (;;)
  {
    do
      fgets(buf, MAX_SLEN, savetime_file);
    while (!feof(savetime_file) && !ferror(savetime_file) && (buf[0] == '\n'));
    if (feof(savetime_file) || ferror(savetime_file))
      break;
    i = strlen(buf);
    if (i && (buf[i - 1] == '\n'))
      buf[i - 1] = '\0';
    n = strtol(buf, &s, 10);
    if (buf == s)
    {
      sprintf(e, "savetime_next: Illegal value \"%s\"", buf);
      do_error(e);
    }
    else if (n <= generation)
    {
      sprintf(e, "savetime_next: Skipping \"%s\", already at generation %ld", buf, generation);
      do_error(e);
    }
    else
      return (n);
  }
  close_savetime_file();
  return (-1);
}

