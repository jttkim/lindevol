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

#include "jklib.h"
#include "lndtypes.h"
#include "lndio.h"
#include "creatarr.h"
#include "lndrandm.h"
#include "lnddispl.h"
#include "lnderror.h"
#include "pixfile.h"

#include "genomelib.h"
#include "gnlib.h"

static char tmpstr[MAX_SLEN];

#define LND_MAIN
#include "lndglobl.h"


#define DO_NEWSIM     1
#define DO_RESUMESIM  2
#define DO_EXIT       3
#define DO_TESTSIM    4


int main_menu(void)
{
  printf("Please choose\n");
  printf("1 - start new simulation\n");
  printf("2 - resume a simulation\n");
  printf("3 - exit program\n");
  printf("4 - test\n");
  fgets(tmpstr, MAX_SLEN, stdin);
  return ((int) strtol(tmpstr, (char **) NULL, 10));
}


int get_control_parameters(void)
{
#ifdef RANDOM_FITNESS
  simname[0] = 'c';
  simname[1] = '-';
  (void) input_str("name of simulation run: ", simname + 2, MAX_SLEN - 2);
#else
#ifdef FUJI_FITNESS
  simname[0] = 'f';
  simname[1] = '-';
  (void) input_str("name of simulation run: ", simname + 2, MAX_SLEN - 2);
#else
  (void) input_str("name of simulation run: ", simname, MAX_SLEN);
#endif
#endif
  prepare_filenames();
  psize = input_ul("population size: ");
  srate = input_dbl("selection rate: ");
  num_select = bas_int(srate * psize);           /* compute number of selected individuals */
  recombination = input_dbl("recombination rate: ");
  num_recombine = bas_int(recombination * psize);
  m_replacement =   input_dbl("replacement rate: ");
  m_insertion = input_dbl("insertion rate  : ");
  m_deletion = input_dbl("deletion rate   : ");
  m_factor =   input_dbl("mutation factor : ");

  gsys_parameters.num_divide = input_long("number of code bytes for divide: ");
  gsys_parameters.num_statebit = 0;
  gsys_parameters.num_broadcastbit = 0;
  gsys_parameters.num_flyingseed = 0;
  gsys_parameters.num_localseed = 0;
  gsys_parameters.num_mutminus = input_long("number of code bytes for mut-: ");
  gsys_parameters.num_mutplus = input_long("number of code bytes for mut+: ");

  random_seed = input_long("random seed: ");

  world_width = input_ul("world width: ");
  world_height = input_ul("world height: ");
  plant_distance = world_width / psize;
  num_days = input_ul("number of days per generation: ");
  glen_init = input_ul("initial genome length (in genes): ");

  num_generations = input_ul("number of generations to simulate: ");
  dmt_savefreq = input_ul("distance matrix every (0 for none): ");
  bp_ceiling = input_ul("ceiling for B&P analysis (0 for none): ");
  ddistr_freq = input_ul("save distance distributions every (0 for none): ");
  phhist_freq = input_ul("save phylogenetic trees every (0 for none): ");

  return (1);
}


int get_test_parameters(void)
{
#ifdef RANDOM_FITNESS
  strcpy(simname, "c-test");
#else
#ifdef FUJI_FITNESS
  strcpy(simname, "f-test");
#else
  strcpy(simname, "test");
#endif
#endif
  prepare_filenames();
  psize = 10;
  srate = 0.8;
  num_select = bas_int(srate * psize);
  recombination = 0.0;
  num_recombine = bas_int(recombination * psize);
  m_replacement = 0.03;
  m_insertion = 0.01;
  m_deletion = 0.01;
  m_factor = 1.0;

  gsys_parameters.num_divide = 256;
  gsys_parameters.num_statebit = 0;
  gsys_parameters.num_broadcastbit = 0;
  gsys_parameters.num_flyingseed = 0;
  gsys_parameters.num_localseed = 0;
  gsys_parameters.num_mutminus = 0;
  gsys_parameters.num_mutplus = 0;

  random_seed = 1;

  world_width = 30;
  world_height = 10;
  plant_distance = world_width / psize;
  num_days = 10;
  glen_init = 20;

  num_generations = 100;
  dmt_savefreq = 0;
  bp_ceiling = 0;
  ddistr_freq = 0;
  phhist_freq = 0;

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

#define par_id_ulong(pname, pvar) \
  if ((s = identify_parameter(line, pname )) != NULL) \
  { \
    pvar = strtoul(s, (char **) NULL, 10); \
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

#ifdef RANDOM_FITNESS
  if ((s = identify_parameter(line, "simname")) != NULL)
  {
    simname[0] = 'c';
    simname[1] = '-';
    i = 0;
    while (!iscntrl(s[i]))
    {
      simname[i + 2] = s[i];
      i++;
    }
    simname[i + 2] = '\0';
    return (0);
  }
#else
#ifdef FUJI_FITNESS
  if ((s = identify_parameter(line, "simname")) != NULL)
  {
    simname[0] = 'f';
    simname[1] = '-';
    i = 0;
    while (!iscntrl(s[i]))
    {
      simname[i + 2] = s[i];
      i++;
    }
    simname[i + 2] = '\0';
    return (0);
  }
#else
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
#endif
#endif
  par_id_long("psize", psize);
  par_id_double("srate", srate);
  par_id_double("recombination", recombination);
  par_id_double("m_replacement", m_replacement);
  par_id_double("m_insertion", m_insertion);
  par_id_double("m_deletion", m_deletion);
  par_id_double("m_factor", m_factor);

  par_id_long("num_divide", gsys_parameters.num_divide);
  gsys_parameters.num_statebit = 0;
  gsys_parameters.num_broadcastbit = 0;
  gsys_parameters.num_flyingseed = 0;
  gsys_parameters.num_localseed = 0;
  par_id_long("num_mutminus", gsys_parameters.num_mutminus);
  par_id_long("num_mutplus", gsys_parameters.num_mutplus);

  par_id_long("random_seed", random_seed);

  par_id_long("world_width", world_width);
  par_id_long("world_height", world_height);
  par_id_long("num_days", num_days);
  par_id_long("glen_init", glen_init);

  par_id_long("num_generations", num_generations);
  par_id_long("dmt_savefreq", dmt_savefreq);
  par_id_long("bp_ceiling", bp_ceiling);
  par_id_long("ddistr_freq", ddistr_freq);
  par_id_long("phhist_freq", phhist_freq);
  return (-1);
}


int load_control_parameters(const char *par_fname)
{
  FILE *f;
  char  errmsg[MAX_SLEN + 15], buf[MAX_SLEN];
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
  prepare_filenames();
  if (num_generations == 0)
    num_generations = input_long("number of generations to simulate: ");
  num_select = bas_int(srate * psize);
  num_recombine = bas_int(recombination * psize);
  plant_distance = world_width / psize;
  return (0);
}


/* init_genomes() initializes the genomes of a population with random
   strings. The usage and B & P counters ara initialized to 0 */

void init_genomes(unsigned long psize, LND_GENOME *lnd_genome, unsigned long len)
{
  unsigned long i, j;
  int ret_code;

  gn_init_tree(&gn_tree);
  for (i = 0; i < psize; i++)
  {
    if (alloc_genome(&(lnd_genome[i].genome), len * 2, len, GNM_USG | GNM_BP) < 0)
    {
      fprintf(stderr, "init_genomes: Failed to allocate genome #%lu\n", i);
      exit (EXIT_FAILURE);
    }
    /* printf("init_genomes: genome %lu has length %lu\n", i, genome[i].len); */
    for (j = 0; j < len; j++)
    {
      lnd_genome[i].genome.g[j * 2] = (unsigned char) lnd_random(256);
      lnd_genome[i].genome.g[j * 2 + 1] = (unsigned char) lnd_random(256);
      lnd_genome[i].genome.bp_count[j] = 0L;
    }

#ifdef testgenome

    lnd_genome[i].genome.g[0] = 0;
    lnd_genome[i].genome.g[1] = 5;
    lnd_genome[i].genome.g[2] = 4;
    lnd_genome[i].genome.g[3] = 5;

#endif

    if ((ret_code = gn_new_treenode(&gn_tree, NULL, -1, NULL, NULL, &(lnd_genome[i].node_id))) < 0)
    {
      fprintf(stderr, "*** error code %d creating node ID for initial genome %lu\n", ret_code, i);
    }
  }
  init_genome(&old_main_species);
}


void shuffle_population(void)
{
  unsigned long i, r1, gi_r1;

  for (i = 0; i < psize; i++)
  {
    r1 = lnd_random(psize);
    /* fprintf(stderr, "shuffle_population: swapping %lu with %lu\n", i, r1); */
    gi_r1 = gi[r1];
    gi[r1] = gi[i];
    gi[i] = gi_r1;
  }
/*
  for (i = 0; i < psize; i++)
    fprintf(stderr, "%3lu: %3lu\n", i, gi[i]);
*/
}


void init_plants(void)
{
  unsigned long i, j;
  long          x, y;

  for (y = 0; y < world_height; y++)
  {
    for (x = 0; x < world_width; x++)
    {
      world[x][y].plant_no = -1;
    }
  }
  for (i = 0; i < psize; i++)
  {
    plant[i].num_cells = 1;
    plant[i].cell[0].x = ((long) (world_width / psize)) * i;
    plant[i].cell[0].y = 0;
    plant[i].cell[0].energy = 0;
    world[plant[i].cell[0].x][plant[i].cell[0].y].plant_no = i;
    world[plant[i].cell[0].x][plant[i].cell[0].y].cell_no = 0;
    for (j = 0; j < lnd_genome[gi[i]].genome.num_genes; j++)
      lnd_genome[gi[i]].genome.usg_count[j] = 0;
    draw_cell(plant[i].cell[0].x, plant[i].cell[0].y);
  }
}


void sunshine()
{
  long x, y;

  for (x = 0; x < world_width; x++)
  {
    for (y = world_height - 1; y >= 0; y--)
    {
      if ((world[x][y].plant_no > -1) && (lnd_random(2) == 1))
      {
        plant[world[x][y].plant_no].cell[world[x][y].cell_no].energy = 1;
        draw_cell(x, y);
        break;
      }
    }
  }
}


void divide(unsigned long plant_no, unsigned long cell_no, int direction)
{
  long x, y;

  /* printf("divide: dividing plant %lu, cell %lu (%lu/%lu) to direction %d\n", plant_no, cell_no, plant[plant_no].cell[cell_no].x, plant[plant_no].cell[cell_no].y, direction); */
/*
  if (cell_no >= plant[plant_no].num_cells)
  {
    printf("  error: dividing non-existent cell #%lu of %lu\n", cell_no, plant[plant_no].num_cells);
  }
*/
  plant[plant_no].cell[cell_no].energy = 0;
  draw_cell(plant[plant_no].cell[cell_no].x, plant[plant_no].cell[cell_no].y);
  if (plant[plant_no].num_cells < MAX_NUM_CELLS)
  {
    x = (plant[plant_no].cell[cell_no].x + x_offset[direction] + world_width) % world_width;
    y = plant[plant_no].cell[cell_no].y + y_offset[direction];
    /* printf("new location: (%ld/%ld) = (%lu + %d/%lu + %d)\n", x, y, plant[plant_no].cell[cell_no].x, x_offset[direction], plant[plant_no].cell[cell_no].y, y_offset[direction]); */
    if ((y >= 0) && (y < world_height))
    {
      if (world[x][y].plant_no == -1)
      {
        /* printf("divide: dividing to (%ld/%ld)\n", x, y); */
        plant[plant_no].cell[plant[plant_no].num_cells].x = x;
        plant[plant_no].cell[plant[plant_no].num_cells].y = y;
        plant[plant_no].cell[plant[plant_no].num_cells].energy = 0;
        world[x][y].plant_no = plant_no;
        world[x][y].cell_no = plant[plant_no].num_cells;
        plant[plant_no].num_cells++;
        draw_cell(x, y);
      }
/*
      else
      {
        printf("divide: site (%ld/%ld) is occupied by plant %ld\n", x, y, world[x][y].plant_no);
      }
*/
    }
  }
  else
  {
    do_error("max. number of cells reached");
  }
}


unsigned char cell_state(unsigned long plant_no, unsigned long cell_no)
{
  unsigned char cs = 0;
  long x, y;
  int i;

  for (i = 0; i < 8; i++)
  {
    x = (plant[plant_no].cell[cell_no].x + x_offset[i] + world_width) % world_width;
    y = plant[plant_no].cell[cell_no].y + y_offset[i];
    if ((y >= 0) && (y < world_height))
    {
      if (world[x][y].plant_no == plant_no)
      {
        cs |= (1 << i);
      }
    }
  }
  return (cs);
}


int lnd_output(unsigned long genome_no, unsigned char input, unsigned long *gene_no)
{
  for (*gene_no = 0; *gene_no < lnd_genome[genome_no].genome.length; *gene_no += 2)
  {
    if (lnd_genome[genome_no].genome.g[*gene_no] == input)
    {
      *gene_no /= 2;
      return (lnd_genome[genome_no].genome.g[*gene_no * 2 + 1]);
    }
  }
  return (-1);
}


void plant_growth(unsigned long plant_no)
{
  unsigned long cell_no, gene_no;
  unsigned char input;
  int           output;
  long          action;

  lnd_genome[gi[plant_no]].genome.mut_flag = 0;
  for (cell_no = 0 ; cell_no < plant[plant_no].num_cells; cell_no++)
  {
    if (plant[plant_no].cell[cell_no].energy)
    {
      input = cell_state(plant_no, cell_no);
      output = lnd_output(gi[plant_no], input, &gene_no);
      if (output > -1)
      {
	action = gene_activity((unsigned char) output);
	if (action == LND_DIVIDE)
	  divide(plant_no, cell_no, output & 7);
	else if (action == LND_MUTMINUS)
	{
	  lnd_genome[gi[plant_no]].genome.mut_flag--;
	  plant[plant_no].cell[cell_no].energy = 0;
	}
	else if (action == LND_MUTPLUS)
	{
	  lnd_genome[gi[plant_no]].genome.mut_flag++;
	  plant[plant_no].cell[cell_no].energy = 0;
	}
        lnd_genome[gi[plant_no]].genome.usg_count[gene_no]++;
        if (lnd_genome[gi[plant_no]].genome.bp_count[gene_no] < bp_ceiling)
        {
          lnd_genome[gi[plant_no]].genome.bp_count[gene_no]++;
        }
      }
    }
  }
}


unsigned long num_energyrich_cells(unsigned long plant_no)
{
  unsigned long cell_no;
  unsigned long f = 0;

  for (cell_no = 0; cell_no < plant[plant_no].num_cells; cell_no++)
  {
    if (plant[plant_no].cell[cell_no].energy)
    {
      f++;
    }
  }
  /* printf("num_energyrich_cells: plant %lu has %lu cells, f=%lu\n", plant_no, plant[plant_no].num_cells, f); */
  return f;
}


void shuffle_pl_index()
{
  unsigned long i, r, tmp;
  for (i = 0; i < psize; i++)
  {
    r = lnd_random(psize);
    tmp = pl_index[r];
    pl_index[r] = pl_index[i];
    pl_index[i] = tmp;
  }
}


/* Note: In fitness(), plant number i has the genome (*genome[gi[i]]).
   In a given day, it is processed as the pl_index[i]-th plant. */

void fitness(void)
{
  unsigned long day;
  unsigned long i;

  init_plants();
  for (day = 0; day < num_days; day++)
  {
    sunshine();
    shuffle_pl_index();
    for (i = 0; i < psize; i++)
    {
      plant_growth(pl_index[i]);
    }
    if (pixmode)
      write_pixfile();
  }

#ifdef RANDOM_FITNESS

  /* printf("Assigning random fitness values\n"); */
  for (i = 0; i < psize; i++)
  {
    lnd_genome[gi[i]].fitness = lnd_random(1000);
  }

#else
#ifdef FUJI_FITNESS

  for (i = 0; i < psize; i++)
  {
    size_t j;

    lnd_genome[gi[i]].fitness = 0;
    for (j = 0; j < lnd_genome[gi[i]].genome.length; j++)
      lnd_genome[gi[i]].fitness += lnd_genome[gi[i]].genome.g[j];
  }

#else

  for (i = 0; i < psize; i++)
  {
    lnd_genome[gi[i]].fitness = num_energyrich_cells(i);
  }

#endif
#endif

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
  unsigned long fitness_sum = 0;
  unsigned long glen_sum = 0;
  unsigned long num_used_sum = 0;
  unsigned long num_used;
  unsigned long i, j;
  long mutflag_sum = 0;

  min_genome_length = lnd_genome[0].genome.num_genes;
  max_genome_length = min_genome_length;
  min_num_used = ULONG_MAX;
  max_num_used = 0;
  min_mutflag = lnd_genome[0].genome.mut_flag;
  max_mutflag = min_mutflag;
  for (i = 0; i < psize; i++)
  {
    num_used = 0;
    for (j = 0; j < lnd_genome[i].genome.num_genes; j++)
    {
      if (lnd_genome[i].genome.usg_count[j] > 0)
      {
        num_used++;
      }
    }
    fitness_sum += lnd_genome[i].fitness;
    glen_sum += lnd_genome[i].genome.num_genes;
    min_genome_length = (min_genome_length < lnd_genome[i].genome.num_genes) ? min_genome_length : lnd_genome[i].genome.num_genes;
    max_genome_length = (max_genome_length > lnd_genome[i].genome.num_genes) ? max_genome_length : lnd_genome[i].genome.num_genes;
    num_used_sum += num_used;
    min_num_used = (min_num_used < num_used) ? min_num_used : num_used;
    max_num_used = (max_num_used > num_used) ? max_num_used : num_used;
    min_mutflag = (lnd_genome[i].genome.mut_flag < min_mutflag) ? lnd_genome[i].genome.mut_flag : min_mutflag;
    max_mutflag = (lnd_genome[i].genome.mut_flag > max_mutflag) ? lnd_genome[i].genome.mut_flag : max_mutflag;
    mutflag_sum += lnd_genome[i].genome.mut_flag;
  }
  average_fitness = (double) fitness_sum / psize;
  average_genome_length = (double) glen_sum / psize;
  average_num_used = (double) num_used_sum / psize;
  average_mutflag = (double) mutflag_sum / psize;
}


unsigned long next_mutpos(double m, unsigned long pos)
{
  double p;
  double r, p1, p2;

  if (m >= 1.0 - DBL_EPSILON)
    return (pos);
  r = lnd_rnd();
  errno = 0;
  p2 = log(1.0 - m);
  if (errno)
  {
    fprintf(stderr, "next_mutpos: error %d (%s) with log(1.0 - %f)\n", errno, strerror(errno), m);
    p2 = 0.0;
  }
  if (p2 == 0.0)
  {
    return (ULONG_MAX);
  }
  errno = 0;
  p1 = log(1.0 - r) / p2;
  if (errno)
  {
    fprintf(stderr, "next_mutpos: error %d (%s) with log(1.0 - %f) / %f\n", errno, strerror(errno), r, p2);
    if (errno != ERANGE)
    {
      perror("next_mutpos failed");
      /* printf("errno = %d\n", errno); */
    }
    return (ULONG_MAX);
  }
  errno = 0;
  p = floor(pos + p1);

#ifdef dreck

  if (errno)
  {
    printf("next_mutpos: error with floor %d, pos=%lu, p1=%f, p=%f, r=%f\n", errno, pos, p1, p, r);
    if (errno != ERANGE)
    {
      perror("next_mutpos failed");
      /* printf("errno = %d\n", errno); */
    }
    /* p = DBL_MAX; */
  }

#endif

  if (p < 0.0)
  {
    p = 0.0;
  }
  /* fprintf(stderr, "next_mutpos: p = %1.12g\n", p); */
  return ((unsigned long) (ULONG_MAX < p) ? ULONG_MAX : p);
}


void mutation(void)
{
  unsigned long i, j;
  unsigned long p;
  double mr, mi, md, mf;

  for (i = 0; i < psize; i++)
  {
    mf = pow(m_factor, lnd_genome[gi[i]].genome.mut_flag);
    mr = m_replacement * mf;
    mi = m_insertion * mf;
    md = m_deletion * mf;
    lnd_genome[gi[i]].genome.mut_flag = 0;
    if (mr > 0.0)
    {
      p = next_mutpos(mr, 0);
      while (p < lnd_genome[gi[i]].genome.length)
      {
        lnd_genome[gi[i]].genome.g[p] = (unsigned char) lnd_random(256);
	/* fprintf(stderr, "g[%ld] = %02x\n", p, lnd_genome[gi[i]].genome.g[p]); */
        lnd_genome[gi[i]].genome.bp_count[p / 2] = 0;
        p = next_mutpos(mr, p) + 1;
      }
    }
    if (mi > 0.0)
    {
      p = next_mutpos(mi, 0);
      while (p <= lnd_genome[gi[i]].genome.length)
      {
	if (resize_genome(&(lnd_genome[gi[i]].genome), lnd_genome[gi[i]].genome.length + 2, lnd_genome[gi[i]].genome.num_genes + 1) < 0)
        {
          do_error("mutation: Error resizing genome during insertion\n");
          break;
        }
        for (j = lnd_genome[gi[i]].genome.length - 1; p + 2 <= j; j--)
        {
          lnd_genome[gi[i]].genome.g[j] = lnd_genome[gi[i]].genome.g[j - 2];
        }
        lnd_genome[gi[i]].genome.g[p] = (unsigned char) lnd_random(256);
        lnd_genome[gi[i]].genome.g[p + 1] = (unsigned char) lnd_random(256);
        for (j = lnd_genome[gi[i]].genome.num_genes - 1; p / 2 + 1 <= j; j--)
        {
          lnd_genome[gi[i]].genome.bp_count[j] = lnd_genome[gi[i]].genome.bp_count[j - 1];
        }
        lnd_genome[gi[i]].genome.bp_count[p / 2] = 0;
        if ((p % 2) == 1)
        {
          lnd_genome[gi[i]].genome.bp_count[p / 2 + 1] = 0;
        }
        p = next_mutpos(mi, p) + 3;
      }
    }
    /* printf("  deleting at random, length = %lu\n", genome[gi[i]].len); */
    if ((lnd_genome[gi[i]].genome.length > 0) && (md > 0.0))
    {
      p = next_mutpos(md, 0);
      while ((lnd_genome[gi[i]].genome.length >= 2) && (p < (lnd_genome[gi[i]].genome.length - 2)))
      {
        for (j = p + 2; j < lnd_genome[gi[i]].genome.length; j++)
        {
          lnd_genome[gi[i]].genome.g[j - 2] = lnd_genome[gi[i]].genome.g[j];
        }
        for (j = p / 2 + 1; j < lnd_genome[gi[i]].genome.num_genes; j++)
        {
          lnd_genome[gi[i]].genome.bp_count[j - 1] = lnd_genome[gi[i]].genome.bp_count[j];
        }
        if ((p % 2) == 1)
        {
          lnd_genome[gi[i]].genome.bp_count[p / 2] = 0;
        }
	if (resize_genome(&(lnd_genome[gi[i]].genome), lnd_genome[gi[i]].genome.length - 2, lnd_genome[gi[i]].genome.num_genes - 1) < 0)
        {
          do_error("mutation: Error resizing genome during deletion\n");
          break;
        }
        p = next_mutpos(md, p);
      }
    }
/*
    fprintf(stderr, "genome %3ld: ", gi[i]);
    for (j = 0; j < lnd_genome[gi[i]].genome.length; j++)
      fprintf(stderr, "%02x", lnd_genome[gi[i]].genome.g[j]);
    fprintf(stderr, "\n");
*/
    for (j = 0; j < lnd_genome[gi[i]].genome.num_genes; j++)
      lnd_genome[gi[i]].genome.usg_count[j] = 0;
  }
}


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

int cmp_lexical(const void *gi1, const void *gi2)
{
  unsigned long g1;
  unsigned long *g1_p;
  unsigned long g2;
  unsigned long *g2_p;
  unsigned long i = 0;
  unsigned long l;

  g1_p = (unsigned long *) gi1;
  g2_p = (unsigned long *) gi2;
  g1 = *g1_p;
  g2 = *g2_p;
  /* printf("comparing genomes %lu and %lu\n", g1, g2); */
  l = (lnd_genome[g1].genome.length < lnd_genome[g2].genome.length) ? lnd_genome[g1].genome.length : lnd_genome[g2].genome.length;
  for (i = 0; i < l; i++)
  {
    if (lnd_genome[g1].genome.g[i] < lnd_genome[g2].genome.g[i])
    {
      return (-1);
    }
    else if (lnd_genome[g1].genome.g[i] > lnd_genome[g2].genome.g[i])
    {
      return (1);
    }
  }
  if (lnd_genome[g1].genome.length < lnd_genome[g2].genome.length)
    return (-1);
  else if (lnd_genome[g1].genome.length > lnd_genome[g2].genome.length)
    return (1);
  else
  {
    if (g1 < g2)
      return (-1);
    else
      return (1);
  }
}

/* cmp_species_d() returns -1 if s1 is greater than s2 and 1 if s1 is less than
   s2, in order to lead to descending qsorting of the species_d[] array */

int cmp_species_d(const void *s1, const void *s2)
{
  unsigned long n1 = ((SPECIES_D *) s1)->num;
  unsigned long n2 = ((SPECIES_D *) s2)->num;

  return (n1 < n2) ? 1 : ((n1 > n2) ? -1 : 0);
}

void count_species(void)
{
  unsigned long i;
  GENOME   *main_species;                      /* the main species */
  double f;

  /* printf("counting species\n"); */
  qsort((void *) gi_tmp, (size_t) psize, sizeof(unsigned long), &cmp_lexical);
  /* printf("sorting done\n"); */
  num_species = 0;
  species_d[0].num = 1;
  species_d[0].index = gi_tmp[0];
  for (i = 1; i < psize; i++)
  {
    if (cmp_lexical((void *) &(gi_tmp[i - 1]), (void *) &(gi_tmp[i])) == 0)
    {
      species_d[num_species].num++;
    }
    else
    {
      num_species++;
      species_d[num_species].num = 1;
      species_d[num_species].index = gi_tmp[i];
    }
  }
  num_species++;
  genetic_diversity = 0.0;
  rel_genetic_diversity = 0.0;
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
  main_species = &(lnd_genome[species_d[0].index].genome);
  if (generation == 0)
    mainspec_distance = 0;
  else
    mainspec_distance = edit_distance(main_species->length, (char *) main_species->g, old_main_species.length, (char *)old_main_species.g);
  if (mainspec_distance == -1)
  {
    do_error("error while computing distance");
  }
  if (generation > 0)
    free_genome(&old_main_species);
  duplicate_genome(&old_main_species, main_species, 0);
}


/*
PROCEDURE sortieren
  ' ************************************************************************
  ' * This procedure sorts the genomes according to their fitness values   *
  ' * such that the genomes with the greatest fitness are the first one in *
  ' * the g$() array.                                                      *
  ' ************************************************************************
*/


/*
 * cmp_fitness returns -1 or 1 randomly if fitness values are equal,
 * and thus never returns 0 (indicating compared elements are equal).
 * This hack is ensures that the sorting order is completely determined
 * by lindevol (and its random number generator) and not dependent on
 * possible peculiarities of randomizes qsort() implementations.
 */

int cmp_fitness(const void *gi1, const void *gi2)
{
  unsigned long i1 = *((unsigned long *) gi1), i2 = *((unsigned long *) gi2);
  unsigned long f1 = lnd_genome[i1].fitness;
  unsigned long f2 = lnd_genome[i2].fitness;

  /* fprintf(stderr, "comparing %lu (%lu) with %lu (%lu)\n", f1, *(unsigned long *) gi1, f2, *(unsigned long *) gi2); */
  if (f1 < f2)
    return (1);
  else if (f1 > f2)
    return (-1);
  else
  {
    if (gi_tmp[i1] < gi_tmp[i2])
      return (-1);
    else
      return (1);
  }
}

void sort_population(void)
{
  unsigned long i, gr, r;

  for (i = 0; i < psize; i++)
    gi_tmp[i] = i;
  for (i = 0; i < psize; i++)
  {
    r = lnd_random(psize);
    gr = gi_tmp[r];
    gi_tmp[r] = gi_tmp[i];
    gi_tmp[i] = gr;
  }
/*
  fprintf(stderr, "before sorting:\n");
  for (i = 0; i < psize; i++)
    fprintf(stderr, "%3lu: genome %3lu, fitness: %lu\n", i, gi[i], lnd_genome[gi[i]].fitness);
*/
  qsort((void *)gi, (size_t) psize, sizeof(unsigned long), cmp_fitness);
/*
  fprintf(stderr, "after sorting:\n");
  for (i = 0; i < psize; i++)
    fprintf(stderr, "%3lu: genome %3lu, fitness: %lu\n", i, gi[i], lnd_genome[gi[i]].fitness);
*/
}


/*
PROCEDURE ersetzen
  ' *******************************************************************
  ' * replace the genomes that are at the end of the array g$() after *
  ' * the array has been sorted according to the fitness values       *
  ' *******************************************************************
*/

void selection(void)
{
  unsigned long i, j;
  unsigned long r1, r2;
  unsigned long l;
  unsigned long p;
  int ret_code;

  for (i = 1; i <= num_recombine; i++)
  {
    /* !!! difference from method used in LNDEVL15.GFA !!! */
    r1 = lnd_random(psize - num_select);
    r2 = (r1 + lnd_random(psize - num_select - 1) + 1) % (psize - num_select);
    l = (lnd_genome[gi[r1]].genome.length < lnd_genome[gi[r2]].genome.length) ? lnd_genome[gi[r1]].genome.length : lnd_genome[gi[r2]].genome.length;
    p = lnd_random(l - 1);
    if ((ret_code = gn_node_death(&gn_tree, &(lnd_genome[gi[psize - i]].node_id), generation)) < 0)
    {
      fprintf(stderr, "*** error %d in gn_node_death() for genome #%lu\n", ret_code, gi[psize - i]);
    }
    free_genome(&(lnd_genome[gi[psize - i]].genome));
    if (duplicate_genome(&(lnd_genome[gi[psize - i]].genome), &(lnd_genome[gi[r1]].genome), GNM_USG | GNM_BP) < 0)
    {
      do_error("selection: Failed to create new genome during recombination");
      break;
    }
    if (resize_genome(&(lnd_genome[gi[psize - i]].genome), lnd_genome[gi[r2]].genome.length, lnd_genome[gi[r2]].genome.num_genes) < 0)
    {
      do_error("selection: Failed to resize newly created genome during recombination");
      break;
    }
    for (j = p + 1; j < lnd_genome[gi[r2]].genome.length; j++)
    {
      lnd_genome[gi[psize - i]].genome.g[j] = lnd_genome[gi[r2]].genome.g[j];
    }
    for (j = (p + 1) / 2; j < lnd_genome[gi[r2]].genome.num_genes; j++)
    {
      lnd_genome[gi[psize - i]].genome.bp_count[j] = lnd_genome[gi[r2]].genome.bp_count[j];
    }
    /* !!! correction from LNDEVL15.GFA: usage counter of recombined gene cleared !!! */
    if ((p % 2) == 0)
    {
      lnd_genome[gi[psize - i]].genome.bp_count[p / 2] = 0;
    }
    if (p < lnd_genome[gi[r2]].genome.length - p)
    {
      if ((ret_code = gn_new_treenode(&gn_tree, &(lnd_genome[gi[r2]].node_id), generation, NULL, NULL, &(lnd_genome[gi[psize - i]].node_id))) < 0)
      {
	fprintf(stderr, "*** error %d in gn_new_treenode() for genome #%lu\n", ret_code, gi[psize - i]);
      }
    }
    else
    {
      if ((ret_code = gn_new_treenode(&gn_tree, &(lnd_genome[gi[r1]].node_id), generation, NULL, NULL, &(lnd_genome[gi[psize - i]].node_id))) < 0)
      {
	fprintf(stderr, "*** error %d in gn_new_treenode() for genome #%lu\n", ret_code, gi[psize - i]);
      }
    }
  }
  for (i = num_recombine + 1; i <= num_select; i++)
  {
    r1 = lnd_random(psize - num_select);
    if ((ret_code = gn_node_death(&gn_tree, &(lnd_genome[gi[psize - i]].node_id), generation)) < 0)
    {
      fprintf(stderr, "*** error %d in gn_node_death() for genome #%lu\n", ret_code, gi[psize - i]);
    }
    free_genome(&(lnd_genome[gi[psize - i]].genome));
    if (duplicate_genome(&(lnd_genome[gi[psize - i]].genome), &(lnd_genome[gi[r1]].genome), GNM_USG | GNM_BP) < 0)
    {
      do_error("selection: Failed to create new genome during vegetative multiplication");
      break;
    }
    if ((ret_code = gn_new_treenode(&gn_tree, &(lnd_genome[gi[r1]].node_id), generation, NULL, NULL, &(lnd_genome[gi[psize - i]].node_id))) < 0)
    {
      fprintf(stderr, "*** error %d in gn_new_treenode() for genome #%lu\n", ret_code, gi[psize - i]);
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


/* evolution() is the main simulation engine. */

void evolution(void)
{
  long finish_generation, save_generation = -1;
  char buf[MAX_SLEN + 32];

  finish_generation = start_generation + num_generations;
  if (savetime_fname[0])
  {
    open_savetime_file(savetime_fname);
    save_generation = savetime_next(start_generation);
  }

  init_world_display(simname);
  for (generation = start_generation; generation < finish_generation; generation++)
  {
    count_species();           /* determine number of different species */
    if (bp_ceiling > 0)
    {
      write_bpe();             /* write B & P histogram for current generation */
    }
    if (phhist_freq && ((generation % phhist_freq) == 0))
    {
      write_stb();             /* write tree */
      write_genomes();
    }
    shuffle_population();      /* arrange gi[] in a random order */
    display_on(generation);    /* entertain the users ... ;-) */
    fitness();                 /* determine fitness values by growth simulation */
    display_off();
    statistics();              /* ... of fitness, genome length, gene usage */
    if (generation == save_generation)
    {
      sprintf(buf, "%s-%07ld.dat", simname, generation);
      write_named_savefile(buf);
      save_generation = savetime_next(generation);
    }
    sort_population();         /* sort gi[] to descending fitness values of genome[] */
    if (ddistr_freq && ((generation % ddistr_freq) == 0))
    {
      write_dst();             /* write the distance distribution */
    }
    write_pro();               /* write out the protocol for current generation */
    if (dmt_savefreq > 0)
    {
      if ((generation % dmt_savefreq) == 0)
      {
        write_dmt();           /* write distance matrix (if desired) */
      }
    }
    selection();               /* selection, assignment of ancestors happens here */
    mutation();                /* mutate genomes in population */
  }
  write_savefile();
  close_world_display();
  gn_free_tree(&gn_tree);
  free_genome(&old_main_species);
}


int get_modes(int argc, char *argv[])
{
  int i, ret_code = 0;

  par_fname[0] = '\0';
  savetime_fname[0] = '\0';
  for (i = 1; i < argc; i++)
  {
    if (!strcmp(argv[i], "-t"))
    {
      testmode = 1;
      continue;
    }
    if (!strcmp(argv[i], "-q"))
    {
      quietmode = 1;
      continue;
    }
    if (!strcmp(argv[i], "-w"))
    {
      worldmode = 1;
      continue;
    }
    if (!strcmp(argv[i], "-p"))
    {
      pixmode = 1;
      continue;
    }
    if (!strcmp(argv[i], "-f"))
    {
      if (i + 1 < argc)
      {
        i++;
        strncpy(par_fname, argv[i], MAX_SLEN);
        ret_code = DO_NEWSIM;
      }
      continue;
    }
    if (!strcmp(argv[i], "-s"))
    {
      if (i + 1 < argc)
      {
        i++;
        strncpy(savetime_fname, argv[i], MAX_SLEN);
      }
      continue;
    }
    if (!strcmp(argv[i], "-r"))
    {
      if (i + 1 < argc)
      {
        i++;
        strncpy(simname, argv[i], MAX_SLEN);
        ret_code = DO_RESUMESIM;
      }
      continue;
    }
    if (!strcmp(argv[i], "-h"))
    {
      printf("\n%s commandline options:\n\n", SIMPRGNAME);
      printf("-q        : start running in quiet mode\n");
      printf("-t        : run test simulation (with builtin test parameter set)\n");
      printf("-w        : create a world file\n");
      printf("-p        : create a pixel file\n");
      printf("-f <file> : read control parameters from <file>\n");
      printf("-r <name> : resume simulation <name>\n");
      printf("-s <file> : read generations at which to save state from <file>\n");
      exit (0);
    }
    fprintf(stderr, "\aCannot understand commandline argument %s\n", argv[i]);
  }
  return (ret_code);
}


int main(int argc, char *argv[])
{
  int what_now;

#ifdef RANDOM_FITNESS
  printf("lndcntrl, RANDOM FITNESS CONTROL for %s, by Jan T. Kim\n", SIMPRGNAME);
#else
#ifdef FUJI_FITNESS
  printf("lndcntrl, FUJI FITNESS CONTROL for %s, by Jan T. Kim\n", SIMPRGNAME);
#else
  printf("%s, written by Jan T. Kim\n", SIMPRGNAME);
#endif
#endif
  printf("compiled %s %s\n", __DATE__, __TIME__);
  par_fname[0] = '\0';
  what_now = get_modes(argc, argv);
  if (what_now == 0)
    what_now = main_menu();
  if ((what_now == DO_NEWSIM) || (what_now == DO_TESTSIM))
  {
    /* get_control_parameters(); */
    if (what_now == DO_NEWSIM)
    {
      if (par_fname[0])
      {
	if (load_control_parameters(par_fname) != 0)
	{
	  do_error("Failed to read parameter file -- exit\n");
	  exit(EXIT_FAILURE);
	}
      }
      else
	get_control_parameters();
    }
    else if (what_now == DO_TESTSIM)
    {
      get_test_parameters();
      what_now = DO_NEWSIM;
    }
    if (create_arrays() < 0)
    {
      do_error("Out of memory error");
      return (-1);
    }
    if (open_pro_file("w") != 0)
    {
      clear_arrays();
      return (-1);
    }
    if (dmt_savefreq > 0)
      if (open_dmt_file("w") != 0)
        dmt_savefreq = 0;
    if (phhist_freq)
    {
      if (open_stb_file("w") != 0)
        phhist_freq = 0;
      if (open_genome_file("w") != 0)
        phhist_freq = 0;
    }
    if (ddistr_freq)
      if (open_dst_file("w") != 0)
        ddistr_freq = 0;
    if (bp_ceiling > 0)
      if (open_bpe_file("w") != 0)
        bp_ceiling = 0;
    if (worldmode)
      open_world_file(simname);
    if (pixmode)
      open_pixfile("w");
    seed_lnd_random(random_seed);
    init_genomes(psize, lnd_genome, glen_init);
  }
  else if (what_now == 2)
  {
    if (load_savefile() < 0)
      return (EXIT_FAILURE);
    if (open_pro_file("a") != 0)
    {
      clear_arrays();
      return (EXIT_FAILURE);
    }
    if (dmt_savefreq > 0)
    {
      if (open_dmt_file("a") != 0)
      {
	close_data_files();
	clear_arrays();
	return (EXIT_FAILURE);
      }
    }
    if (phhist_freq)
    {
      if (open_stb_file("a") != 0)
      {
	close_data_files();
	clear_arrays();
	return (EXIT_FAILURE);
      }
      if (open_genome_file("a") != 0)
      {
	close_data_files();
	clear_arrays();
	return (EXIT_FAILURE);
      }
    }
    if (ddistr_freq)
    {
      if (open_dst_file("a") != 0)
      {
	close_data_files();
	clear_arrays();
	return (EXIT_FAILURE);
      }
    }
    if (bp_ceiling > 0)
    {
      if (open_bpe_file("w") != 0)
      {
	close_data_files();
	clear_arrays();
	return (EXIT_FAILURE);
      }
    }
    if (pixmode)
      open_pixfile("a");
    num_generations = input_long("number of generations to simulate: ");
  }
  if ((what_now == 1) || (what_now == 2))
  {
    calculate_activity_codes(&gsys_parameters);
    evolution();
    close_data_files();
    if (worldmode)
      close_world_file();
    clear_arrays();
  }
  return (0);
}

