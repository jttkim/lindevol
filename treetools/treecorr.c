#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __atarist__
#  include <unistd.h>
#else
#  include <getopt.h>
#endif

#include "gnlib.h"
#include "ptlib.h"
#include "genomelib.h"

#ifdef MEMDEBUG
#  include "memdebug.h"
#endif

typedef struct
{
  GN_NODE_ID node_id;
  char name[80];
  GENOME     genome;
} SEQ;

typedef struct
{
  long size;
  long max_len;
  SEQ *seq;
} POPULATION;


char buf[256];


char *get_line(char *buf, int l, FILE *f)
{
  char *s;
  size_t i;

  do
  {
    if ((s = fgets(buf, l, f)))
    {
      i = strlen(buf);
      if ((i > 0) && (buf[i - 1] == '\n'))
        buf[i - 1] = '\0';
    }
    if (s == NULL)
      break;
  }
  while (buf[0] == '\0');
  return (s);
}


long diffchar_distance(long n1, const unsigned char *s1, long n2, const unsigned char *s2)
{
  long l = (n1 < n2) ? n1 : n2;
  long i, d = 0;

  for (i = 0; i < l; i++)
    d += s1[i] != s2[i];
  return (d);

}


void free_population(POPULATION *pop)
{
  long i;

  for (i = 0; i < pop->size; i++)
    free_genome(&(pop->seq[i].genome));
  free(pop->seq);
}


int read_genomes(FILE *genomefile, POPULATION *pop)
{
  long i, j;

  get_line(buf, 256, genomefile);
  pop->size = strtol(buf, (char **) NULL, 10);
  if ((pop->seq = (SEQ *) malloc(pop->size * sizeof(SEQ))) == NULL)
    return (-1);
  for (i = 0; i < pop->size; i++)
  {
    if (gn_read_id(&(pop->seq[i].node_id), genomefile) < 0)
    {
      fprintf(stderr, "error reading genome ID %ld\n", i);
      for (j = 0; j < i; j++)
        free_genome(&(pop->seq[i].genome));
      free(pop->seq);
      return (-1);
    }
    gn_node_idstring(&(pop->seq[i].node_id), pop->seq[i].name);
    if (read_genome(genomefile, &(pop->seq[i].genome), 0) < 0)
    {
      fprintf(stderr, "error reading genome %ld\n", i);
      for (j = 0; j < i; j++)
        free_genome(&(pop->seq[i].genome));
      free(pop->seq);
      return (-1);
    }
  }
  pop->max_len = 0;
  for (i = 0; i < pop->size; i++)
  {
    if (pop->max_len < pop->seq[i].genome.length)
      pop-> max_len = pop->seq[i].genome.length;
  }
  return (0);
}


int cmp_seq_name(const void *s1, const void *s2)
{
  char *n1 = ((SEQ *) s1)->name;
  char *n2 = ((SEQ *) s2)->name;

  return (strcmp(n1, n2));
}


void sort_population(POPULATION *pop)
{
  qsort(pop->seq, pop->size, sizeof(SEQ *), cmp_seq_name);
}


int linear_regression(long num, const double *x, const double *y, double *a, double *b, double *r)
{
  long i;
  double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
  double x_av, y_av, xdiff, ydiff, sum_xdiff = 0.0, sum_ydiff = 0.0, sum_xydiff = 0.0;
  double det;

  for (i = 0; i < num; i++)
  {
    sum_x += x[i];
    sum_y += y[i];
    sum_xx += x[i] * x[i];
    sum_xy += x[i] * y[i];
  }
  det = num * sum_xx - sum_x * sum_x;
  *a = (num * sum_xy - sum_x * sum_y) / det;
  *b = (sum_xx * sum_y - sum_xy * sum_x) / det;
  x_av = sum_x / num;
  y_av = sum_y / num;
  for (i = 0; i < num; i++)
  {
    xdiff = x[i] - x_av;
    ydiff = y[i] - y_av;
    sum_xydiff += xdiff * ydiff;
    sum_xdiff += xdiff * xdiff;
    sum_ydiff += ydiff * ydiff;
  }
  *r = sum_xydiff / sqrt(sum_xdiff) / sqrt(sum_ydiff);
  return (0);
}


int main(int argc, char **argv)
{
  extern char *optarg;
  int          optchar;
  char *treefile_name = NULL, *genomefile_name = NULL, *outfile_name = NULL, *corrfile_name = NULL, *mutfile_name = NULL;
  FILE *treefile, *genomefile, *outfile, *corrfile, *mutfile;
  long genome_generation, tree_generation, i, j, num_trees;
  double tree_d, rec_d, rec_d1, change_rate, mutation_rate;
  double *tree_dist = NULL, *reconst_dist = NULL;
  double lambda, y0, r;
  long num_data;
  PHYLTREE phyltree;
  POPULATION population;
  int ret_code;

  phyl_init_tree(&phyltree);
  while ((optchar = getopt(argc, argv, "hg:t:o:c:m:")) != -1)
  {
    switch (optchar)
    {
    case 'g':
      genomefile_name = optarg;
      break;
    case 't':
      treefile_name = optarg;
      break;
    case 'c':
      corrfile_name = optarg;
      break;
    case 'm':
      mutfile_name = optarg;
      break;
    case 'o':
      outfile_name = optarg;
      break;
    case 'h':
      printf("treecorr -- check correlation between corrected\n");
      printf("    edit or hamming diatance and true tree distance\n");
      printf("\n");
      printf("Command line usage:\n");
      printf("-g <filename>: Specify genome file (mandatory)\n");
      printf("-t <filename>: Specify tree file (mandatory)\n");
      printf("-o <filename>: Specify base name for individual generation's output file\n");
      printf("-c <filename>: Specify file for treedist vs reconstructed dist correlations\n");
      printf("-m <filename>: Specify file for reconstructed mutation rates\n");
      printf("-h: Print this help and exit\n");
      exit (EXIT_SUCCESS);
    }
  }
  if (treefile_name == NULL)
  {
    fprintf(stderr, "no tree file specified -- exit\n");
    exit (EXIT_FAILURE);
  }
  if ((treefile = fopen(treefile_name, "r")) == NULL)
  {
    fprintf(stderr, "failed to open tree file \"%s\" -- exit\n", treefile_name);
    exit (EXIT_FAILURE);
  }
  if (genomefile_name == NULL)
  {
    fprintf(stderr, "No genome file specified -- exit\n");
    fclose(treefile);
    exit (EXIT_FAILURE);
  }
  if ((genomefile = fopen(genomefile_name, "r")) == NULL)
  {
    fprintf(stderr, "Failed to open genome file \"%s\" -- exit\n", genomefile_name);
    fclose(treefile);
    exit (EXIT_FAILURE);
  }
  if (corrfile_name)
  {
    if ((corrfile = fopen(corrfile_name, "w")) == NULL)
    {
      fprintf(stderr, "Failed to open \"%s\" for correlation output -- exit\n", corrfile_name);
      fclose(treefile);
      fclose(genomefile);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    corrfile = stdout;
    corrfile_name = "stdout";
  }
  if (mutfile_name)
  {
    if ((mutfile = fopen(mutfile_name, "w")) == NULL)
    {
      fprintf(stderr, "Failed to open \"%s\" for mutation output -- exit\n", mutfile_name);
      fclose(treefile);
      fclose(genomefile);
      if (corrfile != stdout)
        fclose(corrfile);
      exit (EXIT_FAILURE);
    }
  }
  else
    mutfile = NULL;
  while (!feof(genomefile) && !ferror(genomefile) && !feof(treefile) && !ferror(treefile))
  {
    fgets(buf, 256, treefile);
    if ((buf[0] != 'g') || (buf[1] != ' '))
    {
      fprintf(stderr, "Treefile corrupt after generation %ld\n", tree_generation);
      break;
    }
    tree_generation = strtol(buf + 2, NULL, 10);
    fgets(buf, 256, genomefile);
    if ((buf[0] != 'g') || (buf[1] != ' '))
    {
      fprintf(stderr, "Genome file corrupt after generation %ld\n", genome_generation);
      break;
    }
    genome_generation = strtol(buf + 2, NULL, 10);
    fgets(buf, 256, treefile);
    num_trees = strtol(buf, (char **) NULL, 10);
    if (num_trees != 1)
    {
      fprintf(stderr, "%ld trees in generation %ld -- skipping\n", num_trees, tree_generation);
      for (i = 0; i < num_trees; i++)
      {
        phyl_read_tree(treefile, &phyltree);
        phyl_free_tree(&phyltree);
      }
    }
    else
      phyl_read_tree(treefile, &phyltree);
    read_genomes(genomefile, &population);
    if ((num_trees == 1) && (tree_generation == genome_generation))
    {
      if ((tree_dist = (double *) malloc((population.size * population.size - 1) / 2 * sizeof(double))) == NULL)
        fprintf(stderr, "Failed to allocate array of tree distances\n");
      else
      {
        if ((reconst_dist = (double *) malloc(population.size * (population.size - 1) / 2 * sizeof(double))) == NULL)
          fprintf(stderr, "Failed to allocate array of reconstructed distances\n");
        else
        {
          if (outfile_name)
          {
            sprintf(buf, "%s-%05ld.gpd", outfile_name, tree_generation);
            outfile = fopen(buf, "w");
	    fprintf(outfile, "# generation %ld\n", tree_generation);
          }
          num_data = 0;
          for (i = 0; i < population.size; i++)
          {
            for (j = 0; j < i; j++)
            {
              tree_d = phyl_leafdistance(&phyltree, population.seq[i].name, population.seq[j].name);
              if (tree_d < 0)
              {
                fprintf(stderr, "Generation %ld: Trees and genomes inconsistent: error #%f\n", tree_generation, tree_d);
                break;
              }
              rec_d = diffchar_distance(population.seq[i].genome.length,
                      population.seq[i].genome.g,
                      population.seq[j].genome.length, population.seq[j].genome.g);
              change_rate = rec_d / (population.seq[i].genome.length);
              rec_d1 = 1.0 - 256.0 / 255.0 * change_rate;
              if (rec_d1 > 0.0)
              {
                rec_d = log(rec_d1);
                tree_dist[num_data] = tree_d;
                reconst_dist[num_data] = rec_d;
                num_data++;
                if (outfile)
                  fprintf(outfile, "%f %f\n", tree_d, rec_d);
              }
              else if (outfile)
                fprintf(outfile, "# Infinite distance (log operand: %f)\n", rec_d1);
            }
            if (j < i)
              break;
          }
          if (outfile_name && outfile)
            fclose(outfile);
          if ((ret_code = linear_regression(num_data, tree_dist, reconst_dist, &lambda, &y0, &r)) < 0)
          {
            fprintf(stderr, "Error #%d in linear regression\n", ret_code);
          }
          else
          {
            mutation_rate = 1.0 - exp(lambda);
            if (corrfile)
              fprintf(corrfile, "%ld %f\n", tree_generation, r);
            if (mutfile)
              fprintf(mutfile, "%ld %f\n", tree_generation, mutation_rate);
          }
          free(reconst_dist);
        }
        free(tree_dist);
      }
    }
    phyl_free_tree(&phyltree);
    free_population(&population);
#ifdef MEMDEBUG
    print_MemdebugStatistics();
#endif
  }
  if (mutfile)
    fclose(mutfile);
  if (corrfile != stdout)
    fclose(corrfile);
  fclose(treefile);
  fclose(genomefile);
  return (EXIT_SUCCESS);
}

