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
#include "jklib.h"
#include "genomelib.h"


typedef struct
{
  char       name[100];
  GN_NODE_ID node_id;
  GENOME     genome;
} SEQ;


char buf[256];
SEQ *seq;
PHYLTREE phyltree;
long psize;


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


int read_population(FILE *f, long generation)
{
  long i, j;

  get_line(buf, 256, f);
  psize = strtol(buf, (char **) NULL, 10);
  if ((seq = (SEQ *) malloc(psize * sizeof(SEQ))) == NULL)
    return (-1);
  for (i = 0; i < psize; i++)
  {
    if (gn_read_id(&(seq[i].node_id), f) < 0)
    {
      fprintf(stderr, "error reading genome ID %ld, generation %ld -- exit\n", i, generation);
      for (j = 0; j < i; j++)
        free_genome(&(seq[i].genome));
      free(seq);
      return (-1);
    }
    gn_node_idstring(&(seq[i].node_id), seq[i].name);
    if (read_genome(f, &(seq[i].genome), 0) < 0)
    {
      fprintf(stderr, "error reading genome %ld, generation %ld -- exit\n", i, generation);
      for (j = 0; j < i; j++)
        free_genome(&(seq[i].genome));
      free(seq);
      return (-1);
    }
  }
  return (0);
}


void free_population(void)
{
  long i;

  for (i = 0; i < psize; i++)
    free_genome(&(seq[i].genome));
  free(seq);
}


int distcorr(FILE *outfile, long generation)
{
  long i, j, editdist;
  double treedist;

  for (i = 0; i < psize; i++)
  {
    for (j = 0; j < i; j++)
    {
      treedist = phyl_leafdistance(&phyltree, seq[i].name, seq[j].name);
      if (treedist < 0.0)
      {
        fprintf(stderr, "distcorr: error #%f determining tree distance between %s and %s\n", treedist, seq[i].name, seq[j].name);
        continue;
      }
      editdist = edit_distance(seq[i].genome.length, seq[i].genome.g,
              seq[j].genome.length, seq[j].genome.g);
      if (editdist < 0)
      {
        fprintf(stderr, "distcorr: error #%ld determining edit distance between %s and %s\n", editdist, seq[i].name, seq[j].name);
        continue;
      }
      fprintf(outfile, "%f %ld\n", treedist, editdist);
    }
  }
  for (i = 0; i < psize; i++)
    free_genome(&(seq[i].genome));
  free(seq);
  phyl_free_tree(&phyltree);
  return (0);
}


int main(int argc, char **argv)
{
  extern int   optind, opterr;
  extern char *optarg;
  int          optchar;
  char *treefile_name = NULL, *genomefile_name = NULL, *outfile_basename = NULL, g_outfile_name[FILENAME_MAX];
  FILE *treefile, *genomefile, *outfile;
  long generation, g, num_trees, i;
  int ret_code;

  genomefile_name = NULL;
  treefile_name = NULL;
  outfile_basename = NULL;
  generation = -1;
  phyl_init_tree(&phyltree);
  while ((optchar = getopt(argc, argv, "t:g:o:h")) != -1)
  {
    switch (optchar)
    {
    case 'g':
      genomefile_name = optarg;
      break;
    case 't':
      treefile_name = optarg;
      break;
    case 'o':
      outfile_basename = optarg;
      break;
    case 'h':
      printf("distcorr -- a hack to visualize correlation between\n");
      printf("    phylogenetic and edit distances\n");
      printf("Usage of command line options:\n");
      printf("-g <filename>: Specify genome file\n");
      printf("-t <filename>: Specify tree file\n");
      printf("-o <outfile>: specify output file basename:");
      printf("-h: Print this info and exit\n");
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
    fprintf(stderr, "no tree file specified -- exit\n");
    exit (EXIT_FAILURE);
  }
  if ((genomefile = fopen(genomefile_name, "r")) == NULL)
  {
    fprintf(stderr, "failed to open genome file \"%s\" -- exit\n", genomefile_name);
    exit (EXIT_FAILURE);
  }
  while (!feof(treefile) && !ferror(treefile) && !feof(genomefile) && !ferror(genomefile))
  {
    get_line(buf, 256, treefile);
    if (feof(treefile) || ferror(treefile))
      break;
    if (buf[0] != 'g')
    {
      fprintf(stderr, "error in header of treefile after generation %ld -- exit\n", generation);
      exit (EXIT_FAILURE);
    }
    generation = strtol(buf + 1, (char **) NULL, 10);
    get_line(buf, 256, genomefile);
    if (feof(genomefile) || ferror(genomefile))
      break;
    if (buf[0] != 'g')
    {
      fprintf(stderr, "error in header of genomefile after generation %ld -- exit\n", g);
      exit (EXIT_FAILURE);
    }
    g = strtol(buf + 1, (char **) NULL, 10);
    if (g != generation)
    {
      fprintf(stderr, "treefile: g=%ld, genomefile: g=%ld -- out of sync\n", generation, g);
      break;
    }
    get_line(buf, 256, treefile);
    num_trees = strtol(buf, NULL, 10);
    if (num_trees < 1)
    {
      fprintf(stderr, "%ld trees in treefile at generation %ld\n", num_trees, generation);
      break;
    }
    if (num_trees > 1)
    {
      fprintf(stderr, "distcorr: multiple trees in treefile\n");
      for (i = 0; i < num_trees; i++)
      {
	if ((ret_code = phyl_read_tree(treefile, &phyltree)) < 0)
	  fprintf(stderr, "error #%d while reading tree #%ld of generation %ld\n", ret_code, i, generation);
	phyl_free_tree(&phyltree);
      }
      read_population(genomefile, g);
      free_population();
      continue;
    }
    if ((ret_code = phyl_read_tree(treefile, &phyltree)) < 0)
    {
      fprintf(stderr, "error #%d while reading tree #%ld of generation %ld\n", ret_code, i, generation);
      break;
    }
    if ((ret_code = read_population(genomefile, g)) < 0)
    {
      fprintf(stderr, "error %d reading genomes of generation %ld\n", ret_code, g);
      phyl_free_tree(&phyltree);
      break;
    }
    if (outfile_basename)
    {
      sprintf(g_outfile_name, "%s-%06ld-dc.gpd", outfile_basename, generation);
      if ((outfile = fopen(g_outfile_name, "w")) == NULL)
      {
        fprintf(stderr, "failed to open \"%s\" for output -- exit\n", g_outfile_name);
        exit (EXIT_FAILURE);
      }
    }
    else
      outfile = stdout;
    ret_code = distcorr(outfile, generation);
    if (ret_code < 0)
      fprintf(stderr, "error %d in distcorr\n", ret_code);
    if (outfile_basename)
      fclose(outfile);
    free_population();
    phyl_free_tree(&phyltree);
  }
  fclose(treefile);
  fclose(genomefile);
  return (EXIT_SUCCESS);
}

