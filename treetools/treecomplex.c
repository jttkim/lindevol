#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef __atarist__
#  include <unistd.h>
#else
#  include <getopt.h>
#endif

#include "jklib.h"
#include "ptlib.h"


#define DIST_RESOLUTION 1000L

int main(int argc, char **argv)
{
  FILE *treefile, *outfile;
  char *treefile_name = NULL, *outfile_name = NULL;
  char buf[256];
  int ret_code;
  PHYLTREE phyltree;
  long generation, num_trees, i, j, k;
  long hsize = DIST_RESOLUTION, *hgram, d;
  double cplx, treeheight;
  extern char *optarg;
  int oc;

  phyl_init_tree(&phyltree);
  while ((oc = getopt(argc, argv, "t:o:r:h")) != -1)
  {
    switch(oc)
    {
    case 't':
      treefile_name = optarg;
      break;
    case 'o':
      outfile_name = optarg;
      break;
    case 'r':
      hsize = strtol(optarg, (char **) NULL, 10);
      if (hsize < 2)
      {
	fprintf(stderr, "Illegal resolution %s -- using default %ld\n", optarg, DIST_RESOLUTION);
	hsize = DIST_RESOLUTION;
      }
      break;
    case 'h':
      printf("treecomplex -- compute distance distribution complexity\n");
      printf("    values of collections of trees\n");
      printf("\n");
      printf("Command line parameters:\n");
      printf("-t <filename>: Specify tree file name\n");
      printf("-o <filename>: Specify output file name\n");
      printf("-r <num>: Specify resolution of distance distribution\n");
      printf("-h: Print this help and exit\n");
      exit (EXIT_SUCCESS);
    }
  }
  if ((hgram = (long *) malloc(hsize * sizeof(long))) == NULL)
  {
    fprintf(stderr, "Failed to allocate distribution array -- exit\n");
    exit (EXIT_FAILURE);
  }
  if (treefile_name)
  {
    if ((treefile = fopen(treefile_name, "r")) == NULL)
    {
      fprintf(stderr, "Failed to open tree file \"%s\" -- exit\n", treefile_name);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    treefile = stdin;
    treefile_name = "stdin";
  }
  if (outfile_name)
  {
    if ((outfile = fopen(outfile_name, "w")) == NULL)
    {
      fprintf(stderr, "Failed to open output file \"%s\" -- exit\n", outfile_name);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    outfile = stdout;
    outfile_name = "stdout";
  }
  while (!feof(treefile) && !ferror(treefile))
  {
    do
      fgets(buf, 256, treefile);
    while (((buf[0] == '\0') || (buf[0] == '\n') || (buf[0] == '#')) && !feof(treefile) && !ferror(treefile));
    if (feof(treefile) || ferror(treefile))
      break;
    if ((buf[0] != 'g') || !isspace(buf[1]))
    {
      fprintf(stderr, "Treefile corrupt after generation #%ld\n", generation);
      break;
    }
    generation = strtod(buf + 2, NULL);
    do
      fgets(buf, 256, treefile);
    while (((buf[0] == '\0') || (buf[0] == '\n') || (buf[0] == '#')) && !feof(treefile) && !ferror(treefile));
    num_trees = strtod(buf, NULL);
    for (i = 0; i < num_trees; i++)
    {
      if ((ret_code = phyl_read_tree(treefile, &phyltree)) < 0)
      {
        fprintf(stderr, "Error #%d reading tree #%ld of generation %ld\n", ret_code, i, generation);
        continue;
      }
      if (!phyltree.lengthinfo_complete)
      {
        fprintf(stderr, "generation %ld, tree #%ld: incomplete tree length info -- ignored\n", generation, i);
        continue;
      }
      for (j = 0; j < hsize; j++)
	hgram[j] = 0;
      treeheight = phyl_treeheight(phyltree.root);
      for (j = 1; j < phyltree.num_leaves; j++)
      {
        for (k = 0; k < j; k++)
        {
          d = phyl_leafdistance(&phyltree, phyltree.leaf[j]->name, phyltree.leaf[k]->name) * hsize / treeheight * 0.5;
          if (d >= hsize)
            d = hsize - 1;
          hgram[d]++;
        }
      }
      cplx = shannon_long(hsize, hgram);
      fprintf(outfile, "%ld %f\n", generation, cplx);
      phyl_free_tree(&phyltree);
    }
  }
  return (EXIT_SUCCESS);
}

