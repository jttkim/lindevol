#include <stdio.h>
#include <stdlib.h>
#ifdef __atarist__
#  include <unistd.h>
#else
#  include <getopt.h>
#endif

#include "jklib.h"
#include "ptlib.h"


int main(int argc, char **argv)
{
  PHYLTREE tree1, tree2;
  PHYL_PSINFO phyl_psinfo;
  FILE *treefile, *infile, *outfile;
  char *treefile_name = NULL, *infile_name = NULL, *outfile_name = NULL;
  int unrooted_style = 0;
  long n_leaves;
  extern char *optarg;
  int oc;

  phyl_init_tree(&tree1);
  phyl_init_tree(&tree2);
  while ((oc = getopt(argc, argv, "t:i:o:uh")) != -1)
  {
    switch (oc)
    {
    case 't':
      treefile_name = optarg;
      break;
    case 'i':
      infile_name = optarg;
      break;
    case 'o':
      outfile_name = optarg;
      break;
    case 'u':
      unrooted_style = 1;
      break;
    case 'h':
      printf("trevstre -- plot two trees in one, using branch length\n");
      printf("    values of second tree as thickness values\n");
      exit (EXIT_SUCCESS);
    }
  }
  if (treefile_name == NULL)
  {
    fprintf(stderr, "No tree file specified -- exit\n");
    exit (EXIT_FAILURE);
  }
  if ((treefile = fopen(treefile_name, "r")) == NULL)
  {
    fprintf(stderr, "Failed to open tree fil \"%s\" -- exit\n", treefile_name);
    exit (EXIT_FAILURE);
  }
  if (infile_name)
  {
    if ((infile = fopen(infile_name, "r")) == NULL)
    {
      fprintf(stderr, "Failed to open input file \"%s\" -- exit\n", infile_name);
      fclose(treefile);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    infile = stdin;
    infile_name = "stdin";
  }
  if (outfile_name)
  {
    if ((outfile = fopen(outfile_name, "w")) == NULL)
    {
      fprintf(stderr, "Failed to open output file \"%s\" -- exit\n", outfile_name);
      fclose(treefile);
      if (infile != stdin)
        fclose(infile);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    outfile = stdout;
    outfile_name = "stdout";
  }
  if (phyl_read_tree(treefile, &tree1) < 0)
  {
    fprintf(stderr, "Error reading tree from \"%s\" -- exit\n", treefile_name);
    fclose(treefile);
    if (infile != stdin)
      fclose(infile);
    exit (EXIT_FAILURE);
  }
  fclose(treefile);
  if (phyl_read_tree(infile, &tree2) < 0)
  {
    fprintf(stderr, "Error reading tree from \"%s\" -- exit\n", infile_name);
    phyl_free_tree(&tree1);
    if (infile != stdin)
      fclose(infile);
    exit (EXIT_FAILURE);
  }
  if (infile != stdin)
    fclose(infile);
  n_leaves = phyl_num_leaves(tree1.root);
  phyl_set_thickness(&tree1, 0.0);
  if (phyl_tree2thick(&tree1, &tree2, 500.0 / n_leaves * 0.5) < 0)
  {
    fprintf(stderr, "Error setting thickness values according to tree 2\n");
  }
  phyl_psinfo.fontname = "Courier";
  phyl_psinfo.min_fontsize = 6.0;
  phyl_psinfo.max_fontsize = 12.0;
  phyl_psinfo.label_fontsize = 8.0;
  phyl_psinfo.label = "units";
  phyl_psinfo.tic_length = 20.0;
  phyl_psinfo.label_start = 0;
  phyl_psinfo.linewidth = 0.2;
  phyl_psinfo.print_leafnames = 1;
  phyl_psinfo.attrlist = NULL;
  if (unrooted_style)
    phyl_ps_utree(outfile, &tree1, 72, 100, 400, 500, &phyl_psinfo);
  else
    phyl_pstree(outfile, &tree1, 72, 100, 400, 500, &phyl_psinfo);
  fprintf(outfile, "showpage\n");
  if (outfile != stdout)
    fclose(outfile);
  phyl_free_tree(&tree1);
  phyl_free_tree(&tree2);
  return (EXIT_SUCCESS);
}

