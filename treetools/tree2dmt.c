#include <stdio.h>
#include <stdlib.h>

#include <ptlib.h>


int main(int argc, char **argv)
{
  char *infile_name = NULL, *outfile_name = NULL;
  FILE *infile, *outfile;
  PHYLTREE ptree;
  long i, j;
  int ret_code;
  double d;

  phyl_init_tree(&ptree);
  if (argc > 1)
    infile_name = argv[1];
  if (infile_name)
  {
    if ((infile = fopen(infile_name, "r")) == NULL)
    {
      fprintf(stderr, "Failed to open %s for input -- exit\n", infile_name);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    infile_name = "stdin";
    infile = stdin;
  }
  if (outfile_name)
  {
    if ((outfile = fopen(outfile_name, "w")) == NULL)
    {
      fprintf(stderr, "Failed to open %s for output -- exit\n", outfile_name);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    outfile_name = "stdout";
    outfile = stdout;
  }
  ret_code = phyl_read_tree(infile, &ptree);
  if (ret_code < 0)
  {
    fprintf(stderr, "Error #%d while reading tree -- exit\n", ret_code);
    phyl_free_tree(&ptree);
    exit (EXIT_FAILURE);
  }
  fprintf(outfile, "%ld\n", ptree.num_leaves);
  for (i = 0; i < ptree.num_leaves; i++)
  {
    fprintf(outfile, "%10s ", ptree.leaf[i]->name);
    for (j = 0; j < ptree.num_leaves; j++)
    {
      d = phyl_leafdistance(&ptree, ptree.leaf[i]->name, ptree.leaf[j]->name);
      fprintf(outfile, " %7.3f", d);
    }
    fprintf(outfile, "\n");
  }
  if (infile != stdin)
    fclose(infile);
  if (outfile != stdout)
    fclose(outfile);
  phyl_free_tree(&ptree);
  return (EXIT_SUCCESS);
}

