#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "pttypes.h"
#include "ptlib.h"


int main(int argc, char **argv)
{
  FILE *f;
  long num_trees = 0, num_errors = 0;
  PHYLTREE phyltree;
  int ret_code;

  if (argc > 1)
  {
    if ((f = fopen(argv[1], "r")) == NULL)
      return (EXIT_FAILURE);
  }
  else
    f = stdin;
  phyl_init_tree(&phyltree);
  while (!feof(f) && !ferror(f))
  {
    if ((ret_code = phyl_read_tree(f, &phyltree)) == 0)
    {
      if (phyltree.num_leaves > 0)
      {
	num_trees++;
	phyl_free_tree(&phyltree);
      }
    }
    else
    {
      num_errors++;
      fprintf(stderr, "errorcode %d while trying to read tree\n", ret_code);
    }
    if (feof(f) || ferror(f))
      break;
  }
  printf("%ld\n", num_trees);
  return (EXIT_SUCCESS);
}

