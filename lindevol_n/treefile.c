#include <stdio.h>
#include <stdlib.h>

#include "gntypes.h"
#include "gnlib.h"
#include "lndglobl.h"
#include "lnderror.h"
#include "lndlib.h"
#include "lndlibin.h"


FILE *jf_file = NULL;


int open_jf_file(const char *mode)
{
  jf_file = fopen(jf_file_name, mode);
  if (jf_file == NULL)
  {
    perror("open_jf_file failed");
    return (-1);
  }
  return (0);
}


void write_jf(long num_samples, const long *sample_index)
{
  long i;
  int ret_code;
  GN_TREE sample_tree;

  gn_init_tree(&sample_tree);
  i = gn_copy_tree(&sample_tree, &gntree);
  if (i < 0)
    fprintf(stderr, "write_jf: Error %ld returned by gn_copy_tree\n", i);
/*
  printf("*** gntree: *****\n");
  gn_print_tree(&gntree);
  printf("\n*** sample_tree: *****\n");
  gn_print_tree(&sample_tree);
  printf("\n");
*/
  i = num_samples;
  while ((i < world_width) && (sample_index[i] > -1))
  {
    ret_code = gn_node_death(&sample_tree, &(plant[sample_index[i]]->node), generation);
    if (ret_code < 0)
      fprintf(stderr, "write_jf: gn_node_death returned %d\n", ret_code);
    i++;
  }
  fprintf(jf_file, "g %ld\n", generation);
  gn_print_jftrees(&sample_tree, generation, jf_file);
/*
  fprintf(jf_file, "*******************************************************\n");
  gn_print_jftrees(&gntree, generation, jf_file);
  fprintf(jf_file, "*******************************************************\n");
*/
  fprintf(jf_file, "\n");
  gn_free_tree(&sample_tree);
}

