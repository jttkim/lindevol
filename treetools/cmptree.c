#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "pttypes.h"
#include "ptlib.h"


void getline(FILE *f, char *buf, int l)
{
  long i;

  buf[0] = '\0';
  while (!strlen(buf))
  {
    fgets(buf, l, f);
    for (i = strlen(buf) - 1; i >= 0; i--)
    {
      if ((buf[i] == '\n') || (buf[i] == '\r'))
        buf[i] = '\0';
    }
  }
}


int main(int argc, char **argv)
{
  FILE *f1, *f2;
  PHYLTREE tree1, tree2;
  long g1, g2, ntr1, ntr2, i;
  char buf[80];

  tree1.root = NULL;
  tree1.num_leaves = 0;
  tree1.leaf = NULL;
  tree2.root = NULL;
  tree2.num_leaves = 0;
  tree2.leaf = NULL;
  if (argc > 2)
  {
    if ((f1 = fopen(argv[1], "r")) == NULL)
    {
      fprintf(stderr, "failed to open %s -- exit.\n", argv[1]);
      exit (EXIT_FAILURE);
    }
    if ((f2 = fopen(argv[2], "r")) == NULL)
    {
      fprintf(stderr, "failed to open %s -- exit.\n", argv[2]);
      exit (EXIT_FAILURE);
    }
    while ((!feof(f1)) && (!feof(f2)))
    {
      getline(f1, buf, 80);
      if (buf[0] != 'g')
      {
        fprintf(stderr, "generation mark not found in %s\n", argv[1]);
	fprintf(stderr, "    %s\n", buf);
        fclose(f1);
        fclose(f2);
        exit (EXIT_FAILURE);
      }
      g1 = strtol(buf + 1, NULL, 10);
      getline(f2, buf, 80);
      if (buf[0] != 'g')
      {
        fprintf(stderr, "generation mark not found in %s\n", argv[2]);
	fprintf(stderr, "    %s\n", buf);
        fclose(f1);
        fclose(f2);
        exit (EXIT_FAILURE);
      }
      g2 = strtol(buf + 1, NULL, 10);
      getline(f1, buf, 80);
      ntr1 = strtol(buf, NULL, 10);
      getline(f2, buf, 80);
      ntr2 = strtol(buf, NULL, 10);
      if ((g1 != g2 ) || (ntr1 != ntr2))
      {
        printf("difference found: %s: g = %ld, n = %ld; %s: g = %ld, n = %ld\n", argv[1], g1, ntr1, argv[2], g2, ntr2);
        fclose(f1);
        fclose(f2);
        exit (EXIT_FAILURE);
      }
      /* printf("g = %ld, n = %ld\n", g1, ntr1); */

      for (i = 0; i < ntr1; i++)
      {
        if (phyl_read_tree(f1, &tree1) < 0)
        {
          fclose(f1);
          fclose(f2);
          fprintf(stderr, "error reading tree from %s, g = %ld, n = %ld\n", argv[1], g1, i);
          exit(EXIT_FAILURE);
        }
/*
        if (tree1.lengthinfo_complete)
          printf("complete branch length info found in %s\n", argv[1]);
        else
          printf("no branch length info found in %s\n", argv[1]);
        printf("***** tree from %s @ g = %ld, n = %ld *****\n", argv[1], g1, i);
        print_tree(&tree1);
*/
        if (phyl_read_tree(f2, &tree2) < 0)
        {
          phyl_free_tree(&tree1);
          fclose(f1);
          fclose(f2);
          fprintf(stderr, "error reading tree from %s, g = %ld, n = %ld\n", argv[2], g2, i);
          exit(EXIT_FAILURE);
        }
/*
        if (tree2.lengthinfo_complete)
          printf("complete branch length info found in %s\n", argv[2]);
        else
          printf("no branch length info found in %s\n", argv[2]);
        printf("***** tree from %s @ g = %ld, n = %ld *****\n", argv[2], g2, i);
        print_tree(&tree2);
*/
        printf("topological tree distance @ g = %ld, n = %ld: %ld\n", g1, i, phyl_topotreedist(&tree1, &tree2));
        phyl_free_tree(&tree1);
        phyl_free_tree(&tree2);
      }
    }
    fclose(f1);
    fclose(f2);
  }
  return (0);
}

