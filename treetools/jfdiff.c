#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "pttypes.h"
#include "ptlib.h"


#define MAX_LINELENGTH 256


void getline(FILE *f, char *buf, int l)
{
  long i;

  buf[0] = '\0';
  while (!strlen(buf))
  {
    fgets(buf, l, f);
    if (feof(f) || ferror(f))
      break;
    for (i = strlen(buf) - 1; i >= 0; i--)
    {
      if ((buf[i] == '\n') || (buf[i] == '\r'))
        buf[i] = '\0';
    }
  }
}


long next_generation(FILE *f)
{
  char buf[MAX_LINELENGTH + 1];

  do
    getline(f, buf, MAX_LINELENGTH);
  while ((buf[0] != 'g' ) && !feof(f) && !ferror(f));
  return (strtol(buf + 1, (char **) NULL, 10));
}


int main(int argc, char **argv)
{
  FILE *f1, *f2;
  PHYLTREE tree1, tree2;
  long g1, g2, ntr1, ntr2, i, d;
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
      g1 = next_generation(f1);
      g2 = next_generation(f2);
      while (g1 < g2)
      {
        fprintf(stderr, "generation %ld is missing in %s\n", g1, argv[2]);
        g1 = next_generation(f1);
        if (feof(f1) || ferror(f1))
          break;
      }
      while (g2 < g1)
      {
        fprintf(stderr, "generation %ld is missing in %s\n", g2, argv[1]);
        g2 = next_generation(f2);
        if (feof(f2) || ferror(f2))
          break;
      }
      if (ferror(f1))
      {
        fprintf(stderr, "error reading from %s\n", argv[1]);
        break;
      }
      if (ferror(f2))
      {
        fprintf(stderr, "error reading from %s\n", argv[2]);
        break;
      }
      if (feof(f1) || feof(f2))
      {
        if (!feof(f2))
          fprintf(stderr, "generation %ld and subsequent ones missing from %s\n", g2, argv[1]);
        else if (!feof(f1))
          fprintf(stderr, "generation %ld and subsequent ones missing from %s\n", g1, argv[2]);
        break;
      }
      getline(f1, buf, 80);
      ntr1 = strtol(buf, NULL, 10);
      getline(f2, buf, 80);
      ntr2 = strtol(buf, NULL, 10);
      if (ntr1 != ntr2)
      {
        fprintf(stderr, "generation %ld: %ld trees in %s, %ld trees in %s\n", g1, ntr1, argv[1], ntr2, argv[2]);
        continue;
      }

      /* printf("g = %ld, n = %ld\n", g1, ntr1); */

      for (i = 0; i < ntr1; i++)
      {
        if (phyl_read_tree(f1, &tree1) < 0)
        {
          fprintf(stderr, "error reading tree from %s, g = %ld, n = %ld\n", argv[1], g1, i);
          continue;
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
          fprintf(stderr, "error reading tree from %s, g = %ld, n = %ld\n", argv[2], g2, i);
          continue;
        }
/*
        if (tree2.lengthinfo_complete)
          printf("complete branch length info found in %s\n", argv[2]);
        else
          printf("no branch length info found in %s\n", argv[2]);
        printf("***** tree from %s @ g = %ld, n = %ld *****\n", argv[2], g2, i);
        print_tree(&tree2);
*/
        d = phyl_topotreedist(&tree1, &tree2);
        if (d >= 0)
          printf("%ld %ld\n", g1, d);
        else
        {
          switch (d)
          {
          case PHYLERR_DIFFNUMLEAVES:
            fprintf(stderr, "different number of leaves, generation %ld, tree #%ld\n", g1, i);
            break;
          case PHYLERR_INCOMPATLEAVES:
            fprintf(stderr, "different sets of leaves: generation %ld, tree #%ld\n", g1, i);
            break;
          default:
            fprintf(stderr, "error code %ld during tree distance computation\n", d);
            break;
          }
        }
        /* printf("topological tree distance @ g = %ld, n = %ld: %ld\n", g1, i, d); */
        phyl_free_tree(&tree1);
        phyl_free_tree(&tree2);
      }
    }
    fclose(f1);
    fclose(f2);
  }
  return (0);
}

