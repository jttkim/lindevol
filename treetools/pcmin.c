#include <ctype.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#ifdef __atarist__
#  include <unistd.h>
#else
#  include <getopt.h>
#endif


#include "ptlib.h"


long num_otus = 0;
char **name = NULL;
PHYL_LEAFSET *set = NULL;


int is_simfile(FILE *f)
{
  long pos;
  char buf[256], *c;

  pos = ftell(f);
  fgets(buf, 256, f);
  fseek(f, pos, SEEK_SET);
  if (buf[0] != 'g')
    return (0);
  strtol(buf + 1, &c, 10);
  return (buf != c);
}


void postscript_init(FILE *f)
{
  fprintf(f, "%%!\n%% prolog\n\n");
  fprintf(f, "/FSD { findfont exch scalefont def } bind def\n");
  fprintf(f, "\n%% end of prolog\n\n");
}


void postscript_pagesetup(FILE *f)
{
  fprintf(f, "%% beginning of page setup\n\n");
  fprintf(f, "save\n");
  fprintf(f, "72 300 div 72 300 div scale\n");
  fprintf(f, "/scratchstr 30 string def\n");
  fprintf(f, "\n%% end of page setup, beginning of page commands\n\n");
}


void postscript_showpage(FILE *f)
{
  fprintf(f, "showpage restore\n\n%% end of page\n\n");
}


void write_pageheader(FILE *f, const char *infile_name, long page_no, long page_width, long page_height, long left_margin)
{
  char buf[250], ps_str[520];

  postscript_pagesetup(f);
  sprintf(buf, "%s, page %ld", infile_name, page_no);
  fprintf(f, "gsave /Helvetica-Bold findfont 50 scalefont setfont\n");
  fprintf(f, "%ld %ld moveto %s show\n", left_margin, page_height - 50, ps_string(buf, ps_str));
  fprintf(f, "0 setlinecap 10 setlinewidth newpath %ld %ld moveto %ld %ld lineto stroke grestore\n",
          left_margin, page_height - 70, page_width, page_height - 70);
}


/*
 * Read the list of species written by consense into the outfile.
 */

int get_namelist(FILE *f)
{
  long i;
  char buf[1024], **n;

  for (i = 0; i < 5; i++)
    fgets(buf, 1024, f);
  num_otus = 0;
  if ((name = (char **) malloc(sizeof(char *))) == NULL)
    return (-1);
  name[0] = NULL;
  fgets(buf, 1024, f);
  while (buf[0] != '\n')
  {
    if ((n = realloc(name, (num_otus + 2) * sizeof(char *))) == NULL)
    {
      for (i = 0; i < num_otus; i++)
        free(name[i]);
      free(name);
      name = NULL;
      return (-1);
    }
    name = n;
    if ((name[num_otus] = malloc((strlen(buf) - 2) * sizeof(char))) == NULL)
    {
      for (i = 0; i < num_otus; i++)
        free(name[i]);
      free(name);
      name = NULL;
      return (-1);
    }
    i = 0;
    while (buf[i + 2] != '\n')
    {
      if (buf[i + 2] == ' ')
        name[num_otus][i] = '_';
      else
        name[num_otus][i] = buf[i + 2];
      i++;
    }
    name[num_otus][i] = '\0';
    printf("%ld: found name: %s\n", num_otus, name[num_otus]);
    name[++num_otus] = NULL;
    fgets(buf, 1024, f);
  }
  for (i = 0; i < 5; i++)
    fgets(buf, 1024, f);
  return (0);
}


void free_names(void)
{
  long i = 0;

  while (name[i])
    free(name[i++]);
  free(name);
  name = NULL;
  num_otus = 0;
}


int prepare_leafsets(const PHYLTREE *ptree)
{
  char *set_flags;
  long i, num_edges;

  if (ptree->num_leaves == 0)
    return (-1);
  if ((set = (PHYL_LEAFSET *) malloc(ptree->num_leaves * sizeof(PHYL_LEAFSET))) == NULL)
    return (-1);
  if ((set_flags = (char *) malloc(ptree->num_leaves * ptree->num_leaves * sizeof(char))) == NULL)
  {
    free(set);
    return (-1);
  }
  memset(set_flags, 0, ptree->num_leaves * ptree->num_leaves);
  for (i = 0; i < ptree->num_leaves; i++)
  {
    set[i].size = ptree->num_leaves;
    set[i].flag = set_flags + ptree->num_leaves * i;
  }
  num_edges = 0;
  phyl_get_leafsets(ptree->root, &num_edges, set);
  return (0);
}


void ltc(PHYLTREE_NODE *node)
{
  long i;

  i = node->length;
  if (i != node->edge_counter)
    printf("counter = %ld, length = %f\n", node->edge_counter, node->length);
  node->edge_counter = node->length;
  node->length = 0.0;
  for (i = 0; i < node->num_descendants; i++)
    ltc(node->descendant[i]);
}


void lengthtocounter(PHYLTREE *ptree)
{
  ltc(ptree->root);
  ptree->lengthinfo_complete = 0;
}


int main(int argc, char **argv)
{
  PHYLTREE phyltree;
  PHYL_LEAFSET consense_set;
  PHYL_PSINFO phyl_psinfo;
  int unrooted_style = 0;
  char *infile_name = NULL, *outfile_name = NULL, *treefile_name = NULL;
  FILE *infile, *outfile, *treefile;
  int ret_code, set_complete;
  char buf[1024];
  char **nameset;
  long page_no, i, j, k, page_width = 2000;
  double percent;

  int oc;
  extern char *optarg;

  if ((ret_code = phyl_strtoset(&phyltree, nameset, &consense_set)) < 0)
  {
    fprintf(stderr, "error %d converting name list to leafset\n", ret_code);
  }
  else
  {
    printf("OK\n");
  }
  return (EXIT_SUCCESS);
}

