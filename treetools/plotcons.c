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


#include "jklib.h"
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

  while ((oc = getopt(argc, argv, "hui:o:t:")) != -1)
  {
    switch (oc)
    {
    case 'i':
      infile_name = optarg;
      break;
    case 'o':
      outfile_name = optarg;
      break;
    case 't':
      treefile_name = optarg;
      break;
    case 'u':
      unrooted_style = 1;
      break;
    case 'h':
      printf("plotcons -- programs to plot consense trees\n");
      printf("Commandline options:\n");
      printf("-t <filename>: Specify treefile\n");
      printf("-i <filename>: Specify consense output file\n");
      printf("-o <filename>: Specify file for postscript output\n");
      printf("-u: Plot trees in \"unrooted style\"\n");
      printf("    default: vertical dendrogram style\n");
      printf("-h: Print this help and exit\n");
      exit (EXIT_SUCCESS);
    }
  }
  if (infile_name)
  {
    if ((infile = fopen(infile_name, "r")) == NULL)
    {
      fprintf(stderr, "failed to open %s for input -- exit\n", infile_name);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    fprintf(stderr, "no consense file specified -- exit\n");
    exit (EXIT_FAILURE);
  }
  if (outfile_name)
  {
    if ((outfile = fopen(outfile_name, "w")) == NULL)
    {
      fprintf(stderr, "failed to open %s for output -- exit\n", outfile_name);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    outfile_name = "stdout";
    outfile = stdout;
  }
  if (treefile_name)
  {
    if ((treefile = fopen(treefile_name, "r")) == NULL)
    {
      fprintf(stderr, "Failed to open treefile %s -- exit\n", treefile_name);
      exit (EXIT_FAILURE);
    }
  }
  else
  {
    infile_name = "stdin";
    infile = stdin;
  }
  phyl_init_tree(&phyltree);
  postscript_init(outfile);
  page_no = 1;
  if (phyl_read_tree(treefile, &phyltree) == 0)
  {
    if (get_namelist(infile) == 0)
    {
      if (prepare_leafsets(&phyltree) == 0)
      {
        write_pageheader(outfile, treefile_name, page_no++, 2000, 3200, 200);
        if (phyl_alloc_set(&phyltree, &consense_set) == 0)
        {
          if ((nameset = (char **) malloc((num_otus + 1) * sizeof(char *))) != NULL)
          {
            fgets(buf, 1024, infile);
            while (buf[0] != '\n')
            {
              printf(buf);
              j = 0;
              k = 0;
              set_complete = 1;
              for (i = 0; j < num_otus; i++)
              {
                if (isspace(buf[i]))
                  continue;
                if ((buf[i] == '\n') || (buf[i] == '\0'))
                {
                  fprintf(stderr, "error: incomplete list of species\n");
                  set_complete = 0;
                  break;
                }
                if (buf[i] == '*')
                  nameset[k++] = name[j++];
                else
                  j++;
              }
              nameset[k] = (char *) NULL;
              if (set_complete)
              {
                if ((ret_code = phyl_strtoset(&phyltree, nameset, &consense_set)) < 0)
                {
                  fprintf(stderr, "error %d converting name list to leafset\n", ret_code);
                }
                else
                {
                  percent = strtod(buf + (i + 1), (char **) NULL);
                  for (i = 0; i < num_otus; i++)
                  {
                    if (!phyl_cmp_leafset(&consense_set, set + i))
                    {
                      set[i].node->edge_counter = percent;
                      break;
                    }
                  }
                  if (i == num_otus)
                  {
                    fprintf(stderr, "error: consense set not found in tree\n");
                  }
                }
              }
              fgets(buf, 1024, infile);
            }
            free(nameset);
	    lengthtocounter(&phyltree);
            phyl_psinfo.fontname = "Courier";
            phyl_psinfo.min_fontsize = 30.0;
            phyl_psinfo.max_fontsize = 60.0;
            phyl_psinfo.label_fontsize = 60.0;
            phyl_psinfo.label = "units";
            phyl_psinfo.tic_length = 20.0;
            phyl_psinfo.label_start = 0.0;
            phyl_psinfo.linewidth = 1.0;
            phyl_psinfo.print_leafnames = 1;
            phyl_psinfo.attrlist = NULL;
            phyl_counter2thick(&phyltree, page_width / phyltree.num_leaves * 0.005);
            if (unrooted_style)
              phyl_ps_utree(outfile, &phyltree, 100.0, 100.0, 2000.0, 3000.0, &phyl_psinfo);
            else
              phyl_pstree(outfile, &phyltree, 100.0, 100.0, 2000.0, 3000.0, &phyl_psinfo);
            postscript_showpage(outfile);
            phyl_free_tree(&phyltree);
          }
          else
          {
            fprintf(stderr, "out of memory: failed to allocate temporary name array\n");
          }
        }
        else
        {
          fprintf(stderr, "failed to allocate test leaf set\n");
        }
      }
      else
      {
        fprintf(stderr, "out of memory -- failed to allocate leaf sets\n");
      }
      free_names();
    }
    else
    {
      fprintf(stderr, "failed to read name list from outfile \"%s\" -- exit\n", infile_name);
      if (treefile != stdin)
        fclose(treefile);
      fclose(infile);
      if (outfile != stdout)
        fclose (outfile);
      return (EXIT_FAILURE);
    }
  }
  if (treefile != stdin)
    fclose(treefile);
  fclose(infile);
  if (outfile != stdout)
    fclose (outfile);
  return (EXIT_SUCCESS);
}

