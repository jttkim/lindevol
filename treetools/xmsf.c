#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __atarist__
#  include <unistd.h>
#else
#  include <getopt.h>
#endif

#include "gnlib.h"
#include "genomelib.h"


#define CODE_NONE 0
#define CODE_DNA  1
#define CODE_GC   2


typedef struct
{
  GN_NODE_ID node_id;
  GENOME     genome;
} SEQ;


char buf[256];


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


int write_generation_phylipdna(FILE *infile, FILE *outfile, long generation, int code)
{
  SEQ *seq;
  long i, j, l, psize, max_len;
  char name[80];

  get_line(buf, 256, infile);
  psize = strtol(buf, (char **) NULL, 10);
  if ((seq = (SEQ *) malloc(psize * sizeof(SEQ))) == NULL)
    return (-1);
  for (i = 0; i < psize; i++)
  {
    if (gn_read_id(&(seq[i].node_id), infile) < 0)
    {
      fprintf(stderr, "error reading genome ID %ld, generation %ld -- exit\n", i, generation);
      for (j = 0; j < i; j++)
        free_genome(&(seq[i].genome));
      free(seq);
      return (-1);
    }
    if (read_genome(infile, &(seq[i].genome), 0) < 0)
    {
      fprintf(stderr, "error reading genome %ld, generation %ld -- exit\n", i, generation);
      for (j = 0; j < i; j++)
        free_genome(&(seq[i].genome));
      free(seq);
      return (-1);
    }
  }
  max_len = 0;
  for (i = 0; i < psize; i++)
  {
    if (max_len < seq[i].genome.length)
      max_len = seq[i].genome.length;
  }
  switch (code)
  {
  case CODE_DNA:
    fprintf(outfile, "%ld %ld\n", psize, max_len * 4);
    for (i = 0; i < psize; i++)
    {
      gn_node_idstring(&(seq[i].node_id), name);
      l = strlen(name);
      for (j = 0; j < 10; j++)
      {
	if (j < l)
	  fputc(name[j], outfile);
	else
	  fputc(' ', outfile);
      }
      for (j = 0; j < ((max_len < 13) ? max_len : 13); j++)
      {
	if (j < seq[i].genome.length)
	  fprint_dnachars(outfile, seq[i].genome.g[j]);
	else
	  fprintf(outfile, "----");
      }
      fprintf(outfile, "\n");
    }
    for (l = j; l < max_len; l += 15)
    {
      fprintf(outfile, "\n");
      for (i = 0; i < psize; i++)
      {
	for (j = l; j < ((l + 15 < max_len) ? l + 15 : max_len); j++)
	{
	  if (j < seq[i].genome.length)
	    fprint_dnachars(outfile, seq[i].genome.g[j]);
	  else
	    fprintf(outfile, "----");
	}
	fprintf(outfile, "\n");
      }
    }
    break;
  case CODE_GC:
    fprintf(outfile, "%ld %ld\n", psize, max_len * 8);
    for (i = 0; i < psize; i++)
    {
      gn_node_idstring(&(seq[i].node_id), name);
      l = strlen(name);
      for (j = 0; j < 10; j++)
      {
	if (j < l)
	  fputc(name[j], outfile);
	else
	  fputc(' ', outfile);
      }
      for (j = 0; j < ((max_len < 7) ? max_len : 7); j++)
      {
	if (j < seq[i].genome.length)
	  fprint_binchars(outfile, seq[i].genome.g[j], "cg");
	else
	  fprintf(outfile, "--------");
      }
      fprintf(outfile, "\n");
    }
    for (l = j; l < max_len; l += 8)
    {
      fprintf(outfile, "\n");
      for (i = 0; i < psize; i++)
      {
	for (j = l; j < ((l + 8 < max_len) ? l + 8 : max_len); j++)
	{
	  if (j < seq[i].genome.length)
	    fprint_binchars(outfile, seq[i].genome.g[j], "cg");
	  else
	    fprintf(outfile, "--------");
	}
	fprintf(outfile, "\n");
      }
    }
    break;
  default:
    fprintf(stderr, "Code type %d unknown\n", code);
    for (i = 0; i < psize; i++)
      free_genome(&(seq[i].genome));
    free(seq);
    return (-1);
  }
  for (i = 0; i < psize; i++)
    free_genome(&(seq[i].genome));
  free(seq);
  return (0);
}


int write_generation_pirdna(FILE *infile, FILE *outfile, const char *fname, long generation, int code)
{
  GENOME genome;
  GN_NODE_ID node;
  long i, psize;
  char name[80], comment[512];

  get_line(buf, 256, infile);
  psize = strtol(buf, (char **) NULL, 10);
  for (i = 0; i < psize; i++)
  {
    if (gn_read_id(&node, infile) < 0)
    {
      fprintf(stderr, "error reading genome ID %ld, generation %ld -- exit\n", i, generation);
      exit (EXIT_FAILURE);
    }
    if (read_genome(infile, &genome, 0) < 0)
    {
      fprintf(stderr, "error reading genome %ld, generation %ld -- exit\n", i, generation);
      exit (EXIT_FAILURE);
    }
    gn_node_idstring(&node, name);
    sprintf(comment, "File: %s, generation: %ld, #%ld, length: %ld", fname, generation, i, genome.length);
    switch (code)
    {
    case CODE_DNA:
      write_pirdna(outfile, &genome, name, comment);
      break;
    case CODE_GC:
      write_pirbin(outfile, &genome, name, comment, "cg");
      break;
    default:
      fprintf(stderr, "Code type %d unknown\n", code);
      free_genome(&genome);
      return (-1);
    }
    free_genome(&genome);
  }
  return (0);
}


int main(int argc, char **argv)
{
  extern int   optind, opterr;
  extern char *optarg;
  int          optchar;
  char *infile_name = NULL, *outfile_name = NULL, *namefile_name = NULL, g_outfile_name[FILENAME_MAX];
  FILE *infile, *outfile, *namefile;
  int pir_format, phylip_format, code = CODE_DNA;
  long generation, x_generation;

  pir_format = 1;
  phylip_format = 0;
  infile_name = NULL;
  outfile_name = NULL;
  generation = -1;
  x_generation = -1;
  while ((optchar = getopt(argc, argv, "pyg:i:o:n:c:h")) != -1)
  {
    switch (optchar)
    {
    case 'p':
      pir_format = 1;
      phylip_format = 0;
      break;
    case 'y':
      pir_format = 0;
      phylip_format = 1;
      break;
    case 'g':
      x_generation = strtol(optarg, (char **) NULL, 10);
      break;
    case 'i':
      infile_name = optarg;
      break;
    case 'o':
      outfile_name = optarg;
      break;
    case 'n':
      namefile_name = optarg;
      break;
    case 'c':
      if (!strcmp(optarg, "dna"))
	code = CODE_DNA;
      else if (!strcmp(optarg, "cg") || !strcmp(optarg, "gc"))
	code = CODE_GC;
      else
	fprintf(stderr, "Unknown code specification \"%s\" -- ignored\n", optarg);
      break;
    case 'h':
      printf("xmsf -- extract multiple sequence files from simulation data\n");
      printf("Usage of command line options:\n");
      printf("-p  - create files in PIR format\n");
      printf("-y  - create files in phylip format\n");
      printf("-g <generation #>: extract from the specified generation\n");
      printf("    the next higher generation is used if the specified\n");
      printf("    generation is not found in savefile\n");
      printf("-i <infile>: specify infile, mandatory\n");
      printf("-o <outfile>: specify output file(s):");
      printf("    when extracting individual generation: specified name\n");
      printf("    is used\n");
      printf("    when extracting multipe generations: specified name is\n");
      printf("    used as base name\n");
      printf("-c {dna|cg}: Specify encoding in DNA alphabed (a, c, g, t) or\n");
      printf("    in binary alphabet (g, c).\n");
      exit (EXIT_SUCCESS);
    }
  }
  if (infile_name == NULL)
  {
    fprintf(stderr, "no input file specified -- exit\n");
    exit (EXIT_FAILURE);
  }
  if ((infile = fopen(infile_name, "r")) == NULL)
  {
    fprintf(stderr, "failed to open \"%s\" for input -- exit\n", infile_name);
    exit (EXIT_FAILURE);
  }
  if (namefile_name)
  {
    if ((namefile = fopen(namefile_name, "w")) == NULL)
    {
      fprintf(stderr, "Failed to open %s as file of filenames -- exit\n", namefile_name);
      fclose(infile);
      exit (EXIT_FAILURE);
    }
  }
  else
    namefile = NULL;
  if (x_generation == -1)
  {
    while (!feof(infile) && !ferror(infile))
    {
      get_line(buf, 256, infile);
      if (feof(infile) || ferror(infile))
        break;
      if (buf[0] != 'g')
      {
        fprintf(stderr, "error in header after generation %ld -- exit\n", generation);
        exit (EXIT_FAILURE);
      }
      generation = strtol(buf + 1, (char **) NULL, 10);
      if (outfile_name)
      {
        if (pir_format)
          sprintf(g_outfile_name, "%s-%06ld.pir", outfile_name, generation);
        else if (phylip_format)
          sprintf(g_outfile_name, "%s-%06ld.phy", outfile_name, generation);
        if ((outfile = fopen(g_outfile_name, "w")) == NULL)
        {
          fprintf(stderr, "failed to open \"%s\" for output -- exit\n", g_outfile_name);
          exit (EXIT_FAILURE);
        }
      }
      else
        outfile = stdout;
      if (pir_format)
        write_generation_pirdna(infile, outfile, infile_name, generation, code);
      else if (phylip_format)
        write_generation_phylipdna(infile, outfile, generation, code);
      if (outfile_name)
        fclose(outfile);
      if (outfile_name && namefile)
	fprintf(namefile, "%s\n", g_outfile_name);
    }
  }
  else
  {
    do
    {
      get_line(buf, 256, infile);
      if (feof(infile) || ferror(infile))
      {
        fprintf(stderr, "failed to find generation %ld in \"%s\"\n", x_generation, infile_name);
        exit (EXIT_FAILURE);
      }
      if (buf[0] == 'g')
        generation = strtol(buf + 1, (char **) NULL, 10);
    }
    while (generation < x_generation);
    if (outfile_name)
    {
      if ((outfile = fopen(outfile_name, "w")) == NULL)
      {
        fprintf(stderr, "failed to open \"%s\" for output -- exit\n", outfile_name);
        exit (EXIT_FAILURE);
      }
    }
    else
      outfile = stdout;
    if (pir_format)
      write_generation_pirdna(infile, outfile, infile_name, generation, code);
    else if (phylip_format)
      write_generation_phylipdna(infile, outfile, generation, code);
  }
  fclose(infile);
  if (namefile)
    fclose(namefile);
  return (EXIT_SUCCESS);
}

