#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "jklib.h"


int main(int argc, char *argv[])
{
  int oc;
  extern char *optarg;
  extern int optind;
  char *infile_name = NULL, *outfile_name = NULL;
  FILE *infile, *outfile;
  unsigned long i, j, n;
  long s16, k, r;
  char buf[1024];

  while ((oc = getopt(argc, argv, "h")) != -1)
  {
    switch(oc)
    {
    case 'h':
      printf("-h: print this help and exit\n");
      exit(EXIT_SUCCESS);
    }
  }
  if (optind < argc)
    infile_name = argv[optind++];
  if (optind < argc)
    outfile_name = argv[optind++];
  if (infile_name)
  {
    if ((infile = fopen(infile_name, "rb")) == NULL)
    {
      fprintf(stderr, "Failed to open \"%s\" for input -- exit\n", infile_name);
      exit(EXIT_FAILURE);
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
      fprintf(stderr, "Failed to open \"%s\" for output -- exit\n", outfile_name);
      if (infile != stdin)
        fclose(infile);
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    outfile = stdout;
    outfile_name = "stdout";
  }
  fgets(buf, 1024, infile);
  for (i = 0; !feof(infile); i++)
  {
    fgets(buf, 1024, infile);
    n = strtoul(buf, NULL, 10);
    fprintf(stderr, "n = %lu\n", n);
    for (j = 0; j < n && !feof(infile); j++)
    {
      s16 = read_int16(infile);
      if (s16 < 0)
      {
	r = -s16;
	s16 = read_int16(infile);
	for (k = 0; k < r; k++)
	  fprintf(outfile, "%lu %lu %ld\n", i, j, s16);
      }
      else
	fprintf(outfile, "%lu %lu %ld\n", i, j, s16);
    }
  }
  if (infile != stdin)
    fclose(infile);
  if (outfile != stdout)
    fclose(outfile);
  exit(EXIT_SUCCESS);
}

