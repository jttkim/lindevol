/* $Id: lndgnups.c,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: lndgnups.c,v $
 * Revision 1.2  2000/01/30 03:11:00  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __atarist__
#  include <gfaport.h>
#else
#  include "gfaport.h"
#endif

#define MAX_SLEN 256

#define print_plain_command(f, bname, ext, text, epsfname) \
  fprintf(f, "set output \"%s.eps\"\n", epsfname); \
  fprintf(f, "plot \'%s.%s\' title \"%s: %s\" with lines\n", bname, ext, bname, text)

#define print_mma_command(f, bname, ext1, ext2, ext3, text, epsfname) \
  fprintf(f, "set output \"%s.eps\"\n", epsfname); \
  fprintf(f, "plot \'%s.%s\' title \"%s: min. %s\" with lines, \'%s.%s\' title \"%s: max. %s\" with lines, \'%s.%s\' title \"%s: avg. %s\" with lines\n", bname, ext1, bname, text, bname, ext2, bname, text, bname, ext3, bname, text)

#define print_twin_command(f, bname, ext1, ext2, text1, text2, epsfname) \
  fprintf(f, "set output \"%s.eps\"\n", epsfname); \
  fprintf(f, "plot \'%s.%s\' title \"%s: %s\" with lines, \'%s.%s\' title \"%s: %s\" with lines\n", bname, ext1, bname, text1, bname, ext2, bname, text2)


char bname[MAX_SLEN], buf[MAX_SLEN];


int main(int argc, char **argv)
{
  FILE *outfile;

  if (argc > 1)
  {
    strncpy(bname, argv[1], MAX_SLEN);
    bname[MAX_SLEN - 1] = '\0';
  }
  else
  {
    input_str("basename of lnd3v04 run: ", bname, MAX_SLEN);
  }
  sprintf(buf, "%s.lep", bname);
  if ((outfile = fopen(buf, "w")) == NULL)
  {
    fprintf(stderr, "failed to open output file %s -- exiting\n", buf);
    return (EXIT_FAILURE);
  }
  fprintf(outfile, "set terminal postscript eps solid 6\n");
  print_plain_command(outfile, bname, "001", "population size", "popsize");
  print_mma_command(outfile, bname, "002", "003", "004", "# of cells", "numcells");
  print_mma_command(outfile, bname, "005", "006", "007", "total energy", "tenergy");
  print_plain_command(outfile, bname, "008", "# of different genomes", "numgenom");
  fprintf(outfile, "set output \"nummnspc.eps\"\n");
  fprintf(outfile, "plot \'%s.009\' title \"%s: copies of most frequent genome\" with lines, \'%s.010\' title \"%s: copies of 2nd most frequent genome\" with lines, \'%s.011\' title \"%s: copies of 3rd most frequent genome\" with lines\n", bname, bname, bname, bname, bname, bname);
  print_plain_command(outfile, bname, "012", "editdist to prev. main species", "editdist");
  print_mma_command(outfile, bname, "013", "014", "015", "genome length", "gnmlen");
  print_mma_command(outfile, bname, "016", "017", "018", "# of used genes", "usedgens");
  print_plain_command(outfile, bname, "019", "distance distribution entropy", "distent");
  print_plain_command(outfile, bname, "020", "rel. distance distribution entropy", "rdistent");
  fprintf(outfile, "set output \"seeds.eps\"\n");
  fprintf(outfile, "plot \'%s.021\' title \"%s: # of seeds\" with lines, \'%s.022\' title \"%s: # of flying seeds\" with lines, \'%s.023\' title \"%s: # of local seeds\" with lines\n", bname, bname, bname, bname, bname, bname);
  print_plain_command(outfile, bname, "024", "# of new plants", "newplnts");
  print_twin_command(outfile, bname, "025", "026", "# of divisions", "# of new cells", "newcells");
  print_plain_command(outfile, bname, "027", "# of attacks", "attacks");
  print_plain_command(outfile, bname, "028", "# of deaths", "deaths");
  print_mma_command(outfile, bname, "029", "030", "031", "age", "age");
  print_twin_command(outfile, bname, "032", "033", "# of mut- actions", "# of mut+ actions", "mut");
  print_mma_command(outfile, bname, "034", "035", "036", "# mutation operations", "nummut");
  print_plain_command(outfile, bname, "037", "# of unmutated genomes", "unmutatd");
  print_plain_command(outfile, bname, "038", "genetic diversity", "gdivers");
  print_plain_command(outfile, bname, "039", "rel. genetic diversity", "rgdivers");

  fclose(outfile);
  return (EXIT_SUCCESS);
}

