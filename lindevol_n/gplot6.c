#include <ctype.h>
#include <errno.h>
#include <float.h>
#ifdef __atarist__
#  include <unistd.h>
#else
#  include <getopt.h>
#endif
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __atarist__
#  include <gfaport.h>
#else
#  include "gfaport.h"
#endif

#include "genomelib.h"
#include "lndtypes.h"
#include "lnderror.h"
#include "lndsignl.h"
#include "lndlib.h"
#include "lnd6.h"

#include "gnlib.h"
#include "jklib.h"

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif

#define DO_NEWSIM         1
#define DO_RESUMESIM      2
#define DO_EXIT           3
#define DO_TESTSIM        4


#include "lndglobl.h"

#ifndef PI
#  define PI 3.1415927
#endif

#ifndef GENE_LENGTH
#  define GENE_LENGTH 2
#endif


long max_numgenes(void)
{
  long plant_no, m = 0;

  for (plant_no = 0; plant_no < world_width; plant_no++)
  {
    if ((plant[plant_no]) && (plant[plant_no]->genome.num_genes > m))
      m = plant[plant_no]->genome.num_genes;
  }
  return (m);
}


long num_usedgenes(const PLANT *plant)
{
  long n = 0, gene_no;

  for (gene_no = 0; gene_no < plant->genome.num_genes; gene_no++)
  {
    if (plant->genome.usg_count[gene_no])
      n++;
  }
  return (n);
}


long max_num_usedgenes(void)
{
  long plant_no, n, m = 0;

  for (plant_no = 0; plant_no < world_width; plant_no++)
  {
    if (plant[plant_no])
    {
      n = num_usedgenes(plant[plant_no]);
      if (n > m)
        m = n;
    }
  }
  return (m);
}


long gene2statespec(const GENOME *genome, long pos, STATE_SPEC *statespec)
{
  return (compute_statespec(genome, pos, statespec, NULL));
}


int calculate_statespecs(const unsigned char *gene, STATE_SPEC *statespec)
{
  int i, j, n, state[3][3], x, y, x1, y1, x2, y2;

  state[1][1] = 1;
  for (i = 0; i < 8; i++)
  {
    if ((*gene & (1 << i)))
      state[x_offset[i] + 1][y_offset[i] + 1] = 1;
    else
      state[x_offset[i] + 1][y_offset[i] + 1] = 0;
  }
  x = x_offset[gene[1] & 7];
  y = y_offset[gene[1] & 7];
  if (state[x + 1][y + 1])
    return (0);
  state[x + 1][y + 1] = 1;
  statespec[0].valid_bits = 0;
  statespec[0].state_bits = 0;
  for (j = 0; j < 8; j++)
  {
    x2 = x + x_offset[j];
    y2 = y + y_offset[j];
    if ((x2 >= -1) && (x2 <= 1) && (y2 >= -1) && (y2 <= 1))
    {
      statespec[0].valid_bits |= (1 << j);
      if (state[x2 + 1][y2 + 1])
        statespec[0].state_bits |= (1 << j);
    }
  }
  n = 1;
  for (i = 0; i < 8; i++)
  {
    x1 = x + x_offset[i];
    y1 = y + y_offset[i];
    if ((x1 >= -1) && (x1 <= 1) && (y1 >= -1) && (y1 <= 1))
    {
      if (state[x1][y1])
      {
        statespec[n].valid_bits = 0;
        statespec[n].state_bits = 0;
        for (j = 0; j < 8; j++)
        {
          x2 = x1 + x_offset[j];
          y2 = y1 + y_offset[j];
          if ((x2 >= -1) && (x2 <= 1) && (y2 >= -1) && (y2 <= 1))
          {
            statespec[n].valid_bits |= (1 << j);
            if (state[x2 + 1][y2 + 1])
              statespec[n].state_bits |= (1 << j);
          }
        }
        n++;
      }
    }
  }
  return (n);
}


int statespec_match(const STATE_SPEC *current_statespec, int n, const STATE_SPEC *statespec)
{
  int i;
  unsigned long valid_bits;

  for (i = 0; i < n; i++)
  {
    valid_bits = current_statespec->valid_bits & statespec[i].valid_bits;
    if (((current_statespec->state_bits ^ statespec[i].state_bits) & valid_bits) == 0)
      return (1);
  }
  return (0);
}


static char *sprint_geneno(char *buf, long gpos, long gene_no, const GENOME *genome, unsigned long flags)
{
  flags &= genome->flags;

  if ((flags & GNM_USG) && (flags & GNM_BP))
    sprintf(buf, "%5ld (p=%5ld) [%5ld, %6ld]", gene_no, gpos, genome->usg_count[gene_no], genome->bp_count[gene_no]);
  else if (flags & GNM_USG)
    sprintf(buf, "%5ld (p=%5ld) [%5ld]", gene_no, gpos, genome->usg_count[gene_no]);
  else if (flags & GNM_BP)
    sprintf(buf, "%5ld (p=%5ld) [BP=%6ld]", gene_no, gpos, genome->bp_count[gene_no]);
  else
    sprintf(buf, "%5ld (p=%5ld)", gene_no, gpos);
  return (buf);
}


static void sprint_statespec(char *buf, const STATE_SPEC *statespec)
{
  int i;

  for (i = 0; i < NUM_STATEBITS; i++)
  {
    if (statespec->valid_bits & (1 << i))
    {
      if (statespec->state_bits & (1 << i))
	buf[NUM_STATEBITS - i - 1] = '1';
      else
	buf[NUM_STATEBITS - i - 1] = '0';
    }
    else
      buf[NUM_STATEBITS - i - 1] = '*';
  }
  buf[NUM_STATEBITS] = '\0';
}


void list_genome(FILE *f, long plant_no, int used_only)
{
  long gpos, gstart, ga, glpos, i;
  unsigned char output;
  char str[200], gnum[200], astr[200], gline[1000];
  long gene_no = 0;
  STATE_SPEC statespec;
  const GENOME *genome = &(plant[plant_no]->genome);

  gene_no = 0;
  gpos = 0;
  while (gpos < genome->length)
  {
    while ((gpos < genome->length) && !l4_promoter(genome->g[gpos]))
      gpos++;
    if (gpos >= genome->length)
      break;
    gstart = gpos;
    sprint_geneno(gnum, gpos, gene_no, genome, 0xffffffff);
    gpos++;
    if (gpos >= genome->length)
      statespec.valid_bits = 0;
    else
      gpos = compute_statespec(genome, gpos, &statespec, NULL);
    sprint_statespec(str, &statespec);
    if ((gpos >= genome->length) || !l4_terminator(genome->g[gpos]))
    {
      sprintf(astr, "no action");
    }
    else
    {
      output = genome->g[gpos];
      ga = gene_activity(output);
      switch (ga)
      {
      case LND_DIVIDE:
	sprintf(astr, "divide %u", (output & 0x7));
        break;
      case LND_STATEBIT:
        sprintf(astr, "set bit #%u", (output & 0xf));
        break;
      case LND_FLYINGSEED:
        sprintf(astr, "flying seed");
        break;
      case LND_LOCALSEED:
        sprintf(astr, "local seed");
        break;
      case LND_MUTMINUS:
        sprintf(astr, "mut-");
        break;
      case LND_MUTPLUS:
        sprintf(astr, "mut+");
        break;
      case LND_TO_EPOOL:
        sprintf(astr, "to epool");
        break;
      case LND_TO_NPOOL:
        sprintf(astr, "to npool");
        break;
      case LND_FROM_EPOOL:
        sprintf(astr, "from epool");
        break;
      case LND_FROM_NPOOL:
        sprintf(astr, "from npool");
        break;
      default:
        sprintf(astr, "***** unknown activity code %ld *****", ga);
        break;
      }
    }
    sprintf(gline, "%s: %s -> %s", gnum, str, astr);
    for (glpos = strlen(gline); glpos < 43 + NUM_STATEBITS; glpos++)
      gline[glpos] = ' ';
    gline[glpos] = '\0';
    for (i = gstart; i <= gpos; i++)
    {
      if (i == genome->length)
      {
	sprintf(gline + strlen(gline), " <END>");
	break;
      }
      sprintf(gline + strlen(gline), " %02x", genome->g[i]);
    }
    if (!used_only || ((genome->flags & GNM_USG) && (genome->usg_count[gene_no])))
      fprintf(f, "%s\n", gline);
    gene_no++;
  }
  fprintf(f, "\n");
}


int postscript_counterbox(FILE *f, const GENOME *genome, long gene_no, double x0, double y0, double width, double height, double fontheight, unsigned long flags)
{
  char buf[256], ps_str[520], *s;

  fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath stroke\n",
          x0, y0,
          x0 + width, y0,
          x0 + width, y0 + height,
          x0, y0 + height);
  sprintf(buf, "%ld", gene_no);
  if ((genome->flags & GNM_USG) && (flags & GNM_USG))
  {
    s = buf + strlen(buf);
    sprintf(s, ": u=%ld", genome->usg_count[gene_no]);
  }
  if ((genome->flags & GNM_BP) && (flags & GNM_BP))
  {
    s = buf + strlen(buf);
    if ((genome->flags & GNM_USG) && (flags & GNM_USG))
      sprintf(s, ", c=%ld", genome->bp_count[gene_no]);
    else
      sprintf(s, ": c=%ld", genome->bp_count[gene_no]);
  }
    

  fprintf(f, "%f %f moveto %s show\n", x0 + 0.3 * fontheight, y0 + 0.3 * fontheight, ps_string(buf, ps_str)); 
  return (0);
}


int postscript_genebox(FILE *f, const GENOME *genome, long gene_no, double x0, double y0, double width, double height)
{
  double c, cx0, cy0;
  int i;
  long activity;
  unsigned long gene_lhs = genome->g[gene_no * GENE_LENGTH], gene_rhs = genome->g[gene_no * GENE_LENGTH + 1];

  fprintf(f, "gsave newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath stroke\n",
          x0, y0, x0 + width, y0, x0 + width, y0 + height, x0, y0 + height);
  c = width / 11.0;
  cx0 = x0 + c;
  cy0 = y0 + (height - 3.0 * c) * 0.5;
  if (c > (height / 5.0))
  {
    c = height / 5.0;
    cx0 = x0 + (width - 9.0 * c) * 0.5;
    cy0 = y0 + c;
  }
  activity = gene_activity(gene_rhs);
  fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath gsave 0 setgray fill grestore stroke\n",
          cx0 + c, cy0 + c,
          cx0 + c + c, cy0 + c,
          cx0 + c + c, cy0 + c + c,
          cx0 + c, cy0 + c + c);
  for (i = 0; i < 8; i++)
  {
    if ((gene_lhs & (1 << i)))
    {
      fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath stroke\n",
              cx0 + (x_offset[i] + 1) * c, cy0 + (y_offset[i] + 1) * c,
              cx0 + (x_offset[i] + 1) * c + c, cy0 + (y_offset[i] + 1) * c,
              cx0 + (x_offset[i] + 1) * c + c, cy0 + (y_offset[i] + 1) * c + c,
              cx0 + (x_offset[i] + 1) * c, cy0 + (y_offset[i] + 1) * c + c);
    }
  }
  switch (activity)
  {
  case LND_DIVIDE:
    fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath gsave 0 setgray fill grestore stroke\n",
            cx0 + 7 * c, cy0 + c,
            cx0 + 8 * c, cy0 + c,
            cx0 + 8 * c, cy0 + c + c,
            cx0 + 7 * c, cy0 + c + c);
    for (i = 0; i < 8; i++)
    {
      if ((gene_lhs & (1 << i)))
      {
        fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath stroke\n",
                cx0 + (x_offset[i] + 7) * c, cy0 + (y_offset[i] + 1) * c,
                cx0 + (x_offset[i] + 7) * c + c, cy0 + (y_offset[i] + 1) * c,
                cx0 + (x_offset[i] + 7) * c + c, cy0 + (y_offset[i] + 1) * c + c,
                cx0 + (x_offset[i] + 7) * c, cy0 + (y_offset[i] + 1) * c + c);
      }
    }
    fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath gsave .5 setgray fill grestore stroke\n",
            cx0 + (x_offset[gene_rhs & 0x07] + 7) * c, cy0 + (y_offset[gene_rhs & 0x07] + 1) * c,
            cx0 + (x_offset[gene_rhs & 0x07] + 7) * c + c, cy0 + (y_offset[gene_rhs & 0x07] + 1) * c,
            cx0 + (x_offset[gene_rhs & 0x07] + 7) * c + c, cy0 + (y_offset[gene_rhs & 0x07] + 1) * c + c,
            cx0 + (x_offset[gene_rhs & 0x07] + 7) * c, cy0 + (y_offset[gene_rhs & 0x07] + 1) * c + c);
    break;
  case LND_FLYINGSEED:
    fprintf(f, "newpath %f %f %f 0 360 arc stroke\n", cx0 + 7.5 * c, cy0 + 1.5 * c, c);
    break;
  case LND_LOCALSEED:
    fprintf(f, "newpath %f %f %f 0 360 arc gsave 0 setgray fill grestore stroke\n", cx0 + 7.5 * c, cy0 + 1.5 * c, c);
    break;
  case LND_MUTMINUS:
    fprintf(f, "newpath %f %f moveto %f %f lineto stroke\n", cx0 + 6.5 * c, cy0 + 1.5 * c, cx0 + 8.5 * c, cy0 + 1.5 * c);
    break;
  case LND_MUTPLUS:
    fprintf(f, "newpath %f %f moveto %f %f lineto stroke\n", cx0 + 6.5 * c, cy0 + 1.5 * c, cx0 + 8.5 * c, cy0 + 1.5 * c);
    fprintf(f, "newpath %f %f moveto %f %f lineto stroke\n", cx0 + 7.5 * c, cy0 + 0.5 * c, cx0 + 7.5 * c, cy0 + 2.5 * c);
    break;
  }
  fprintf(f, "grestore\n");
  return (0);
}

