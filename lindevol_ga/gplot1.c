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

#include "gnlib.h"
#include "jklib.h"

#include "lndio.h"


#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif

#include "gplot.h"


#define DO_NEWSIM         1
#define DO_RESUMESIM      2
#define DO_EXIT           3
#define DO_TESTSIM        4


#include "lndglobl.h"

#ifndef PI
#  define PI 3.1415927
#endif


long max_numgenes(void)
{
  long plant_no, m = 0;

  for (plant_no = 0; plant_no < psize; plant_no++)
  {
    if (lnd_genome[plant_no].genome.num_genes > m)
      m = lnd_genome[plant_no].genome.num_genes;
  }
  return (m);
}


long num_usedgenes(const LND_GENOME *lnd_genome)
{
  long n = 0, gene_no;

  for (gene_no = 0; gene_no < lnd_genome->genome.num_genes; gene_no++)
  {
    if (lnd_genome->genome.usg_count[gene_no])
      n++;
  }
  return (n);
}


long max_num_usedgenes(void)
{
  long plant_no, n, m = 0;

  for (plant_no = 0; plant_no < psize; plant_no++)
  {
    n = num_usedgenes(lnd_genome + plant_no);
    if (n > m)
      m = n;
  }
  return (m);
}


long gene2statespec(const GENOME *genome, long pos, STATE_SPEC *statespec)
{
  if ((pos < 0) || (pos >= genome->length))
    return (-1);
  statespec->valid_bits = 0xff;
  statespec->state_bits = genome->g[pos];
  return (0);
}


int calculate_statespecs(const unsigned char *gene, STATE_SPEC *statespec)
{
  int i, j, n, state[3][3], x, y, x1, y1, x2, y2;

  state[1][1] = 1;
  for (i = 0; i < 8; i++)
  {
    if ((*gene & (1U << i)))
      state[x_offset[i] + 1][y_offset[i] + 1] = 1;
    else
      state[x_offset[i] + 1][y_offset[i] + 1] = 0;
  }
  x = x_offset[gene[1] & 7];
  y = y_offset[gene[1] & 7];
  if (state[x + 1][y + 1])  /* Nonsense gene: target pos already occupied, no cell production, no activation */
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
      if (state[x1 + 1][y1 + 1])
      {
        statespec[n].valid_bits = 0;
        statespec[n].state_bits = 0;
        for (j = 0; j < 8; j++)
        {
          x2 = x1 + x_offset[j];
          y2 = y1 + y_offset[j];
          if ((x2 >= -1) && (x2 <= 1) && (y2 >= -1) && (y2 <= 1))
          {
            statespec[n].valid_bits |= (1U << j);
            if (state[x2 + 1][y2 + 1])
              statespec[n].state_bits |= (1U << j);
          }
        }
        n++;
      }
    }
  }
  return (n);
}


int calculate_statespecs_outside(const unsigned char *gene, STATE_SPEC *statespec)
{
  int i, j, n, state[3][3], x, y, x1, y1, x2, y2;

  state[1][1] = 1;
  for (i = 0; i < 8; i++)
  {
    if ((*gene & (1U << i)))
      state[x_offset[i] + 1][y_offset[i] + 1] = 1;
    else
      state[x_offset[i] + 1][y_offset[i] + 1] = 0;
  }
  x = x_offset[gene[1] & 7];
  y = y_offset[gene[1] & 7];
  if (state[x + 1][y + 1])  /* Nonsense gene: target pos already occupied, no cell production, no activation */
    return (0);
  state[x + 1][y + 1] = 1;
  statespec[0].valid_bits = 0;
  statespec[0].state_bits = 0;
  n = 0;
  for (i = 0; i < 8; i++)
  {
    x1 = x + x_offset[i];
    y1 = y + y_offset[i];
    if ((x1 < -1) || (x1 > 1) || (y1 < -1) || (y1 > 1))
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


static char *sprint_geneno(char *buf, long gene_no, const GENOME *genome, unsigned long flags)
{
  flags &= genome->flags;

  if ((flags & GNM_USG) && (flags & GNM_BP))
    sprintf(buf, "%5ld [%5ld, %6ld]", gene_no, genome->usg_count[gene_no], genome->bp_count[gene_no]);
  else if (flags & GNM_USG)
    sprintf(buf, "%5ld [%5ld]", gene_no, genome->usg_count[gene_no]);
  else if (flags & GNM_BP)
    sprintf(buf, "%5ld [BP=%6ld]", gene_no, genome->bp_count[gene_no]);
  else
    sprintf(buf, "%5ld", gene_no);
  return (buf);
}


static void sprint_statespec(char *buf, const STATE_SPEC *statespec)
{
  int i;

  for (i = 0; i < 8; i++)
  {
    if (statespec->valid_bits & (1 << i))
    {
      if (statespec->state_bits & (1 << i))
	buf[7 - i] = '1';
      else
	buf[7 - i] = '0';
    }
    else
      buf[7 - i] = '*';
  }
  buf[8] = '\0';
}


void list_genome(FILE *f, long plant_no, int used_only)
{
  long gpos, ga, glpos, i;
  unsigned char output;
  char str[200], gnum[200], astr[200], gline[1000];
  STATE_SPEC statespec;
  const GENOME *genome = &(lnd_genome[gi[plant_no]].genome);
  unsigned long flags;

  if (bp_ceiling)
    flags = 0xffffffff;
  else
    flags = 0xffffffff ^ GNM_BP;
  gpos = 0;
  for (i = 0; i < genome->num_genes; i++)
  {
    gpos = i * GENE_LENGTH;
    sprint_geneno(gnum, i, genome, flags);
    gene2statespec(genome, gpos, &statespec);
    sprint_statespec(str, &statespec);
    output = genome->g[GENE_LENGTH * i + 1];
    ga = gene_activity(output);
    switch (ga)
    {
    case LND_DIVIDE:
      sprintf(astr, "divide %u", (output & 0x7));
      break;
/*
    case LND_STATEBIT:
      sprintf(astr, "set bit #%u", (output & 0xf));
      break;
    case LND_FLYINGSEED:
      sprintf(astr, "flying seed");
      break;
    case LND_LOCALSEED:
      sprintf(astr, "local seed");
      break;
*/
    case LND_MUTMINUS:
      sprintf(astr, "mut-");
      break;
    case LND_MUTPLUS:
      sprintf(astr, "mut+");
      break;
    default:
      sprintf(astr, "***** unknown activity code %ld *****", ga);
      break;
    }
    sprintf(gline, "%s: %s -> %s", gnum, str, astr);
    for (glpos = strlen(gline); glpos < 51; glpos++)
      gline[glpos] = ' ';
    gline[glpos] = '\0';
    sprintf(gline + strlen(gline), " %02x %02x", genome->g[gpos], genome->g[gpos + 1]);
    if (!used_only || ((genome->flags & GNM_USG) && (genome->usg_count[i])))
      fprintf(f, "%s\n", gline);
  }
  fprintf(f, "\n");
}


void list_genome_verbose(FILE *f, long plant_no, int used_only)
{
  long i, j, n;
  unsigned char gene_output;
  char input_conf[3][6] = {"     ", "  X  ", "     "}, output_conf[3][6] = {"     ", "  x  ", "     "};
  char str[200];
  long activity;
  STATE_SPEC current_statespec, statespec[9];

  fprintf(f, "raw listing:\n");
  for (i = 0; i < lnd_genome[gi[plant_no]].genome.num_genes; i++)
  {
    if (!used_only || lnd_genome[gi[plant_no]].genome.usg_count[i])
    {
      for (j = 7; j >= 0; j--)
      {
        if ((lnd_genome[gi[plant_no]].genome.g[GENE_LENGTH * i] & (1 << j)))
        {
          input_conf[1 - y_offset[j]][2 + x_offset[j] * 2] = '*';
          output_conf[1 - y_offset[j]][2 + x_offset[j] * 2] = '*';
        }
        else
        {
          input_conf[1 - y_offset[j]][2 + x_offset[j] * 2] = ' ';
          output_conf[1 - y_offset[j]][2 + x_offset[j] * 2] = ' ';
        }
      }
      gene_output = lnd_genome[gi[plant_no]].genome.g[GENE_LENGTH * i + 1];
      activity = gene_activity(gene_output);
      *str = '\0';
      if (activity == LND_DIVIDE)
      {
        j = gene_output & 0x07;
        if (output_conf[1 - y_offset[j]][2 + x_offset[j] * 2] == '*')
          output_conf[1 - y_offset[j]][2 + x_offset[j] * 2] = 'o';
        else
          output_conf[1 - y_offset[j]][2 + x_offset[j] * 2] = 'x';
      }
      else if (activity == LND_FLYINGSEED)
        sprintf(str, "flying seed  ");
      else if (activity == LND_LOCALSEED)
        sprintf(str, "local seed   ");
      else if (activity == LND_MUTMINUS)
        sprintf(str, "mut-         ");
      else if (activity == LND_MUTPLUS)
        sprintf(str, "mut+         ");
      if (lnd_genome[gi[plant_no]].genome.bp_count)
        sprintf(str + strlen(str), "[%4ld / %6ld]", lnd_genome[gi[plant_no]].genome.usg_count[i], lnd_genome[gi[plant_no]].genome.bp_count[i]);
      else
        sprintf(str + strlen(str), "[%4ld / ------]", lnd_genome[gi[plant_no]].genome.usg_count[i]);
      if (activity == LND_DIVIDE)
      {
        fprintf(f, "      %s    %s\n", input_conf[0], output_conf[0]);
        fprintf(f, "%4ld: %s -> %s        %s\n", i, input_conf[1], output_conf[1], str);
        fprintf(f, "      %s    %s\n\n", input_conf[2], output_conf[2]);
      }
      else
      {
        fprintf(f, "      %s\n", input_conf[0]);
        fprintf(f, "%4ld: %s -> %s\n", i, input_conf[1], str);
        fprintf(f, "      %s\n", input_conf[2]);
      }
      if (activity == LND_DIVIDE)
      {
        n = calculate_statespecs(lnd_genome[gi[plant_no]].genome.g + GENE_LENGTH * i, statespec);
/*
        fprintf(f, "gene: %02x %02x: ", lnd_genome[gi[plant_no]].genome.g[GENE_LENGTH * i],  lnd_genome[gi[plant_no]].genome.g[GENE_LENGTH * i + 1]);
        for (j = 0; j < n; j++)
          fprintf(f, "(%02x %02x) ", statespec[j].valid_bits, statespec[j].state_bits);
        fprintf(f, "\n");
*/
        fprintf(f, "      activator of:");
        for (j = 0; j < lnd_genome[gi[plant_no]].genome.num_genes; j++)
        {
          if (!used_only || lnd_genome[gi[plant_no]].genome.usg_count[j])
          {
            gene2statespec(&(lnd_genome[gi[plant_no]].genome), j * GENE_LENGTH, &current_statespec);
            if (statespec_match(&current_statespec, n, statespec))
              fprintf(f, " %ld,", j);
          }
        }
        fprintf(f, "\n");
      }
      fprintf(f, "\n");
    }
  }
  fprintf(f, "\n");
}


static int postscript_counterbox(FILE *f, const GENOME *genome, long gene_no, double x0, double y0, double width, double height, unsigned long flags)
{
  char buf[256], ps_str[520], *s;
  int  len;
  double fh;

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
  len = strlen(buf);
  fh = height * 0.8;
  if (fh > (width / (len + 1) / 0.7))
    fh = width / (len + 1) / 0.7;
  ps_string(buf, ps_str);
  fprintf(f, "CFB %f scalefont setfont %f %s stringwidth pop .5 mul sub %f moveto %s show\n", fh, x0 + 0.5 * width, ps_str, y0 + 0.3 * fh, ps_str); 
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
  fprintf(f, "gsave %f setlinewidth newpath %f %f moveto %f 0 rlineto %f %f rmoveto %f %f rlineto %f %f rlineto stroke grestore\n",
	  0.2 * c, cx0 + 3.5 * c, cy0 + 1.5 * c, 2 * c, -0.5 * c, 0.5 * c, 0.5 * c, -0.5 * c, -0.5 * c, -0.5 * c);
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


int postscript_connectgraph(FILE *f, const LND_GENOME *lnd_genome, int used_only,
        double x0, double y0, double width, double height, double gbox_width, double gbox_height, unsigned long counter_flags)
{
  long gene_no, j, n, n_out, gnum;
  double gene_x0, arr_x0, arrbase_offset, d, d_abs, m, h, gw;
  double spline_x0, spline_y0, spline_x1, spline_y1, spline_x2;
  double counter_height = gbox_height * 0.3;
  STATE_SPEC current_statespec, statespec[9], statespec_outside[9];
  long *gene_pos = NULL;
  long activity;

  if ((gene_pos = (long *) malloc(lnd_genome->genome.num_genes * sizeof(long))) == NULL)
    return (-1);
  fprintf(f, "save gsave\n");
  gene_x0 = x0;
  gnum = 0;
  if (counter_flags)
  {
    counter_height = gbox_height * 0.3;
    counter_flags &= lnd_genome->genome.flags;
    fprintf(f, "/CFB /Courier-Bold findfont def\n");
  }
  for (gene_no = 0; gene_no < lnd_genome->genome.num_genes; gene_no++)
  {
    if (!used_only || lnd_genome->genome.usg_count[gene_no])
    {
      if (counter_flags)
      {
        postscript_counterbox(f, &(lnd_genome->genome), gene_no, gene_x0, y0, gbox_width, counter_height, counter_flags);
        postscript_genebox(f, &(lnd_genome->genome), gene_no, gene_x0, y0 + counter_height, gbox_width, gbox_height - counter_height);
      }
      else
        postscript_genebox(f, &(lnd_genome->genome), gene_no, gene_x0, y0, gbox_width, gbox_height);
      gene_pos[gene_no] = gnum++;
      gene_x0 += gbox_width;
    }
  }
  if ((gbox_width * gnum) < width)
    gw = gbox_width * gnum;
  else
    gw = width;
/*
  printf("width = %f, gw = %f\n", width, gw);
  fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath stroke\n", x0, y0, x0 + width, y0, x0 + width, y0 + height, x0, y0 + height);
*/
  gene_x0 = x0;
  for (gene_no = 0; gene_no < lnd_genome->genome.num_genes; gene_no++)
  {
    if (!used_only || lnd_genome->genome.usg_count[gene_no])
    {
      activity = gene_activity(lnd_genome->genome.g[GENE_LENGTH * gene_no + 1]);
      if (activity == LND_DIVIDE)
      {
        n = calculate_statespecs(lnd_genome->genome.g + GENE_LENGTH * gene_no, statespec);
        n_out = calculate_statespecs_outside(lnd_genome->genome.g + GENE_LENGTH * gene_no, statespec_outside);
        for (j = 0; j < lnd_genome->genome.num_genes; j++)
        {
          if (!used_only || lnd_genome->genome.usg_count[j])
          {
            gene2statespec(&(lnd_genome->genome), j * GENE_LENGTH, &current_statespec);
            if (!statespec_match(&current_statespec, n, statespec) && statespec_match(&current_statespec, n_out, statespec_outside))
            {
	      fprintf(f, "0.5 setgray\n");
              d = (gene_pos[j] - gene_pos[gene_no]) * gbox_width;
              d += 0.4 * gbox_width;
              d_abs = (d < 0.0) ? -d : d;
              m = 0.1 * gbox_width * 0.4 * gbox_width / d / 0.2 / gbox_height;
              h = (height - gbox_height * 1.2) / (gw - 0.6 * gbox_width);
              arr_x0 = gene_x0 + 0.3 * gbox_width;
              arrbase_offset = 0.15 * gbox_height;
              if (arrbase_offset > (0.08 * gbox_width))
                arrbase_offset = 0.08 * gbox_width;
              spline_x0 = arr_x0 + m * 0.2 * gbox_height;
              spline_x1 = spline_x0 + d_abs * h * m;
              spline_x2 = arr_x0 + d - d_abs * h * m;
              if (((d > 0.0) && (spline_x1 > spline_x2)) || ((d < 0.0) && (spline_x1 < spline_x2)))
              {
                spline_x1 = (spline_x1 + spline_x2) * 0.5;
                spline_x2 = spline_x1;
              }
              spline_y0 = y0 + gbox_height * 1.2;
              spline_y1 = spline_y0 + d_abs * h;
              spline_y1 = y0 + gbox_height + (spline_x1 - arr_x0) / m;
              fprintf(f, "newpath %f %f moveto %f %f lineto\n", arr_x0, y0 + gbox_height, spline_x0, spline_y0);
#ifdef NOSPLINES
              fprintf(f, "%f %f lineto %f %f lineto %f %f lineto stroke\n", spline_x1, spline_y1,
                      spline_x2, spline_y1,
                      arr_x0 + d, spline_y0);
#else
              fprintf(f, "%f %f %f %f %f %f curveto stroke\n", spline_x1, spline_y1,
                      spline_x2, spline_y1,
                      arr_x0 + d, spline_y0);
#endif
              fprintf(f, "0 setgray newpath %f %f moveto %f %f lineto %f %f lineto closepath gsave fill grestore stroke\n",
                      arr_x0 + d - arrbase_offset, spline_y0,
                      arr_x0 + d, y0 + gbox_height,
                      arr_x0 + d + arrbase_offset, spline_y0);
            }
          }
        }
      }
      gene_x0 += gbox_width;
    }
  }
/*
 * Pass for connections to neighborhoods with mother cells inside local neighborhoods
 */
  gene_x0 = x0;
  for (gene_no = 0; gene_no < lnd_genome->genome.num_genes; gene_no++)
  {
    if (!used_only || lnd_genome->genome.usg_count[gene_no])
    {
      activity = gene_activity(lnd_genome->genome.g[GENE_LENGTH * gene_no + 1]);
      if (activity == LND_DIVIDE)
      {
        n = calculate_statespecs(lnd_genome->genome.g + GENE_LENGTH * gene_no, statespec);
        for (j = 0; j < lnd_genome->genome.num_genes; j++)
        {
          if (!used_only || lnd_genome->genome.usg_count[j])
          {
            gene2statespec(&(lnd_genome->genome), j * GENE_LENGTH, &current_statespec);
            if (statespec_match(&current_statespec, n, statespec))
            {
              d = (gene_pos[j] - gene_pos[gene_no]) * gbox_width;
              d += 0.4 * gbox_width;
              d_abs = (d < 0.0) ? -d : d;
              m = 0.1 * gbox_width * 0.4 * gbox_width / d / 0.2 / gbox_height;
              h = (height - gbox_height * 1.2) / (gw - 0.6 * gbox_width);
              arr_x0 = gene_x0 + 0.3 * gbox_width;
              arrbase_offset = 0.15 * gbox_height;
              if (arrbase_offset > (0.08 * gbox_width))
                arrbase_offset = 0.08 * gbox_width;
              spline_x0 = arr_x0 + m * 0.2 * gbox_height;
              spline_x1 = spline_x0 + d_abs * h * m;
              spline_x2 = arr_x0 + d - d_abs * h * m;
              if (((d > 0.0) && (spline_x1 > spline_x2)) || ((d < 0.0) && (spline_x1 < spline_x2)))
              {
                spline_x1 = (spline_x1 + spline_x2) * 0.5;
                spline_x2 = spline_x1;
              }
              spline_y0 = y0 + gbox_height * 1.2;
              spline_y1 = spline_y0 + d_abs * h;
              spline_y1 = y0 + gbox_height + (spline_x1 - arr_x0) / m;
              fprintf(f, "newpath %f %f moveto %f %f lineto\n", arr_x0, y0 + gbox_height, spline_x0, spline_y0);
#ifdef NOSPLINES
              fprintf(f, "%f %f lineto %f %f lineto %f %f lineto stroke\n", spline_x1, spline_y1,
                      spline_x2, spline_y1,
                      arr_x0 + d, spline_y0);
#else
              fprintf(f, "%f %f %f %f %f %f curveto stroke\n", spline_x1, spline_y1,
                      spline_x2, spline_y1,
                      arr_x0 + d, spline_y0);
#endif
              fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto closepath gsave fill grestore stroke\n",
                      arr_x0 + d - arrbase_offset, spline_y0,
                      arr_x0 + d, y0 + gbox_height,
                      arr_x0 + d + arrbase_offset, spline_y0);
            }
          }
        }
      }
      gene_x0 += gbox_width;
    }
  }
  fprintf(f, "grestore restore\n");
  return (0);
}

