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
#include "gplot.h"

#include "lndio.h"

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


#define LND_MAIN
#include "lndglobl.h"

#ifndef PI
#  define PI 3.1415927
#endif


double page_x0, page_y0, page_width, page_height;

static char buf[MAX_SLEN + 1];


void set_pageglobals(int eps, int landscape)
{
  if (eps)
  {
    page_x0 = 0.0;
    page_y0 = 0.0;
    if (landscape)
    {
      page_width = 500.0;
      page_height = 354.0;
    }
    else
    {
      page_width = 354.0;
      page_height = 500.0;
    }
  }
  else
  {
    page_x0 = 30.0;
    page_y0 = 30.0;
    if (landscape)
    {
      page_width = 707.0;
      page_height = 500.0;
    }
    else
    {
      page_width = 500.0;
      page_height = 707.0;
    }
  }
}


void fprint_dnachars(FILE *f, unsigned char d)
{
  int i;
  char dna[5];

  for (i = 0; i < 4; i++)
  {
    switch (d & 0x03)
    {
    case 0:
      dna[i] = 'a';
      break;
    case 1:
      dna[i] = 'c';
      break;
    case 2:
      dna[i] = 'g';
      break;
    case 3:
      dna[i] = 't';
      break;
    }
    d >>= 2;
  }
  dna[4] = '\0';
  fprintf(f, dna);
}


void dna_dump(FILE *f, long plant_no)
{
  long i;

  gn_node_idstring(&(lnd_genome[plant_no].node_id), buf);
  fprintf(f, ">DL;%s\n", buf);
  fprintf(f, "run %s, at generation %ld, plant #%ld, length %ld\n", simname, start_generation, plant_no, lnd_genome[gi[plant_no]].genome.length);
  for (i = 0; i < lnd_genome[plant_no].genome.length; i++)
  {
    fprint_dnachars(f, lnd_genome[plant_no].genome.g[i]);
    if ((i % 16) == 15)
      fprintf(f, "\n");
  }
  if ((i % 16) != 0)
    fprintf(f, "\n");
  fprintf(f, "*\n");
}


int plant_boundingbox(const PLANT *plant, long *xmin, long *xmax, long *ymax)
{
  long i, dx, ww2 = world_width / 2;

  *xmin = 0;
  *xmax = 0;
  *ymax = 0;
  for (i = 0; i < plant->num_cells; i++)
  {
    dx = (plant->cell[i].x - plant->cell[0].x + world_width) % world_width;
    if (dx >= ww2)
      dx -= world_width;
    if (dx < *xmin)
      *xmin = dx;
    if (dx > *xmax)
      *xmax = dx;
    if (plant->cell[i].y > *ymax)
      *ymax = plant->cell[i].y;
  }
  return (0);
}


void print_plant(FILE *f, long plant_no, int used_only)
{
  fprintf(f, "******** plant #%1ld ******************************\n\n", plant_no);
  list_genome(f, plant_no, used_only);
}


int postscript_plant(FILE *f, const PLANT *plant, long offset, double x0, double y0, double u)
{
  long i, x, y, ww2 = world_width / 2;

  fprintf(f, "%% offset = %ld, x0 = %f, y0 = %f, u = %f\n", offset, x0, y0, u);
  fprintf(f, "gsave %f setlinewidth\n", u * 0.2);
  for (i = 0; i < plant->num_cells; i++)
  {
    x = (plant->cell[i].x - plant->cell[0].x + offset + world_width) % world_width;
    if (x >= ww2)
      x -= world_width;
    y = plant->cell[i].y;
    fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath\n",
            x0 + x * u, y0 + y * u,
            x0 + x * u + u, y0 + y * u,
            x0 + x * u + u, y0 + y * u + u,
            x0 + x * u, y0 + y * u + u);
    if (plant->cell[i].energy)
      fprintf(f, "gsave 0.2 setgray fill grestore\n");
    else
      fprintf(f, "gsave 1.0 setgray fill grestore\n");
    fprintf(f, "stroke\n");
  }
  fprintf(f, "grestore\n");
  return (ferror(f));
}


int postscript_worldheader(FILE *f, int color, int num_steps)
{
  int i;
  double r, g, b, xr, xg, xb, yr, yg, yb, radius, x, y, angle;

  fprintf(f, "%%%%Creator: stplot2 of %s\n", SIMPRGNAME);
  fprintf(f, "%%%%Title: %s-%ld\n", simname, generation);
  fprintf(f, "%%%%EndComments\n");
  fprintf(f, "%%%%BeginProlog\n");
  if (color)
  {
    fprintf(f, "/lndcolor0 { 1 0 0 setrgbcolor } bind def\n");
    fprintf(f, "/lndcolor0x { 0.5 0 0 setrgbcolor } bind def\n");
    fprintf(f, "/lndcolor1 { 0 1 0 setrgbcolor } bind def\n");
    fprintf(f, "/lndcolor1x { 0 0.5 0 setrgbcolor } bind def\n");
    fprintf(f, "/lndcolor2 { 0 0 1 setrgbcolor } bind def\n");
    fprintf(f, "/lndcolor2x { 0 0 0.5 setrgbcolor } bind def\n");
    fprintf(f, "/lndcolor3 { 1 1 0 setrgbcolor } bind def\n");
    fprintf(f, "/lndcolor3x { 0.5 0.5 0 setrgbcolor } bind def\n");
    fprintf(f, "/lndcolor4 { 0 1 1 setrgbcolor } bind def\n");
    fprintf(f, "/lndcolor4x { 0 0.5 0.5 setrgbcolor } bind def\n");
    fprintf(f, "/lndcolor5 { 1 0 1 setrgbcolor } bind def\n");
    fprintf(f, "/lndcolor5x { 0.5 0 0.5 setrgbcolor } bind def\n");
    if (num_steps > 6)
    {
      num_steps -= 6;
      radius = sqrt(0.5) * tan(0.1666667 * PI);
      xr = 1.0 / sqrt(2.0);
      xg = -xr;
      xb = 0.0;
      yr = -1.0 / sqrt(6.0);
      yg = yr;
      yb = -2.0 * yr;
      angle = 0;
      for (i = 0; i < num_steps; i++)
      {
        angle = ((double) i) / num_steps * PI * 2.0;
        x = radius * cos(angle);
        y = radius * sin(angle);
        r = x * xr + y * yr + 0.333333;
        g = x * xg + y * yg + 0.333333;
        b = x * xb + y * yb + 0.333333;
        fprintf(f, "/lndcolor%1d { %f %f %f setrgbcolor } bind def\n", i + 6, r, g, b);
        fprintf(f, "/lndcolor%1dx { %f %f %f setrgbcolor } bind def\n", i + 6, r * 0.5, g * 0.5, b * 0.5);
      }
    }
  }
  else
  {
    for (i = 0; i < num_steps; i++)
    {
      fprintf(f, "lndcolor%1d { %f setgray } bind def\n", i, 0.7 + 0.3 * i / num_steps);
      fprintf(f, "lndcolor%1dx { %f setgray } bind def\n", i, 0.3 * i / num_steps);
    }
  }
  fprintf(f, "%%%%EndProlog\n");
  return (0);
}


int postscript_world(FILE *f, double x0, double y0, double width, double height, int steps)
{
  long plant_no, cell_no, xpos, color_step;
  double x, y, u = ((width / world_width) < (height / world_height)) ? width / world_width : height / world_height;

  color_step = 0;
  fprintf(f, "%f setlinewidth\n", u * 0.2);
  for (xpos = 0; xpos < world_width; xpos++)
  {
    if (world[xpos][0].plant_no > -1)
    {
      if (world[xpos][0].cell_no == 0)
      {
        plant_no = world[xpos][0].plant_no;
	fprintf(f, "%% plant #%ld, has %ld cells\n", plant_no, plant[plant_no].num_cells);
        for (cell_no = 0; cell_no < plant[plant_no].num_cells; cell_no++)
        {
          if (plant[plant_no].cell[cell_no].energy)
            fprintf(f, "lndcolor%ldx\n", color_step);
          else
            fprintf(f, "lndcolor%ld\n", color_step);
          x = x0 + plant[plant_no].cell[cell_no].x * u;
          y = y0 + plant[plant_no].cell[cell_no].y * u;
          fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath\n", x, y, x + u, y, x + u, y + u, x, y + u);
          fprintf(f, "gsave 0 setgray stroke grestore\n");
          fprintf(f, "fill\n");
        }
        color_step = (color_step + 1) % steps;
      }
    }
  }
  return (0);
}


int postscript_singleplant(FILE *f, long plant_no, int used_only, double x0, double y0, double width, double pbox_height, double genome_height)
{
  double u, gbw, gbh, x_off;
  long xmin, xmax, ymax, n;

  plant_boundingbox(plant + plant_no, &xmin, &xmax, &ymax);
  u = pbox_height / world_height;
  if (u > width / (xmax - xmin + 1))
    u = width / (xmax - xmin + 1);
  postscript_plant(f, plant + plant_no, -xmin, x0 + (width - u * (xmax - xmin)) * 0.5, y0 + genome_height, u);
  if (used_only)
    n = num_usedgenes(lnd_genome + gi[plant_no]);
  else
    n = lnd_genome[gi[plant_no]].genome.num_genes;
  gbh = genome_height * 0.4;
  if (lnd_genome[gi[plant_no]].genome.flags & GNM_USG)
    gbw = (gbh * 0.7) * 4.0;
  else
    gbw = gbh * 4.0;
  if ((gbw * n) > width)
    gbw = width / n;
/*
  if (width / n < 72.0)
    gbw = width / n;
  else
    gbw = 72.0;
  if (lnd_genome[gi[plant_no]].genome.flags & GNM_USG)
    gbh = (gbw - 10.0) / 11.0 * 5.0;
  else
    gbh = gbw / 11.0 * 5.0;
  if (gbh > genome_height * 0.2)
    gbh = genome_height * 0.2;
  if (gbh > gbw)
    gbh = gbw;
*/
  if ((gbw * n) < width)
    x_off = 0.5 * (width - gbw * n);
  else
    x_off = 0.0;
  postscript_connectgraph(f, lnd_genome + gi[plant_no], used_only, x0 + x_off, y0, n * gbw, genome_height, gbw, gbh, GNM_USG);
  return (0);
}


int write_worldfile(FILE *f, int eps)
{
  double worldbox_width, worldbox_height, u;
  char ps_str[MAX_SLEN * 2 + 3];

  set_pageglobals(eps, 1);
  u = page_width / world_width;
  if (page_height / world_height < u)
  {
    u = page_height / world_height;
    worldbox_width = world_width * u;
    worldbox_height = page_height;
  }
  else
  {
    worldbox_height = world_height * u;
    worldbox_width = page_width;
  }
  if (eps)
  {
    fprintf(f, "%%!PS-Adobe-3.0 EPSF-3.0\n");
    fprintf(f, "%%%%BoundingBox: 0 0 %f %f\n", worldbox_width, worldbox_height);
  }
  else
    fprintf(f, "%%!\n");
  postscript_worldheader(f, 1, 6);
  if (!eps)
  {
    fprintf(f, "90 rotate 0 %f translate\n", -page_y0 - page_height);
    fprintf(f, "/Courier-Bold findfont 12 scalefont setfont\n");
    sprintf(buf, "run %s, generation %ld", simname, start_generation);
    fprintf(f, "%f %f moveto %s show\n", page_x0, page_y0 + worldbox_height + 8.0, ps_string(buf, ps_str));
  }
/*
  fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath stroke\n",
	  page_x0, page_y0,
	  page_x0 + worldbox_width, page_y0,
	  page_x0 + worldbox_width, page_y0 + worldbox_height,
	  page_x0, page_y0 + worldbox_height);
*/
  postscript_world(f, page_x0, page_y0, worldbox_width, worldbox_height, 6);
  if (eps)
    fprintf(f, "%%%%EOF\n");
  else
    fprintf(f, "showpage\n");
  return (0);
}


int write_singleplant(FILE *f, long select_plantno, int eps, int used_only)
{
  double plantbox_height, genomebox_height;
  long i;

  for (i = 0; i < psize; i++)
  {
    if (gi[i] == select_plantno)
      break;
  }
  set_pageglobals(eps, 1);
  plantbox_height = page_width * 0.5;
  if (plantbox_height > (page_height * 0.5))
    plantbox_height = page_height * 0.5;
  genomebox_height = page_height - plantbox_height;
  if (eps)
  {
    fprintf(f, "%%!PS-Adobe-3.0 EPSF-3.0\n");
    fprintf(f, "%%%%BoundingBox: 0 0 %f %f\n", page_width, page_height);
    fprintf(f, "%%%%Creator: stplot of %s\n", SIMPRGNAME);
    fprintf(f, "%%%%Title: %s-%ld\n", simname, generation);
    fprintf(f, "%%%%EndComments\n");
  }
  else
  {
    fprintf(f, "%%!\n");
    fprintf(f, "90 rotate 0 %f translate\n", -page_y0 - page_height);
    fprintf(f, "/Courier findfont 10 scalefont setfont\n");
    fprintf(f, "%f %f moveto (plant #%3ld, gi=%3ld, gi_rev=%3ld) show\n", page_x0, page_height - 10, select_plantno, gi[select_plantno], i);
  }
/*
  fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath stroke\n",
	  page_x0, page_y0,
	  page_x0 + page_width, page_y0,
	  page_x0 + page_width, page_y0 + page_height,
	  page_x0, page_y0 + page_height);
*/
  postscript_singleplant(f, select_plantno, used_only, page_x0, page_y0, page_width, plantbox_height, genomebox_height);
  if (eps)
    fprintf(f, "%%%%EOF\n");
  else
    fprintf(f, "showpage\n");
  return (0);
}


int main(int argc, char *argv[])
{
  char *datfile_name = NULL, *listfile_name = NULL, *worldfile_name = NULL, *genomefile_name = NULL, *connectfile_name = NULL, *plantfile_name = NULL;
  char *singlefile_name = "plant.ps";
  FILE *listfile, *f;
  char buf[FILENAME_MAX];
  long plant_no, m, n, xmax, xmin, ymax;
  long min_age = 0, select_plantno = -1, select_xpos = -1;
  double x0, p, gbw, cg_h = 80.0, u, plantbox_height;
  int used_only = 0, eps = 0;
  int oc;
  extern int optind;
  extern char *optarg;

  while ((oc = getopt(argc, argv, "w:g:c:o:p:a:n:x:s:euh")) != -1)
  {
    switch (oc)
    {
    case 'e':
      eps = 1;
      break;
    case 'w':
      worldfile_name = optarg;
      break;
    case 'g':
      genomefile_name = optarg;
      break;
    case 'c':
      connectfile_name = optarg;
      break;
    case 'o':
      listfile_name = optarg;
      break;
    case 'p':
      plantfile_name = optarg;
      break;
    case 'a':
      min_age = strtol(optarg, (char **) NULL, 10);
      if (min_age < 0)
      {
        fprintf(stderr, "Illegal minimal age \"%s\" -- ignored\n", optarg);
        min_age = 0;
      }
      break;
    case 'n':
      select_plantno = strtol(optarg, (char **) NULL, 10);
      if (select_plantno < 0)
      {
        fprintf(stderr, "Illegal plant number \"%s\" -- ignored\n", optarg);
        select_plantno = -1;
      }
      break;
    case 'x':
      select_xpos = strtol(optarg, (char **) NULL, 10);
      if (select_xpos < 0)
      {
        fprintf(stderr, "Illegal plant position \"%s\" -- ignored\n", optarg);
        select_xpos = -1;
      }
      break;
    case 's':
      singlefile_name = optarg;
      break;
    case 'u':
      used_only = 1;
      break;
    case 'h':
      printf("stplot1 -- make plots from LindEvol-1 data files\n\n");
      printf("Command line usage:\n");
      printf("-e: activate EPS mode\n");
      printf("-o <filename>: Specify output file\n");
      printf("-g <filename>: Specify genome file\n");
      printf("-p <filename>: Specify plant file\n");
      printf("-w <filename>: Specify world file\n");
      printf("-c <filename>: Specify files for genomes with connection graph\n");
      printf("-n <num>: Select individual plant by number\n");
      printf("-x <num>: Select individual plant by x position\n");
      printf("-s <filename>: Specify file for single plant output\n");
      printf("-u: Print / display only used genes\n");
      exit(EXIT_SUCCESS);
    }
  }
  if (argc > optind)
    datfile_name = argv[optind++];
  else
  {
    fprintf(stderr, "No simulation specified -- exit\n");
    exit(EXIT_FAILURE);
  }
  if (argc > optind)
    listfile_name = argv[optind++];
/*
  else
  {
    listfile = stdout;
    listfile_name = "stdout";
  }
*/
  if (load_named_savefile(datfile_name) != 0)
  {
    fprintf(stderr, "failed to read data for \"%s\n", simname);
    return (-1);
  }
/*
  for (plant_no = 0; plant_no < psize; plant_no++)
    printf("%3ld: %3ld\n", plant_no, gi[plant_no]);
*/
  calculate_activity_codes(&gsys_parameters);
  if (listfile_name)
  {
    if ((listfile = fopen(listfile_name, "w")) == NULL)
      fprintf(stderr, "Failed to open \"%s\" for output", listfile_name);
    else
    {
      for (plant_no = 0; plant_no < psize; plant_no++)
      {
	print_plant(listfile, plant_no, used_only);
      }
      fclose(listfile);
    }
  }
  if (genomefile_name)
  {
    if ((f = fopen(genomefile_name, "w")) == NULL)
      fprintf(stderr, "Failed to open genome file \"%s\"\n", genomefile_name);
    else
    {
    }
  }
  if (worldfile_name)
  {
    if ((f = fopen(worldfile_name, "w")) == NULL)
      fprintf(stderr, "Failed to open world file \"%s\"\n", worldfile_name);
    else
    {
      write_worldfile(f, eps);
      fclose(f);
    }
  }
  if (connectfile_name)
  {
    if ((f = fopen(connectfile_name, "w")) == NULL)
      fprintf(stderr, "Failed to open connect file \"%s\"\n", connectfile_name);
    else
    {
      if (used_only)
      {
        m = max_num_usedgenes();
        printf("max number of used genes = %ld\n", m);
      }
      else
      {
        m = max_numgenes();
        printf("max number of genes = %ld\n", m);
      }
      fprintf(f, "/usedconn_line\n");
      fprintf(f, "{0 setgray 0.1 setlinewidth} bind def\n");
      fprintf(f, "/unusedconn_line\n");
      fprintf(f, "{0.3 setgray 0.1 setlinewidth} bind def\n");
      fprintf(f, "/genebox_line\n");
      fprintf(f, "{0 setgray 0.1 setlinewidth} bind def\n");
      fprintf(f, "/Courier findfont 10 scalefont setfont\n");
      p = 760.0 - cg_h;
      x0 = 40.0;
      for (plant_no = 0; plant_no < psize; plant_no++)
      {
	if (used_only)
	  n = num_usedgenes(lnd_genome + gi[plant_no]);
	else
	  n = lnd_genome[gi[plant_no]].genome.length;
	if (n == 0)
	  continue;
	if (lnd_genome[gi[plant_no]].genome.flags & GNM_USG)
	  gbw = 44.0;
	else
	  gbw = 66.0;
	if (gbw > (450.0 / m))
	  gbw = 450.0 / m;
	if ((x0 + n * gbw) > 490.0)
	{
	  p -= cg_h + 20.0;
	  x0 = 40.0;
	}
	if (p < 30.0)
	{
	  fprintf(f, "showpage\n");
	  p = 760.0 - cg_h;
	}
	fprintf(f, "%f %f moveto (#%ld) show\n", x0, p + cg_h - 7.0, plant_no);
	postscript_connectgraph(f, lnd_genome + gi[plant_no], used_only, x0, p, n * gbw, cg_h - 10.0, gbw, 30.0, GNM_USG);
	x0 += (n + 1) * gbw;
      }
      if ((x0 > 40.0) || (p < 700.0))
        fprintf(f, "showpage\n");
      fclose(f);
    }
  }
  if (plantfile_name)
  {
    if (eps)
    {
      for (plant_no = 0; plant_no < psize; plant_no++)
      {
	sprintf(buf, "%s%03ld.eps", plantfile_name, plant_no);
        if ((f = fopen(buf, "w")) == NULL)
          fprintf(stderr, "Failed to open \"%s\" for single plant output\n", buf);
        else
        {
	  write_singleplant(f, plant_no, eps, used_only);
          fclose(f);
        }
      }
    }
    else
    {
      if ((f = fopen(plantfile_name, "w")) == NULL)
	fprintf(stderr, "error opening %s as plant output file\n", plantfile_name);
      else
      {
	set_pageglobals(0, 0);
	plantbox_height = (page_height - page_y0) / 6.0;
	fprintf(f, "/Courier findfont 10 scalefont setfont\n");
	p = page_height - plantbox_height;
	for (plant_no = 0; plant_no < psize; plant_no++)
	{
	  plant_boundingbox(plant + plant_no, &xmin, &xmax, &ymax);
	  u = (plantbox_height - 15.0) / world_height;
	  if (u > page_width / (xmax - xmin + 1))
	    u = page_width / (xmax - xmin + 1);
	  if (p < page_y0)
	  {
	    fprintf(f, "showpage\n");
	    p = page_height - plantbox_height;
	  }
	  fprintf(f, "%f %f moveto (plant #%ld) show\n", page_x0, p + plantbox_height - 11.0, plant_no);
	  postscript_plant(f, plant + plant_no, -xmin, page_x0, p, u);
	  p -= plantbox_height;
	}
	if (p < 760.0 - 200.0)
	  fprintf(f, "showpage\n");
      }
    }
  }
  if (singlefile_name)
  {
    if (select_plantno >= 0)
    {
      if (select_plantno < psize)
      {
        if ((f = fopen(singlefile_name, "w")) == NULL)
          fprintf(stderr, "Failed to open \"%s\" for single plant output\n", singlefile_name);
        else
        {
	  write_singleplant(f, select_plantno, eps, used_only);
          fclose(f);
        }
      }
      else
        fprintf(stderr, "Plant #%ld does not exist\n", select_plantno);
    }
    else if (select_xpos >= 0)
    {
      if ((select_xpos < world_width) && (world[select_xpos][0].plant_no > -1) && (world[select_xpos][0].cell_no == 0))
      {
        if ((f = fopen(singlefile_name, "w")) == NULL)
          fprintf(stderr, "Failed to open \"%s\" for single plant output\n", singlefile_name);
        else
        {
          select_plantno = world[select_xpos][0].plant_no;
	  write_singleplant(f, select_plantno, eps, used_only);
          fclose(f);
        }
      }
      else
        fprintf(stderr, "No plant rooted at (%ld, 0)\n", select_xpos);
    }
  }
  return (EXIT_SUCCESS);
}

