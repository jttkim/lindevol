/* $Id: stplot6.c,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: stplot6.c,v $
 * Revision 1.2  2000/01/30 03:11:00  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

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
#include "lnddispl.h"
#include "lnderror.h"
#include "lndsignl.h"
#include "lndlib.h"
#include "lnd6.h"

#include "gnlib.h"
#include "jklib.h"
#include "gplot.h"

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

  gn_node_idstring(&(plant[plant_no]->node), buf);
  fprintf(f, ">DL;%s\n", buf);
  fprintf(f, "run %s, at generation %ld, plant #%ld, length %ld\n", simname, start_generation, plant_no, plant[plant_no]->genome.length);
  for (i = 0; i < plant[plant_no]->genome.length; i++)
  {
    fprint_dnachars(f, plant[plant_no]->genome.g[i]);
    if ((i % 16) == 15)
      fprintf(f, "\n");
  }
  if ((i % 16) != 0)
    fprintf(f, "\n");
  fprintf(f, "*\n");
}


void print_plant(FILE *f, long plant_no, int used_only)
{
  long i;
  double death_p;

   death_p = death_probability(plant_no, p_random_death, rdeath_f_energy, rdeath_f_numcells, leanover_penalty);
  fprintf(f, "******** plant #%1ld ******************************\n\n", plant_no);
  fprintf(f, "genealogy node id: ");
  gn_save_id(&(plant[plant_no]->node), f);
  fprintf(f, "number of cells:   %ld\n", plant[plant_no]->num_cells);
  fprintf(f, "cellular energy:   %ld\n", plant[plant_no]->cellular_energy);
  fprintf(f, "cellular nutrient: %ld\n", plant[plant_no]->cellular_nutrient);
  fprintf(f, "energy pool:       %ld\n", plant[plant_no]->energy_pool);
  fprintf(f, "nutrient pool:     %ld\n", plant[plant_no]->nutrient_pool);
  fprintf(f, "death probability: %1.12g\n", death_p);
  fprintf(f, "age:               %ld\n", plant[plant_no]->age);
  fprintf(f, "genome length:     %ld\n", plant[plant_no]->genome.length);
  fprintf(f, "genome dump:\n");
  for (i = 0; i < plant[plant_no]->genome.length; i++)
  {
    fprintf(f, "%02x ", plant[plant_no]->genome.g[i]);
    if ((i % 25) == 24)
      fprintf(f, "\n");
  }
  fprintf(f, "\nGenome listing:\n\n");
  list_genome(f, plant_no, used_only);
}


int postscript_plant(FILE *f, const PLANT *plant, long offset, double x0, double y0, double u)
{
  long i, x, y, ww2 = world_width / 2;

  fprintf(f, "%% offset = %ld, x0 = %f, y0 = %f, u = %f\n", offset, x0, y0, u);
  fprintf(f, "gsave\n");
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
      fprintf(f, "gsave 0 setgray fill grestore\n");
    else
      fprintf(f, "gsave 0.7 setgray fill grestore\n");
    fprintf(f, "stroke\n");
  }
  fprintf(f, "grestore\n");
  return (ferror(f));
}


int postscript_worldheader(FILE *f, int color, int num_steps)
{
  int i;
  double r, g, b, xr, xg, xb, yr, yg, yb, radius, x, y, angle;

  fprintf(f, "%%%%Creator: stplot of %s\n", SIMPRGNAME);
  fprintf(f, "%%%%Title: %s-%ld\n", simname, start_generation);
  fprintf(f, "%%%%EndComments\n");
  fprintf(f, "%%%%BeginProlog\n");
  fprintf(f, "/B { newpath moveto rlineto rlineto rlineto closepath } bind def\n");
  fprintf(f, "/S { newpath moveto rlineto rlineto rlineto closepath stroke} bind def\n");
  fprintf(f, "/F { newpath moveto rlineto rlineto rlineto closepath fill} bind def\n");
  fprintf(f, "/G { 0 setgray } bind def\n");
  fprintf(f, "/SC1 { 0.6 0.3 0.3 setrgbcolor } bind def\n");
  fprintf(f, "/SC2 { 0.3 0.3 0.3 setrgbcolor } bind def\n");
  fprintf(f, "/SC3 { 0.6 0.0 0.0 setrgbcolor } bind def\n");
  fprintf(f, "/SC4 { 0.0 0.0 0.0 setrgbcolor } bind def\n");
  if (color)
  {
    fprintf(f, "/LC0 { 1 0 0 setrgbcolor } bind def\n");
    fprintf(f, "/LX0 { 0.5 0 0 setrgbcolor } bind def\n");
    fprintf(f, "/LC1 { 0 1 0 setrgbcolor } bind def\n");
    fprintf(f, "/LX1 { 0 0.5 0 setrgbcolor } bind def\n");
    fprintf(f, "/LC3 { 0 0 1 setrgbcolor } bind def\n");
    fprintf(f, "/LX3 { 0 0 0.5 setrgbcolor } bind def\n");
    fprintf(f, "/LC2 { 1 1 0 setrgbcolor } bind def\n");
    fprintf(f, "/LX2 { 0.5 0.5 0 setrgbcolor } bind def\n");
    fprintf(f, "/LC5 { 0 1 1 setrgbcolor } bind def\n");
    fprintf(f, "/LX5 { 0 0.5 0.5 setrgbcolor } bind def\n");
    fprintf(f, "/LC4 { 1 0 1 setrgbcolor } bind def\n");
    fprintf(f, "/LX4 { 0.5 0 0.5 setrgbcolor } bind def\n");
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
        fprintf(f, "/LC%1d { %f %f %f setrgbcolor } bind def\n", i + 6, r, g, b);
        fprintf(f, "/LX%1d { %f %f %f setrgbcolor } bind def\n", i + 6, r * 0.5, g * 0.5, b * 0.5);
      }
    }
  }
  else
  {
    for (i = 0; i < num_steps; i++)
    {
      fprintf(f, "LC%1d { %f setgray } bind def\n", i, 0.7 + 0.3 * i / num_steps);
      fprintf(f, "LX%1d { %f setgray } bind def\n", i, 0.3 * i / num_steps);
    }
  }
  fprintf(f, "%%%%EndProlog\n");
  return (0);
}


int postscript_world(FILE *f, double x0, double y0, double width, double height, int steps)
{
  long plant_no, cell_no, xpos, ypos, color_step;
  double x, y;
  double u = ((width / world_width) < (height / world_height)) ? width / world_width : height / world_height;
  double u1 = u * 0.08, u2 = u - 2.0 * u1;

  fprintf(f, "%f setlinewidth\n", 2.0 * u1);
  for (xpos = 0; xpos < world_width; xpos++)
  {
    if (world[xpos][world_soil].plant_no > -1)
    {
      if (world[xpos][world_soil].cell_no == 0)
      {
        plant_no = world[xpos][world_soil].plant_no;
        for (cell_no = 0; cell_no < plant[plant_no]->num_cells; cell_no++)
        {
	  color_step = plant_no % steps;
          x = x0 + plant[plant_no]->cell[cell_no].x * u;
          y = y0 + plant[plant_no]->cell[cell_no].y * u;
          fprintf(f, "%1.4f 0.0 0.0 %1.4f %1.4f 0.0  %1.4f %1.4f S\n", -u, u, u, x, y);
          if (plant[plant_no]->cell[cell_no].energy)
            fprintf(f, "LX%ld ", color_step);
          else
            fprintf(f, "LC%ld ", color_step);
          fprintf(f, "%1.4f 0.0 0.0 %1.4f %1.4f 0.0  %1.4f %1.4f F G\n", -u2, u2, u2, x + u1, y + u1);
        }
      }
    }
  }
  for (xpos = 0; xpos < world_width; xpos++)
  {
    for (ypos = 0; ypos <= world_soil; ypos++)
    {
      if (world[xpos][ypos].plant_no == -1)
      {
	x = x0 + xpos * u;
	y = y0 + ypos * u;
        if (world[xpos][ypos].nutrient)
	{
	  if (world[xpos][ypos].organic_nutrient)
	    fprintf(f, "SC1 ");
	  else
	    fprintf(f, "SC2 ");
	}
	else
	{
	  if (world[xpos][ypos].organic_nutrient)
	    fprintf(f, "SC3 ");
	  else
	    fprintf(f, "SC4 ");
	}
	fprintf(f, "%1.4f 0.0 0.0 %1.4f %1.4f 0.0  %1.4f %1.4f F G\n", -u2, u2, u2, x + u1, y + u1);
      }
    }
  }
  return (0);
}


int postscript_singleplant(FILE *f, const PLANT *plant, int used_only, double x0, double y0, double width, double pbox_height, double genome_height)
{
  double u, gbw, gbh;
  long xmin, xmax, ymax, n;

  plant_boundingbox(plant, &xmin, &xmax, &ymax);
  u = pbox_height / world_height;
  if (u > width / (xmax - xmin + 1))
    u = width / (xmax - xmin + 1);
  postscript_plant(f, plant, -xmin, x0 +(width - u * (xmax - xmin)) * 0.5, y0 + genome_height, u);
  if (used_only)
    n = num_usedgenes(plant);
  else
    n = plant->genome.num_genes;
  if (width / n < 72.0)
    gbw = width / n;
  else
    gbw = 72.0;
  gbh = genome_height * 0.2;
  if (gbh > gbw)
    gbh = gbw;
  return (0);
}


int write_worldfile(FILE *f, int eps)
{
  double worldbox_width, worldbox_height, u;
  char ps_str[MAX_SLEN * 2 + 3];

  set_pageglobals(eps, 1);
  u = page_width / world_width;
  if (eps)
  {
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
  }
  else
  {
    if ((page_height - 20.0) / world_height < u)
    {
      u = (page_height - 20.0) / world_height;
      worldbox_width = world_width * u;
      worldbox_height = (page_height - 20.0);
    }
    else
    {
      worldbox_height = world_height * u;
      worldbox_width = page_width;
    }
  }
  if (eps)
  {
    fprintf(f, "%%!PS-Adobe-3.0 EPSF-3.0\n");
    fprintf(f, "%%%%BoundingBox: 0 0 %f %f\n", worldbox_width, worldbox_height);
  }
  else
  {
    fprintf(f, "%%!\n");
  }
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

  set_pageglobals(eps, 1);
  plantbox_height = page_width * 0.5;
  genomebox_height = page_height - plantbox_height;
  if (eps)
  {
    fprintf(f, "%%!PS-Adobe-3.0 EPSF-3.0\n");
    fprintf(f, "%%%%BoundingBox: 0 0 %f %f\n", page_width, page_height);
    fprintf(f, "%%%%Creator: stplot of %s\n", SIMPRGNAME);
    fprintf(f, "%%%%Title: %s-%ld\n", simname, start_generation);
    fprintf(f, "%%%%EndComments\n");
  }
  else
  {
    fprintf(f, "%%!\n");
    fprintf(f, "90 rotate 0 %f translate\n", -page_y0 - page_height);
  }
/*
  fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath stroke\n",
          page_x0, page_y0,
          page_x0 + page_width, page_y0,
          page_x0 + page_width, page_y0 + page_height,
          page_x0, page_y0 + page_height);
*/
  postscript_singleplant(f, plant[select_plantno], used_only, page_x0, page_y0, page_width, plantbox_height, genomebox_height);
  if (eps)
    fprintf(f, "%%%%EOF\n");
  else
    fprintf(f, "showpage\n");
  return (0);
}


int write_ppmfile(const char *fname, int pixel_size)
{
  FILE *f;
  int i, j;
  char rgb[6] = "     ", bgchar;
  unsigned int p;
  long x, y, b, plant_no;

  if ((f = fopen(fname, "w")) == NULL)
    return (-1);
  fprintf(f, "P3\n");
  fprintf(f, "# simulator: %s, run: %s, generation: %ld\n", SIMPRGNAME, simname, generation);
  fprintf(f, "%ld %ld\n", world_width * pixel_size, world_height * pixel_size);
  fprintf(f, "3\n");
  for (y = world_height - 1; y >= 0; y--)
  {
    for (i = 0; i < pixel_size; i++)
    {
      for (x = 0; x < world_width; x++)
      {
	plant_no = world[x][y].plant_no;
	if (plant_no > -1)
	{
	  if (world[x][y].cell_no == 0)
	    bgchar = '0';
	  else
	    bgchar = '1';
	  p = plant_no % 6 + 1;
/*
	  if (plant[plant_no]->cell[world[x][y].cell_no].energy)
	  {
	    if (p & 1U)
	      rgb[0] = '2';
	    else
	      rgb[0] = bgchar;
	    if (p & 2U)
	      rgb[2] = '2';
	    else
	      rgb[2] = bgchar;
	    if (p & 4U)
	      rgb[4] = '2';
	    else
	      rgb[4] = bgchar;
	  }
	  else
*/
	  {
	    if (p & 1U)
	      rgb[0] = '3';
	    else
	      rgb[0] = bgchar;
	    if (p & 2U)
	      rgb[2] = '3';
	    else
	      rgb[2] = bgchar;
	    if (p & 4U)
	      rgb[4] = '3';
	    else
	      rgb[4] = bgchar;
	  }
	}
	else
	{
	  if (y < world_soil)
	  {
	    if (world[x][y].nutrient)
	      sprintf(rgb, "1 1 1");
	    else
	      sprintf(rgb, "0 0 0");
	  }
	  else if (y == world_soil)
	  {
	    if (world[x][y].nutrient)
	      sprintf(rgb, "1 1 1");
	    else
	      sprintf(rgb, "3 3 3");
	  }
	  else
	    sprintf(rgb, "3 3 3");
	}
	for (j = 0; j < pixel_size; j++)
	  fprintf(f, "%s\n", rgb);
      }
    }
  }
  fclose(f);
  return (0);
}


int main(int argc, char *argv[])
{
  char *datfile_name = NULL, *listfile_name = NULL, *worldfile_name = NULL, *genomefile_name = NULL, *plantfile_name = NULL;
  char *singlefile_name = "plant.ps", *ppmfile_name = NULL;
  FILE *listfile, *f;
  long plant_no, m, n, xmax, xmin, ymax, i;
  long min_age = 0, select_plantno = -1, select_xpos = -1;
  double x0, p, gbw, cg_h = 80.0, u, plantbox_height;
  int pixel_size = 1;
  int used_only = 0, eps = 0;
  int oc;
  extern int optind;
  extern char *optarg;

  while ((oc = getopt(argc, argv, "w:g:o:p:a:n:x:s:m:z:euh")) != -1)
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
    case 'o':
      listfile_name = optarg;
      break;
    case 'p':
      plantfile_name = optarg;
      break;
    case 'm':
      ppmfile_name = optarg;
      break;
    case 'z':
      pixel_size = strtol(optarg, NULL, 10);
      if (pixel_size < 1)
	pixel_size = 1;
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
  calculate_activity_codes(&gsys_parameters);
  if (listfile_name)
  {
    if ((listfile = fopen(listfile_name, "w")) == NULL)
      fprintf(stderr, "Failed to open \"%s\" for output", listfile_name);
    else
    {
      for (plant_no = 0; plant_no < world_width; plant_no++)
      {
        if (plant[plant_no] != NULL)
        {
          print_plant(listfile, plant_no, used_only);
        }
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
      for (plant_no = 0; plant_no < world_width; plant_no++)
      {
        if (plant[plant_no] != NULL)
        {
	  fprintf(f, "%ld\n", plant[plant_no]->genome.length);
	  for (i = 0; i < plant[plant_no]->genome.length; i++)
	    fprintf(f, "%02x", plant[plant_no]->genome.g[i]);
	  fprintf(f, "\n");
        }
      }
      fclose(f);
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
  if (plantfile_name)
  {
    if ((f = fopen(plantfile_name, "w")) == NULL)
      fprintf(stderr, "error opening %s as plant output file\n", plantfile_name);
    else
    {
      set_pageglobals(0, 0);
      plantbox_height = (page_height - page_y0) / 6.0;
      fprintf(f, "/Courier findfont 10 scalefont setfont\n");
      p = page_height - plantbox_height;
      for (plant_no = 0; plant_no < world_width; plant_no++)
      {
        if (plant[plant_no])
        {
          plant_boundingbox(plant[plant_no], &xmin, &xmax, &ymax);
          u = (plantbox_height - 15.0) / world_height;
          if (u > page_width / (xmax - xmin + 1))
            u = page_width / (xmax - xmin + 1);
          if (p < page_y0)
          {
            fprintf(f, "showpage\n");
            p = page_height - plantbox_height;
          }
          fprintf(f, "%f %f moveto (plant #%ld) show\n", page_x0, p + plantbox_height - 11.0, plant_no);
          postscript_plant(f, plant[plant_no], -xmin, page_x0, p, u);
          p -= plantbox_height;
        }
      }
      if (p < 760.0 - 200.0)
        fprintf(f, "showpage\n");
    }
  }
  if (ppmfile_name)
  {
    write_ppmfile(ppmfile_name, pixel_size);
  }
  if (singlefile_name)
  {
    if (select_plantno >= 0)
    {
      if ((select_plantno < world_width) && plant[select_plantno])
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

