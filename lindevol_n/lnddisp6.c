/* This module contains the functions for displaying the LindEvol world. */

#ifndef __atarist__ /* ASCII display for non-Atari systems */

#  include "lnddasc.c"

#else /* For the Atari system, a graphical, windowed user interface is used. */
      /* Use this at your own risk, this code was inherited from old versions
       * without being checked for functionality
       */

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <osbind.h>
#include <aesbind.h>
#include <vdibind.h>
#include <gemfast.h>

#include <gfaport.h>
#include <wndhndlr.h>

#ifdef MEMDEBUG
#  include <memdebug.h>
#endif

#include "lnderror.h"
#include "lndvals.h"
#include "lndtypes.h"
#include "lnd6.h"
#include "lndglobl.h"


#define SINGLEPLANT_CELLSIZE   8
#define HISTORY_BUFSIZE      200
#define HWINDOW_MAX          100
#define LANE_WIDTH           (HWINDOW_MAX + 15)
#define NUM_LANES              4

#ifdef DEBUG
#  define report(x) printf("\033Y\040\160%s, l=%1d: %s\n", __FILE__, __LINE__, x)
#else
#  define report(x)
#endif


/* struct for holding plant history in a ring buffer */

typedef struct
{
  int  bufsize;          /* size of the buffer */
  int  num_generations;  /* number of generations stored in buffer */
  char buffer_full;      /* flag indicating whether buffer is full */
  long first_generation; /* first generation stored in buffer */
  char num_cells[HISTORY_BUFSIZE];
  char energy[HISTORY_BUFSIZE];
  char num_mutations[HISTORY_BUFSIZE];
  char genome_length[HISTORY_BUFSIZE];
} HISTORY_BUFFER;

/* descriptor for a single plant window: */

typedef struct
{
  int whandle;
  char title[80];
  char info[80];
  int cell_size;
  enum {pw_plant, pw_genome, pw_history} contents;
  char single_step;
  HISTORY_BUFFER hbuf;
} PLANT_WDESCR;


PLANT_WDESCR **plant_wdescr = NULL;

char display_active = 0;
char display_shutdown = 0;

FILE  *world_file = NULL;
char   world_file_name[MAX_SLEN];

int ap_id, aes_screenhandle, v_screenhandle, dw_x, dw_y, dw_w, dw_h;
int work_in[11], work_out[57];

static int world_whandle = 0;
static char title[MAX_SLEN], info[MAX_SLEN];

static int celldisplay_size;
static int default_char_height;
static long highlit_plant = -1;
static char pause_flag = 0;

int open_world_file(const char *simname, const char *mode)
{
  sprintf(world_file_name, "%s%swld", simname, ".");
  world_file = fopen(world_file_name, mode);
  if (world_file == NULL)
  {
    perror("open_world_file failed");
    return (-1);
  }
  return (0);
}


void close_world_file(void)
{
  if (world_file != NULL)
    fclose(world_file);
}

static void alloc_plant_wdescr(void)
{
  long i;

  plant_wdescr = (PLANT_WDESCR **) malloc(world_width * sizeof(PLANT_WDESCR *));
  if (plant_wdescr != NULL)
  {
    for (i = 0; i < world_width; i++)
    {
      plant_wdescr[i] = NULL;
    }
  }
  else
  {
    form_alert(1,"[1][Not enough memory|for plant window|descriptors. Single plant|windows are disabled][ Ok ]");
  }
}


/*
 * find the plant at the screen coordinates (x/y). If no
 * plant is there, -1 is returned, otherwise the index of
 * the plant in the plant[] array.
 */

long find_plant(int x, int y)
{
  long xpos, ypos, cx, cy;
  int w_handle, x0, y0, w0, h0;

  w_handle = wind_find(x, y);
  if ((w_handle > 0) && (w_handle == world_whandle))
  {
    get_slider_positions(w_handle, &xpos, &ypos);
    wind_get(w_handle, WF_WORKXYWH, &x0, &y0, &w0, &h0);
    cx = (xpos + x - x0) / celldisplay_size;
    cy = world_height - 1 - (ypos + y - y0) / celldisplay_size;
    return (world[cx][cy].plant_no);
  }
  else
  {
    return (-1);
  }
}


/*
 * returns the number of the plant displayed in the single plant window
 * specified by w_handle. If this window is not a single plant window,
 * whandle_to_plantno returns a negative value
 */

long whandle_to_plantno(int w_handle)
{
  long plant_no;

  if (plant_wdescr == NULL)
  {
    return (-1);
  }
  for (plant_no = 0; plant_no < world_width; plant_no++)
  {
    if (plant_wdescr[plant_no] != NULL)
    {
      if (plant_wdescr[plant_no]->whandle == w_handle)
      {
        return (plant_no);
      }
    }
  }
  return (-1);
}


/*
 * find the next plant to the right (direction = 1) or to the
 * left of plant plant_no. If plant_no equals -1, the search
 * to the left starts at the right end of the world, and a search
 * to the right starts at the left side of the world.
 * The function is not intended for values other than -1 and 1 for
 * the direction parameter.
 */

long next_plant(long plant_no, int direction)
{
  long x, i;

  if (plant_no > -1)
  {
    x = (plant[plant_no]->cell[0].x + direction + world_width) % world_width;
  }
  else
  {
    if (direction == -1)
    {
      x = world_width - 1;
    }
    else
    {
      x = 0;
    }
  }
  for (i = 0; i < world_width; i+= direction)
  {
    if (world[(x + i + world_width) % world_width][0].plant_no > -1)
    {
      if (world[(x + i + world_width) % world_width][0].cell_no == 0)
      {
        return (world[(x + i + world_width) % world_width][0].plant_no);
      }
    }
  }
  return (-1);
}


/*
 * update the info bar of single plant window of plant_no. It is assumed
 * that this plant does have a single plant window, otherwise be prepared
 * for a crash...
 */

void update_single_info(long plant_no)
{
  sprintf(plant_wdescr[plant_no]->info, "%1ld cells, energy=%1ld, age=%1ld, gl=%1ld",
          plant[plant_no]->num_cells, plant[plant_no]->cellular_energy, plant[plant_no]->age, plant[plant_no]->genome.len);
  wind_set(plant_wdescr[plant_no]->whandle, WF_INFO, plant_wdescr[plant_no]->info);
}


void redraw_history_lane(int bufsize, int count, long first_generation, int base_x, int base_y,
                         char buffer_full, char *buffer,
                         int g_start, int g_end, int x, int y, int w, int h)
{
  int xy[10], g;

  xy[0] = base_x - 1;
  xy[1] = base_y - HWINDOW_MAX - 1;
  xy[2] = xy[0];
  xy[3] = base_y + 1;
  xy[4] = base_x + bufsize + 1;
  xy[5] = xy[3];
  xy[6] = xy[4];
  xy[7] = xy[1];
  xy[8] = xy[0];
  xy[9] = xy[1];
  
  if ((x <= xy[4]) && (xy[0] < x + w) && (y <= xy[5]) && (xy[1] < y + h))
  {
    v_pline(v_screenhandle, 5, xy);
    /* printf("plotting for positions %d ... %d\n", g_start - 1, g_end); */
    if (count > -1)
    {
      if (buffer_full)
      {
        /* printf("buffer is full\n"); */
        for (g = g_start; g < g_end; g++)
        {
          xy[0] = base_x + g - 1;
          xy[1] = base_y - buffer[(count + g) % bufsize];
          xy[2] = base_x + g;
          xy[3] = base_y - buffer[(count + g + 1) % bufsize];
          v_pline(v_screenhandle, 2, xy);
        }
      }
      else
      {
        if (count == 0)
        {
          /* printf("there's just one point\n"); */
          xy[0] = base_x;
          xy[1] = base_y - buffer[0];
          xy[2] = xy[0];
          xy[3] = xy[1];
          v_pline(v_screenhandle, 2, xy);
        }
        else
        {
          /* printf("buffer not yet full\n"); */
          for (g = g_start; g < g_end; g++)
          {
            xy[0] = base_x + g - 1;
            xy[1] = base_y - buffer[g - 1];
            xy[2] = base_x + g;
            xy[3] = base_y - buffer[g];
            v_pline(v_screenhandle, 2, xy);
          }
        }
      }
    }
  }
}


void redraw_single_plant(long plant_no, int x, int y, int w, int h, long xpos, long ypos)
{
  int x0, y0, w0, h0, xy[4];
  long xp, yp, xmin, ymin, xmax, ymax;

  /*printf("redrawing single plant #%1ld\n", plant_no); */
  wind_get(plant_wdescr[plant_no]->whandle, WF_WORKXYWH, &x0, &y0, &w0, &h0);
  vswr_mode(v_screenhandle, MD_REPLACE);
  vsl_width(v_screenhandle, 1);
  vsl_type(v_screenhandle, LT_SOLID);
  vsl_color(v_screenhandle, BLACK);
  vsf_color(v_screenhandle, BLACK);
  xmin = (xpos + x - x0) / plant_wdescr[plant_no]->cell_size;
  ymin = world_height - (ypos + y - y0 + h - 1) / plant_wdescr[plant_no]->cell_size - 1;
  xmax = xmin + (w - 1) / plant_wdescr[plant_no]->cell_size + 1;
  ymax = ymin + (h - 1) / plant_wdescr[plant_no]->cell_size + 1;
  /* printf("redrawing x=%1ld ... %1ld, y=%1ld ... %1ld\n", xmin, xmax, ymin, ymax); */
  for (yp = ymin; yp < ymax; yp++)
  {
    xy[1] = y0 + (world_height - yp - 1) * plant_wdescr[plant_no]->cell_size - ypos;
    xy[3] = xy[1] + plant_wdescr[plant_no]->cell_size - 1;
    for (xp = xmin; xp < xmax; xp++)
    {
      if (world[(xp + plant[plant_no]->cell[0].x + world_width / 2 + world_width) % world_width][yp].plant_no > -1)
      {
        xy[0] = x0 + xp * plant_wdescr[plant_no]->cell_size - xpos;
        xy[2] = xy[0] + plant_wdescr[plant_no]->cell_size - 1;
        if (world[(xp + plant[plant_no]->cell[0].x + world_width / 2 + world_width) % world_width][yp].plant_no == plant_no)
        {
          if (plant[plant_no]->cell[world[(xp + plant[plant_no]->cell[0].x + world_width / 2 + world_width) % world_width][yp].cell_no].energy)
          {
            vsf_interior(v_screenhandle, FIS_SOLID);
          }
          else
          {
            vsf_interior(v_screenhandle, FIS_PATTERN);
            vsf_style(v_screenhandle, 4);
          }
        }
        else
        {
          vsf_interior(v_screenhandle, FIS_HOLLOW);
        }
        v_bar(v_screenhandle, xy);
      }
    }
  }
  /* printf("done\n"); */
}


void redraw_genome(long plant_no, int x, int y, int w, int h, long xpos, long ypos)
{
  int  char_width, char_height, cell_width, cell_height, x0, y0, w0, h0, j, txtpos;
  long first_gene, last_gene, i;
  unsigned char gene_output;

  vst_color(v_screenhandle, BLACK);
  vst_height(v_screenhandle, default_char_height, &char_width, &char_height, &cell_width, &cell_height);
  wind_get(plant_wdescr[plant_no]->whandle, WF_WORKXYWH, &x0, &y0, &w0, &h0);
  first_gene = (ypos + y - y0) / cell_height;
  last_gene = (ypos + y + h - y0) / cell_height;
  last_gene = (last_gene < plant[plant_no]->genome.len) ? last_gene : plant[plant_no]->genome.len;
  for (i = first_gene; i <= last_gene; i++)
  {
    txtpos = 5;
    sprintf(buf, "%4ld: ", i);
    if (plant[plant_no]->genome.usg_count[i])
    {
      vst_effects(v_screenhandle, TF_NORMAL);
    }
    else
    {
      vst_effects(v_screenhandle, TF_LIGHTENED);
    }
    v_gtext(v_screenhandle, x0 - xpos + txtpos, y0 - ypos + i * cell_height + char_height, buf);
    txtpos += strlen(buf) * cell_width;
    for (j = 7; j >= 0; j--)
    {
      if (plant[plant_no]->genome.g[GENE_LENGTH * i] & (1 << j))
      {
        buf[7 - j] = '1';
      }
      else
      {
        buf[7 - j] = '0';
      }
    }
    for (j = 7; j >= 0; j--)
    {
      if (!(plant[plant_no]->genome.g[GENE_LENGTH * i + 1] & (1 << j)))
      {
        buf[7 - j] = '*';
      }
    }
    sprintf(buf + 8, " -> ");
    vst_effects(v_screenhandle, TF_NORMAL);
    v_gtext(v_screenhandle, x0 - xpos + txtpos, y0 - ypos + i * cell_height + char_height, buf);
    txtpos += strlen(buf) * cell_width;
    gene_output = plant[plant_no]->genome.g[GENE_LENGTH * i + 2];
    if ((gsys_parameters.divide_code <= gene_output) && (gene_output < gsys_parameters.statebit_code))
    {
      vst_effects(v_screenhandle, TF_THICKENED);
      sprintf(buf, "divide %1d", (int) gene_output & 0x07);
    }
    if ((gsys_parameters.statebit_code <= gene_output) && (gene_output < gsys_parameters.broadcastbit_code))
    {
      sprintf(buf, "........");
      buf[7 - (gene_output & 0x07)] = '1';
    }
    if ((gsys_parameters.broadcastbit_code <= gene_output) && (gene_output < gsys_parameters.flyingseed_code))
    {
      vst_color(v_screenhandle, BLUE);
      sprintf(buf, "\'\'\'\'\'\'\'\'");
      buf[7 - (gene_output & 0x07)] = 'b';
    }
    if ((gsys_parameters.flyingseed_code <= gene_output) && (gene_output < gsys_parameters.localseed_code))
    {
      vst_color(v_screenhandle, RED);
      vst_effects(v_screenhandle, TF_THICKENED | TF_UNDERLINED);
      sprintf(buf, "flying seed");
    }
    if ((gsys_parameters.localseed_code <= gene_output) && (gene_output < gsys_parameters.mutminus_code))
    {
      vst_color(v_screenhandle, RED);
      vst_effects(v_screenhandle, TF_THICKENED | TF_UNDERLINED);
      sprintf(buf, "local seed");
    }
    if ((gsys_parameters.mutminus_code <= gene_output) && (gene_output < gsys_parameters.mutplus_code))
    {
      vst_color(v_screenhandle, GREEN);
      vst_effects(v_screenhandle, TF_THICKENED | TF_SLANTED);
      sprintf(buf, "mut-");
    }
    if (gsys_parameters.mutplus_code <= gene_output)
    {
      vst_color(v_screenhandle, GREEN);
      vst_effects(v_screenhandle, TF_THICKENED | TF_SLANTED);
      sprintf(buf, "mut+");
    }
    v_gtext(v_screenhandle, x0 - xpos + txtpos, y0 - ypos + i * cell_height + char_height, buf);
    txtpos += strlen(buf) * cell_width;
    if (plant[plant_no]->genome.bp_count)
    {
      sprintf(buf, "[%4ld / %5ld]", plant[plant_no]->genome.usg_count[i], plant[plant_no]->genome.bp_count[i]);
    }
    else
    {
      sprintf(buf, "[%4ld / -----]", plant[plant_no]->genome.usg_count[i]);
    }
    txtpos = 5 + 31 * cell_width;
    vst_color(v_screenhandle, BLACK);
    vst_effects(v_screenhandle, TF_NORMAL);
    v_gtext(v_screenhandle, x0 - xpos + 32 * cell_width, y0 - ypos + i * cell_height + char_height, buf);
  }
}


void redraw_plant_history(long plant_no, int x, int y, int w, int h, long xpos, long ypos)
{
  int  char_width, char_height, cell_width, cell_height, x0, y0, w0, h0;
  int  count, g_start, g_end;

  wind_get(plant_wdescr[plant_no]->whandle, WF_WORKXYWH, &x0, &y0, &w0, &h0);
  vst_height(v_screenhandle, default_char_height, &char_width, &char_height, &cell_width, &cell_height);
  vst_effects(v_screenhandle, TF_NORMAL);
  vsl_width(v_screenhandle, 1);
  vsl_type(v_screenhandle, LT_SOLID);
  vsl_color(v_screenhandle, BLACK);
  count = plant_wdescr[plant_no]->hbuf.num_generations;
  g_start = xpos + x - x0 - 12 * cell_width - 1;
  g_start = (g_start > 1) ? g_start : 1;
  g_end = xpos + x + w - x0 - 12 * cell_width;
  if (plant_wdescr[plant_no]->hbuf.buffer_full)
  {
    g_end = (g_end < HISTORY_BUFSIZE) ? g_end : HISTORY_BUFSIZE;
  }
  else
  {
    g_end = (g_end < plant_wdescr[plant_no]->hbuf.num_generations + 1) ? g_end : plant_wdescr[plant_no]->hbuf.num_generations + 1;
  }
  if (xpos + x - x0 < 11 * cell_width)
  {
    v_gtext(v_screenhandle, x0 - xpos + 5, y0 - ypos + 5 + 0 * LANE_WIDTH + HWINDOW_MAX / 2 + char_height, "       #cells");
    v_gtext(v_screenhandle, x0 - xpos + 5, y0 - ypos + 5 + 1 * LANE_WIDTH + HWINDOW_MAX / 2 + char_height, "       energy");
    v_gtext(v_screenhandle, x0 - xpos + 5, y0 - ypos + 5 + 2 * LANE_WIDTH + HWINDOW_MAX / 2 + char_height, "   #mutations");
    v_gtext(v_screenhandle, x0 - xpos + 5, y0 - ypos + 5 + 3 * LANE_WIDTH + HWINDOW_MAX / 2 + char_height, "genome length");
  }
  redraw_history_lane(plant_wdescr[plant_no]->hbuf.bufsize, plant_wdescr[plant_no]->hbuf.num_generations,
          plant_wdescr[plant_no]->hbuf.first_generation, x0 - xpos + 15 * cell_width + 1, y0 - ypos + 6 + HWINDOW_MAX,
          plant_wdescr[plant_no]->hbuf.buffer_full,
          plant_wdescr[plant_no]->hbuf.num_cells, g_start, g_end, x, y, w, h);
  redraw_history_lane(plant_wdescr[plant_no]->hbuf.bufsize, plant_wdescr[plant_no]->hbuf.num_generations,
          plant_wdescr[plant_no]->hbuf.first_generation, x0 - xpos + 15 * cell_width + 1, y0 - ypos + LANE_WIDTH + 6 + HWINDOW_MAX,
          plant_wdescr[plant_no]->hbuf.buffer_full,
          plant_wdescr[plant_no]->hbuf.energy, g_start, g_end, x, y, w, h);
  redraw_history_lane(plant_wdescr[plant_no]->hbuf.bufsize, plant_wdescr[plant_no]->hbuf.num_generations,
          plant_wdescr[plant_no]->hbuf.first_generation, x0 - xpos + 15 * cell_width + 1, y0 - ypos + 2 * LANE_WIDTH + 6 + HWINDOW_MAX,
          plant_wdescr[plant_no]->hbuf.buffer_full,
          plant_wdescr[plant_no]->hbuf.num_mutations, g_start, g_end, x, y, w, h);
  redraw_history_lane(plant_wdescr[plant_no]->hbuf.bufsize, plant_wdescr[plant_no]->hbuf.num_generations,
          plant_wdescr[plant_no]->hbuf.first_generation, x0 - xpos + 15 * cell_width + 1, y0 - ypos + 3 * LANE_WIDTH + 6 + HWINDOW_MAX,
          plant_wdescr[plant_no]->hbuf.buffer_full,
          plant_wdescr[plant_no]->hbuf.genome_length, g_start, g_end, x, y, w, h);
}


/*
 * redraw a plant window.
 * Note: redraw_pwindow() expects the existence of a single plant window
 * for the specified plant and does not verify that. redraw_pwindow() is
 * called only through wndhndlr redraw functions, and its address is only
 * passed to create_window() upon creation of a single plant window, hence
 * the expectation should always hold.
 */

void redraw_pwindow(int w_handle, int x, int y, int w, int h, long xpos, long ypos)
{
  long plant_no;

#ifdef UNLOCK_SCREEN
  wind_update(END_UPDATE);
#endif
  /* printf("redrawing plant window #%1d\n", w_handle); */
  for (plant_no = 0; plant_no < world_width; plant_no++)
  {
    if (plant_wdescr[plant_no] != NULL)
    {
      if (plant_wdescr[plant_no]->whandle == w_handle)
      {
        break;
      }
    }
  }
  if (plant_no == world_width)
  {
    sprintf(buf, "[1][Error: Attempt to|redraw plant window|#%1d, for which no plant|can be found][ Ok ]",
            w_handle);
    form_alert(1, buf);
    return;
  }
  switch (plant_wdescr[plant_no]->contents)
  {
    case pw_plant:
      redraw_single_plant(plant_no, x, y, w, h, xpos, ypos);
      break;
    case pw_genome:
      redraw_genome(plant_no, x, y, w, h, xpos, ypos);
      break;
    case pw_history:
      redraw_plant_history(plant_no, x, y, w, h, xpos, ypos);
      break;
  }
#ifdef UNLOCK_SCREEN
  wind_update(BEG_UPDATE);
#endif
}


void redraw_world(int w_handle, int x, int y, int w, int h, long xpos, long ypos)
{
  long xp, yp;
  int xy[4], outline_xy[10], x0, y0, w0, h0;

#ifdef UNLOCK_SCREEN
  wind_update(END_UPDATE);
#endif
  if (display_active)
  {
    wind_get(w_handle, WF_WORKXYWH, &x0, &y0, &w0, &h0);
    vswr_mode(v_screenhandle, MD_REPLACE);
    vsl_width(v_screenhandle, 1);
    vsl_type(v_screenhandle, LT_SOLID);
    vsl_color(v_screenhandle, BLACK);
    vsf_color(v_screenhandle, BLACK);
    for (yp = 0; yp < world_height; yp++)
    {
      xy[1] = y0 - ypos + (world_height - yp - 1) * celldisplay_size;
      xy[3] = xy[1] + celldisplay_size - 1;
      for (xp = 0; xp < world_width; xp++)
      {
        if (world[xp][yp].plant_no > -1)
        {
          xy[0] = x0 - xpos + xp * celldisplay_size;
          xy[2] = xy[0] + celldisplay_size - 1;
          if (world[xp][yp].plant_no == highlit_plant)
          {
            vsl_width(v_screenhandle, 2);
            outline_xy[0] = xy[0] + 1;
            outline_xy[1] = xy[1] + 1;
            outline_xy[2] = xy[2] - 1;
            outline_xy[3] = xy[1] + 1;
            outline_xy[4] = xy[2] - 1;
            outline_xy[5] = xy[3] - 1;
            outline_xy[6] = xy[0] + 1;
            outline_xy[7] = xy[3] - 1;
            outline_xy[8] = outline_xy[0];
            outline_xy[9] = outline_xy[1];
          }
          if (plant[world[xp][yp].plant_no]->cell[world[xp][yp].cell_no].energy)
          {
            vsf_interior(v_screenhandle, FIS_PATTERN);
            vsf_style(v_screenhandle, 4);
          }
          else
          {
            vsf_interior(v_screenhandle, FIS_HOLLOW);
          }
          v_bar(v_screenhandle, xy);
          if ((world[xp][yp].plant_no == highlit_plant) && (highlit_plant > -1))
          {
            v_pline(v_screenhandle, 5, outline_xy);
          }
        }
      }
    }
  }
#ifdef UNLOCK_SCREEN
  wind_update(BEG_UPDATE);
#endif
}


/*
 * close a single plant window. The corresponding window descriptor is free()'d.
 * Note: plant_window_closed() expects the existence of the plant_wdescr array.
 */

int plant_window_closed(int w_handle)
{
  long i;

  for (i = 0; i < world_width; i++)
  {
    if (plant_wdescr[i] != NULL)
    {
      if (plant_wdescr[i]->whandle == w_handle)
      {
        wind_close(w_handle);
        free(plant_wdescr[i]);
        plant_wdescr[i] = NULL;
        return (1);
      }
    }
  }
  return (0);
}


int world_window_closed(int w_handle)
{
  if (display_shutdown)
  {
    wind_close(w_handle);
    world_whandle = 0;
    display_active = 0;
    return (1);
  }
  sprintf(buf, "[1][Do you really|want to finish|simulation|%s][ Yes | No ]", simname);
  if (form_alert(2, buf) == 1)
  {
    finish_flag = 1;
  }
  return (0);
}


/* forward declaration of poll_user_interface, which is in this module */

extern void poll_user_interface(void);


void create_plant_window(long plant_no)
{
  static int pw_xpos = INT_MAX - 40, pw_ypos = INT_MAX;
  int x, y, w, h, x0, y0, w0, h0;

  if (plant_wdescr == NULL)
  {
    if (form_alert(1, "[2][Memory for internal|single plant window|descriptors could not|be allocated. Do you|want to...][Retry|Abort]") == 1)
    {
      alloc_plant_wdescr();
    }
  }
  if (plant_wdescr != NULL)
  {
    if (plant_wdescr[plant_no] != NULL)
    {
      wind_set(plant_wdescr[plant_no]->whandle, WF_TOP);
      return;
    }
    plant_wdescr[plant_no] = (PLANT_WDESCR *) malloc(sizeof(PLANT_WDESCR));
    if (plant_wdescr[plant_no] == NULL)
    {
      form_alert(1,"[1][Sorry: Out of|memory. Cannot open|a plant window][ Ok ]");
      return;
    }
    plant_wdescr[plant_no]->cell_size = SINGLEPLANT_CELLSIZE;
    plant_wdescr[plant_no]->whandle = create_window(0xfff, world_width * plant_wdescr[plant_no]->cell_size, world_height * plant_wdescr[plant_no]->cell_size,
                                             v_screenhandle, redraw_pwindow, plant_window_closed);
    if (plant_wdescr[plant_no]->whandle <= 0)
    {
      form_alert(1, "[1][Sorry: No more|windows available.|Cannot open a|plant window.][ Ok ]");
      free(plant_wdescr[plant_no]);
      plant_wdescr[plant_no] = NULL;
      return;
    }
    plant_wdescr[plant_no]->contents = pw_plant;
    plant_wdescr[plant_no]->single_step = 0;
    plant_wdescr[plant_no]->hbuf.bufsize = HISTORY_BUFSIZE;
    plant_wdescr[plant_no]->hbuf.num_generations = -1;
    plant_wdescr[plant_no]->hbuf.first_generation = generation;
    plant_wdescr[plant_no]->hbuf.buffer_full = 0;
    sprintf(plant_wdescr[plant_no]->title, "plant #%1ld", plant_no);
    wind_set(plant_wdescr[plant_no]->whandle, WF_NAME, plant_wdescr[plant_no]->title);
    update_single_info(plant_no);
    wind_calc(WC_WORK, 0xfff, dw_x, dw_y, dw_w, dw_h, &x0, &y0, &w0, &h0);
    w0 = w0 / 2;
    w0 = (world_width * plant_wdescr[plant_no]->cell_size < w0) ? world_width * plant_wdescr[plant_no]->cell_size : w0;
    h0 = h0 / 2;
    h0 = (world_height * plant_wdescr[plant_no]->cell_size < h0) ? world_height * plant_wdescr[plant_no]->cell_size : h0;
    set_wsize_snap(plant_wdescr[plant_no]->whandle, plant_wdescr[plant_no]->cell_size, plant_wdescr[plant_no]->cell_size);
    set_arrow_step(plant_wdescr[plant_no]->whandle, plant_wdescr[plant_no]->cell_size, plant_wdescr[plant_no]->cell_size);
    set_vplane_alignment(plant_wdescr[plant_no]->whandle, 1, 1);
    set_hslider(plant_wdescr[plant_no]->whandle, (world_width * plant_wdescr[plant_no]->cell_size - w0) / 2, w0);
    set_vslider(plant_wdescr[plant_no]->whandle, world_height * plant_wdescr[plant_no]->cell_size - h0, h0);
    wind_calc(WC_BORDER, 0xfff, x0, y0, w0, h0, &x, &y, &w, &h);
    if (pw_ypos > dw_h - h)
    {
      pw_ypos = dw_y;
      pw_xpos += 40;
      if (pw_xpos > dw_w - w)
      {
        pw_xpos = dw_x;
      }
    }
    wind_open(plant_wdescr[plant_no]->whandle, pw_xpos, pw_ypos, w, h);
    pw_ypos += 40;
    poll_user_interface();
  }
}


void create_world_window(void)
{
  int x, y, w, h, x0, y0, w0, h0;

  wind_calc(WC_WORK, 0xfff, dw_x, dw_y, dw_w, dw_h, &x, &y, &w, &h);
  celldisplay_size = w / world_width;
  celldisplay_size = (h / world_height < celldisplay_size) ? h / world_height : celldisplay_size;
  celldisplay_size = (celldisplay_size > 1) ? celldisplay_size : 1;
  world_whandle = create_window(0xfff, world_width * celldisplay_size, world_height * celldisplay_size,
                                v_screenhandle, redraw_world, world_window_closed);
  if (world_whandle > 0)
  {
    wind_get(world_whandle, WF_FULLXYWH, &x, &y, &w, &h);
    wind_calc(WC_WORK, 0xfff, x, y, w, h, &x0, &y0, &w0, &h0);
    w0 = (w0 < world_width * celldisplay_size) ? w0 : world_width * celldisplay_size;
    h0 = (h0 < world_height * celldisplay_size) ? h0 : world_height * celldisplay_size;
    sprintf(title, " %s ", simname);
    wind_set(world_whandle, WF_NAME, title);
    info[0] = '\0';
    wind_set(world_whandle, WF_INFO, info);
    set_hslider(world_whandle, 0, w0);
    set_vslider(world_whandle, world_height * celldisplay_size - h0, h0);
    set_wsize_snap(world_whandle, celldisplay_size, celldisplay_size);
    set_arrow_step(world_whandle, celldisplay_size, celldisplay_size);
    set_vplane_alignment(world_whandle, 1, 1);
    wind_calc(WC_BORDER, 0xfff, x0, y0, w0, h0, &x, &y, &w, &h);
    wind_open(world_whandle, x, y, w, h);
    poll_user_interface();
  }
}


void draw_singleplant_cell(long xp, long yp)
{
  long plant_no, xpos, ypos;
  int x0, y0, w0, h0, x, y, w, h, cell_xy[4], clip_xy[4];

  if (plant_wdescr != NULL)
  {
    graf_mouse(M_OFF, NULL);
#ifndef UNLOCK_SCREEN
    wind_update(BEG_UPDATE);
#endif
    for (plant_no = 0; plant_no < world_width; plant_no++)
    {
      if (plant_wdescr[plant_no] != NULL)
      {
        if (plant_wdescr[plant_no]->contents == pw_plant)
        {
          wind_get(plant_wdescr[plant_no]->whandle, WF_WORKXYWH, &x0, &y0, &w0, &h0);
          get_slider_positions(plant_wdescr[plant_no]->whandle, &xpos, &ypos);
          cell_xy[0] = x0 - xpos + ((xp - plant[plant_no]->cell[0].x + world_width / 2 + world_width) % world_width) * plant_wdescr[plant_no]->cell_size;
          cell_xy[1] = y0 - ypos + (world_height - yp - 1) * plant_wdescr[plant_no]->cell_size;
          cell_xy[2] = cell_xy[0] + plant_wdescr[plant_no]->cell_size - 1;
          cell_xy[3] = cell_xy[1] + plant_wdescr[plant_no]->cell_size - 1;
          if ((x0 <= cell_xy[2]) && (cell_xy[0] < x0 + w0) && (y0 <= cell_xy[3]) && (cell_xy[1] < y0 + h0))
          {
            if (world[xp][yp].plant_no > -1)
            {
              vsf_color(v_screenhandle, BLACK);
              if (world[xp][yp].plant_no == plant_no)
              {
                if (plant[plant_no]->cell[world[xp][yp].cell_no].energy)
                {
                  vsf_interior(v_screenhandle, FIS_SOLID);
                }
                else
                {
                  vsf_interior(v_screenhandle, FIS_PATTERN);
                  vsf_style(v_screenhandle, 4);
                }
              }
              else
              {
                vsf_interior(v_screenhandle, FIS_HOLLOW);
              }
            }
            else
            {
              vsf_interior(v_screenhandle, FIS_SOLID);
              vsf_color(v_screenhandle, WHITE);
            }
            wind_get(plant_wdescr[plant_no]->whandle, WF_FIRSTXYWH, &x, &y, &w, &h);
            while ((w > 0) || (h > 0))
            {
              if ((x <= cell_xy[2]) && (cell_xy[0] < x + w) && (y <= cell_xy[3]) && (cell_xy[1] < y + h))
              {
                if (rc_inter(dw_x, dw_y, dw_w, dw_h, &x, &y, &w, &h))
                {
                  clip_xy[0] = x;
                  clip_xy[1] = y;
                  clip_xy[2] = x + w - 1;
                  clip_xy[3] = y + h - 1;
                  vs_clip(v_screenhandle, 1, clip_xy);
                  v_bar(v_screenhandle, cell_xy);
                }
              }
              wind_get(world_whandle, WF_NEXTXYWH, &x, &y, &w, &h);
            }
          }
        }
      }
    }
#ifndef UNLOCK_SCREEN
    wind_update(END_UPDATE);
#endif
    graf_mouse(M_ON, NULL);
    if ((world[xp][yp].plant_no > -1) && (plant_wdescr[world[xp][yp].plant_no] != NULL))
    {
      update_single_info(world[xp][yp].plant_no);
    }
  }
}


void draw_cell(long xp, long yp)
{
  int x, y, w, h, x0, y0, w0, h0, cell_xy[4], outline_xy[10], clip_xy[4];
  long xpos, ypos;

  /* VDI attibutes set also for single plant windows --
   * don't move under any conditions */
  vswr_mode(v_screenhandle, MD_REPLACE);
  vsl_width(v_screenhandle, 1);
  vsl_type(v_screenhandle, LT_SOLID);
  vsl_color(v_screenhandle, BLACK);
  if (world_whandle > 0)
  {
    wind_get(world_whandle, WF_WORKXYWH, &x0, &y0, &w0, &h0);
    get_slider_positions(world_whandle, &xpos, &ypos);
    cell_xy[0] = x0 - xpos + xp * celldisplay_size;
    cell_xy[2] = cell_xy[0] + celldisplay_size - 1;
    cell_xy[1] = y0 - ypos + (world_height - yp - 1) * celldisplay_size;
    cell_xy[3] = cell_xy[1] + celldisplay_size - 1;
    if ((x0 <= cell_xy[2]) && (cell_xy[0] < x0 + w0) && (y0 <= cell_xy[3]) && (cell_xy[1] < y0 + h0))
    {
      if (world[xp][yp].plant_no > -1)
      {
        vsf_color(v_screenhandle, BLACK);
        if (world[xp][yp].plant_no == highlit_plant)
        {
          vsl_width(v_screenhandle, 2);
          outline_xy[0] = cell_xy[0] + 1;
          outline_xy[1] = cell_xy[1] + 1;
          outline_xy[2] = cell_xy[2] - 1;
          outline_xy[3] = cell_xy[1] + 1;
          outline_xy[4] = cell_xy[2] - 1;
          outline_xy[5] = cell_xy[3] - 1;
          outline_xy[6] = cell_xy[0] + 1;
          outline_xy[7] = cell_xy[3] - 1;
          outline_xy[8] = outline_xy[0];
          outline_xy[9] = outline_xy[1];
        }
        if (plant[world[xp][yp].plant_no]->cell[world[xp][yp].cell_no].energy)
        {
          vsf_interior(v_screenhandle, FIS_PATTERN);
          vsf_style(v_screenhandle, 4);
        }
        else
        {
          vsf_interior(v_screenhandle, FIS_HOLLOW);
        }
      }
      else
      {
        vsf_color(v_screenhandle, WHITE);
        vsf_interior(v_screenhandle, FIS_SOLID);
      }
#ifndef UNLOCK_SCREEN
      wind_update(BEG_UPDATE);
#endif
      graf_mouse(M_OFF, NULL);
      wind_get(world_whandle, WF_FIRSTXYWH, &x, &y, &w, &h);
      while ((w > 0) || (h > 0))
      {
        if ((x <= cell_xy[2]) && (cell_xy[0] < x + w) && (y <= cell_xy[3]) && (cell_xy[1] < y + h))
        {
          if (rc_inter(dw_x, dw_y, dw_w, dw_h, &x, &y, &w, &h))
          {
            clip_xy[0] = x;
            clip_xy[1] = y;
            clip_xy[2] = x + w - 1;
            clip_xy[3] = y + h - 1;
            vs_clip(v_screenhandle, 1, clip_xy);
            v_bar(v_screenhandle, cell_xy);
            if ((world[xp][yp].plant_no == highlit_plant) && (highlit_plant > -1))
            {
              v_pline(v_screenhandle, 5, outline_xy);
            }
          }
        }
        wind_get(world_whandle, WF_NEXTXYWH, &x, &y, &w, &h);        
      }
      graf_mouse(M_ON, NULL);
#ifndef UNLOCK_SCREEN
      wind_update(END_UPDATE);
#endif
    }
  }
  draw_singleplant_cell(xp, yp);
}


void draw_plant(long plant_no)
{
  long i;

  for (i = 0; i < plant[plant_no]->num_cells; i++)
  {
    draw_cell(plant[plant_no]->cell[i].x, plant[plant_no]->cell[i].y);
  }
}


/*
 * update the history of values during the last couple of generations
 * that is displayed in the history window.
 * Note: update_history() expects the existence of a single plant window
 * for the specified plant and does not verify that.
 */

void update_history(long plant_no)
{
  long tmp;
  int  g;

  if (++plant_wdescr[plant_no]->hbuf.num_generations == plant_wdescr[plant_no]->hbuf.bufsize)
  {
    plant_wdescr[plant_no]->hbuf.num_generations = 0;
    plant_wdescr[plant_no]->hbuf.buffer_full = 1;
  }
  if (plant_wdescr[plant_no]->hbuf.buffer_full)
  {
    plant_wdescr[plant_no]->hbuf.first_generation++;
  }
  g = plant_wdescr[plant_no]->hbuf.num_generations;
  tmp = plant[plant_no]->num_cells;
  plant_wdescr[plant_no]->hbuf.num_cells[g] = (tmp < HWINDOW_MAX) ? tmp : HWINDOW_MAX;
  tmp = plant[plant_no]->cellular_energy;
  plant_wdescr[plant_no]->hbuf.energy[g] = (tmp < HWINDOW_MAX) ? tmp : HWINDOW_MAX;
  tmp = plant[plant_no]->genome.num_mutations;
  plant_wdescr[plant_no]->hbuf.num_mutations[g] = (tmp < HWINDOW_MAX) ? tmp : HWINDOW_MAX;
  tmp = plant[plant_no]->genome.len * 0.01;
  plant_wdescr[plant_no]->hbuf.genome_length[g] = (tmp < HWINDOW_MAX) ? tmp : HWINDOW_MAX;
  if (plant_wdescr[plant_no]->contents == pw_history)
  {
    w_redraw_all(plant_wdescr[plant_no]->whandle);
  }
}


/*
 * sets plant number plant_no as the plant highlit in the main world window.
 */

void set_highlit_plant(long plant_no)
{
  long old_highlit_plant;

  if (plant_no != highlit_plant)
  {
    old_highlit_plant = highlit_plant;
    highlit_plant = plant_no;
    if (old_highlit_plant > -1)
    {
      draw_plant(old_highlit_plant);
    }
    draw_plant(highlit_plant);
  }
}


/*
 * cycle_contents() changes the contents of a single plant window. If the
 * window currently displays the plant, the contens are changed to the genome.
 * If the current contents are the genome, the history is displayed, and if
 * the history is currently displayed, the window contents are set to the
 * plant.
 * Note: cycle_contents() assumes that plant_wdescr[plant_no] points to a valid
 * plant window descriptor.
 */

void cycle_contents(long plant_no)
{
  int  char_width, char_height, cell_width, cell_height, x0, y0, w0, h0;

  wind_get(plant_wdescr[plant_no]->whandle, WF_WORKXYWH, &x0, &y0, &w0, &h0);
  switch (plant_wdescr[plant_no]->contents)
  {
  case pw_plant:
    plant_wdescr[plant_no]->contents = pw_genome;
    vst_height(v_screenhandle, default_char_height, &char_width, &char_height, &cell_width, &cell_height);
    set_vplane_size(plant_wdescr[plant_no]->whandle, 50 * cell_width, plant[plant_no]->genome.len * cell_height);
    set_wsize_snap(plant_wdescr[plant_no]->whandle, cell_width, cell_height);
    set_arrow_step(plant_wdescr[plant_no]->whandle, cell_width, cell_height);
    set_vplane_alignment(plant_wdescr[plant_no]->whandle, 1, 1);
    set_hslider(plant_wdescr[plant_no]->whandle, 0, w0);
    set_vslider(plant_wdescr[plant_no]->whandle, 0, h0);
    break;
  case pw_genome:
    plant_wdescr[plant_no]->contents = pw_history;
    vst_height(v_screenhandle, default_char_height, &char_width, &char_height, &cell_width, &cell_height);
    set_vplane_size(plant_wdescr[plant_no]->whandle, 15 * cell_width + HISTORY_BUFSIZE + 10, LANE_WIDTH * NUM_LANES + 10);
    set_wsize_snap(plant_wdescr[plant_no]->whandle, 1, 1);
    set_arrow_step(plant_wdescr[plant_no]->whandle, 10, 10);
    set_vplane_alignment(plant_wdescr[plant_no]->whandle, 0, 0);
    set_hslider(plant_wdescr[plant_no]->whandle, 0, w0);
    set_vslider(plant_wdescr[plant_no]->whandle, 0, h0);
    break;
  case pw_history:
    plant_wdescr[plant_no]->contents = pw_plant;
    set_vplane_size(plant_wdescr[plant_no]->whandle, world_width * plant_wdescr[plant_no]->cell_size, world_height * plant_wdescr[plant_no]->cell_size);
    set_wsize_snap(plant_wdescr[plant_no]->whandle, plant_wdescr[plant_no]->cell_size, plant_wdescr[plant_no]->cell_size);
    set_arrow_step(plant_wdescr[plant_no]->whandle, plant_wdescr[plant_no]->cell_size, plant_wdescr[plant_no]->cell_size);
    set_vplane_alignment(plant_wdescr[plant_no]->whandle, 1, 1);
    set_hslider(plant_wdescr[plant_no]->whandle, (world_width * plant_wdescr[plant_no]->cell_size - w0) / 2, w0);
    set_vslider(plant_wdescr[plant_no]->whandle, world_height * plant_wdescr[plant_no]->cell_size - h0, h0);
    break;
  }
  wind_get(plant_wdescr[plant_no]->whandle, WF_CURRXYWH, &x0, &y0, &w0, &h0);
  w_sized(plant_wdescr[plant_no]->whandle, x0, y0, w0, h0);
  w_redraw_all(plant_wdescr[plant_no]->whandle);
}


/*
 * respond to a mouseclick at (x, y).
 * If the mouseclick is in the main world window, the plant clicked to is
 * highlit (if any is clicked to).
 */

void do_mouseclick(int x, int y)
{
  long plant_no;
  int  w_handle, scratch;

  wind_get(0, WF_TOP, &w_handle, &scratch, &scratch, &scratch);
  if (w_handle == world_whandle)
  {
    plant_no = find_plant(x, y);
    if (plant_no > -1)
    {
      set_highlit_plant(plant_no);
    }
  }
}


void do_keypress(unsigned int ret_keybd, long kbshift)
{
  long plant_no;
  int  w_handle, scratch;
  char tmp[200];


  /* switch for cursor keys -- emulate window arrows */
  switch (ret_keybd & 0xff00)
  {
  case 0x4b00:  /* left arrow */
    wind_get(0, WF_TOP, &w_handle, &scratch, &scratch, &scratch);
    if (kbshift & 3)
      w_arrow(w_handle, WA_LFPAGE);
    else
      w_arrow(w_handle, WA_LFLINE);
    return;
  case 0x4d00:  /* right arrow */
    wind_get(0, WF_TOP, &w_handle, &scratch, &scratch, &scratch);
    if (kbshift & 3)
      w_arrow(w_handle, WA_RTPAGE);
    else
      w_arrow(w_handle, WA_RTLINE);
    return;
  case 0x4800:  /* up arrow */
    wind_get(0, WF_TOP, &w_handle, &scratch, &scratch, &scratch);
    if (kbshift & 3)
      w_arrow(w_handle, WA_UPPAGE);
    else
      w_arrow(w_handle, WA_UPLINE);
    return;
  case 0x5000:  /* down arrow */
    wind_get(0, WF_TOP, &w_handle, &scratch, &scratch, &scratch);
    if (kbshift & 3)
      w_arrow(w_handle, WA_DNPAGE);
    else
      w_arrow(w_handle, WA_DNLINE);
    return;
  }

  if (ret_keybd & 0xff)
  {
    /* switch for keys in ASCII orientated fashion */
    switch (ret_keybd & 0xff)
    {
    case 17: /* ^Q */
      sprintf(buf, "[1][Do you really|want to finish|simulation|%s][ Yes | No ]", simname);
      if (form_alert(2, buf) == 1)
      {
        finish_flag = 1;
      }
      return;
    case 'f': case 'F':
      wind_get(0, WF_TOP, &w_handle, &scratch, &scratch, &scratch);
      w_fulled(w_handle);
      return;
    case 'c': case 'C':
      wind_get(0, WF_TOP, &w_handle, &scratch, &scratch, &scratch);
      remove_window(w_handle);
      return;
    case 'p': case 'P':
      pause_flag = !pause_flag;
      return;
    case 'o': case 'O':
      if (highlit_plant > -1)
      {
        create_plant_window(highlit_plant);
      }
      return;
    case '?':
      sprintf(buf, "[1][Plant #%1ld: %1ld cells,|total energy = %1ld,|age = %1ld][ Ok ]",
              highlit_plant, plant[highlit_plant]->num_cells, plant[highlit_plant]->cellular_energy, plant[highlit_plant]->age);
      form_alert(1, buf);
      return;
    case '+':
      wind_get(0, WF_TOP, &w_handle, &scratch, &scratch, &scratch);
      if ((plant_no = whandle_to_plantno(w_handle)) > -1)
      {
        cycle_contents(plant_no);
      }
      return;
    default:
      sprintf(tmp, "[1][scancode + ascii|of key pressed: %x][ Ok ]", ret_keybd);
      form_alert(1, tmp);
      return;
    }
  }
  else
  {
    /* switch for keys that don't have an ASCII code */
    switch (ret_keybd)
    {
    case 0x3b00:  /* F1 */
      if (world_whandle > 0)
      {
        wind_set(world_whandle, WF_TOP);
      }
      else
      {
        create_world_window();
      }
      return;
    case 0x6200:  /* Help */
      form_alert(1,"[0][F1 - open/top main window|\004\003 - move selection| p  - toggle pause|^Q - quit][ Ok ]");
      return;
    case 0x5200:  /* Insert */
      plant_no = next_plant(highlit_plant, -1);
      set_highlit_plant(plant_no);
      return;
    case 0x4700:  /* ClrHome */
      plant_no = next_plant(highlit_plant, 1);
      set_highlit_plant(plant_no);
      return;
    default:
      sprintf(tmp, "[1][scancode of key|pressed: %x][ Ok ]", ret_keybd);
      form_alert(1, tmp);
      return;
    }
  }
}


void poll_user_interface(void)
{
  static unsigned long last_time = 0;

  int  mbuf[8];
  unsigned int mx, my, ks, ret_keybd;
  int  mb, ret_button, event;

  if ((hz200timer() - last_time) > 40)
  {
    last_time = hz200timer();
    do
    {
      if (pause_flag)
      {
        event = evnt_multi(MU_KEYBD | MU_MESAG | MU_BUTTON,
                           1, 1, 1,             /* wait for 1 click with left mouse button */
                           0, 0, 0, 0, 0,       /* dummy values for mouse events */
                           0, 0, 0, 0, 0,
                           mbuf,                /* buffer for AES messages */
                           (unsigned long) 0,   /* 0 timer ticks (return immediately) */
                           &mx, &my, &mb, &ks, &ret_keybd, &ret_button);
      }
      else
      {
        event = evnt_multi(MU_KEYBD | MU_MESAG | MU_TIMER | MU_BUTTON,
                           1, 1, 1,             /* wait for 1 click with left mouse button */
                           0, 0, 0, 0, 0,       /* dummy values for mouse events */
                           0, 0, 0, 0, 0,
                           mbuf,                /* buffer for AES messages */
                           (unsigned long) 0,   /* 0 timer ticks (return immediately) */
                           &mx, &my, &mb, &ks, &ret_keybd, &ret_button);
      }
      if (event & MU_KEYBD)
      {
        do_keypress(ret_keybd, Kbshift(-1));
        if (finish_flag)
        {
          pause_flag = 0;
          break;
        }
      }
      if (event & MU_MESAG)
      {
        process_message(mbuf);
      }
      if (event & MU_BUTTON)
      {
        evnt_button(1, 1, 0, &mx, &my, &mb, &ks);
        do_mouseclick(mx, my);
      }
    }
    while (((event != MU_TIMER) && (event != 0)) || pause_flag);
  }
}


/*
 * Update the generation display in the info bar of the main world window
 */

void display_start(long generation)
{
  if (world_whandle > 0)
  {
    display_active = 1;
    sprintf(info, "generation %5ld", generation);
    wind_set(world_whandle, WF_INFO, info);
    poll_user_interface();
    w_redraw_all(world_whandle);
  }
}


/*
 * a call to generation_done() indicates that no more changes in plant
 * structures etc. are to be expected before the next start_generation()
 */

void display_done(void)
{
  long plant_no;

  if (world_whandle > 0)
  {
    poll_user_interface();
    display_active = 0;
    /* w_redraw_all(world_whandle); */
  }
  if (plant_wdescr != NULL)
  {
    for (plant_no = 0; plant_no < world_width; plant_no++)
    {
      if (plant_wdescr[plant_no] != NULL)
      {
        update_history(plant_no);
      }
    }
  }
}


void display_cell(long xp, long yp)
{
  draw_cell(xp, yp);
  poll_user_interface();
}


void display_plant_killed(long plant_no)
{
  if (plant_no == highlit_plant)
  {
    highlit_plant = -1;
  }
  if (plant_wdescr != NULL)
  {
    if (plant_wdescr[plant_no] != NULL)
    {
      remove_window(plant_wdescr[plant_no]->whandle);
    }
  }
}


void display_plant_mutated(long plant_no)
{
  int char_width, char_height, cell_width, cell_height, x0, y0, w0, h0;

  if (plant_wdescr != NULL)
  {
    if (plant_wdescr[plant_no] != NULL)
    {
      update_single_info(plant_no);
      if (plant_wdescr[plant_no]->contents == pw_genome)
      {
        vst_height(v_screenhandle, default_char_height, &char_width, &char_height, &cell_width, &cell_height);
        set_vplane_size(plant_wdescr[plant_no]->whandle, 50 * cell_width, plant[plant_no]->genome.len * cell_height);
        wind_get(plant_wdescr[plant_no]->whandle, WF_CURRXYWH, &x0, &y0, &w0, &h0);
        w_sized(plant_wdescr[plant_no]->whandle, x0, y0, w0, h0);
        w_redraw_all(plant_wdescr[plant_no]->whandle);
      }
    }
  }
}


void init_world_display(const char *simname)
{
  int scratch;

  gem_init(&dw_x, &dw_y, &dw_w, &dw_h, &aes_screenhandle, &v_screenhandle, work_in, work_out);
  if (work_out[1] < 399)
  {
    default_char_height = 6;
  }
  else
  {
    default_char_height = 13;
  }
  vst_alignment(v_screenhandle, TA_LEFT, TA_BASELINE, &scratch, &scratch);
  alloc_plant_wdescr();
  create_world_window();
}


void close_world_display(void)
{
  display_shutdown = 1;
  kill_windows();
  if (plant_wdescr != NULL)
  {
    free(plant_wdescr);
    plant_wdescr = NULL;
  }
}

#endif /* __atarist__,  Atari specific user interface */

