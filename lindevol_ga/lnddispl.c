/* This module contains the functions for displaying the LindEvol world. */

#ifndef __atarist__
#  include "lnddasc.c"

#else /* For the Atari system, a graphical user interface is used. */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <osbind.h>
#include <aesbind.h>
#include <vdibind.h>
#include <gemfast.h>

#include <gfaport.h>

#include "lnderror.h"
#include "lndvals.h"
#include "lndtypes.h"
#include "wndhndlr.h"


/***** These variables are imported from lndevlxx.c *****/

extern char quietmode;
extern char worldmode;

extern char finish_flag;

extern unsigned long world_width, world_height; /* these are defined and */
extern unsigned long psize;
extern PLANT *plant;                            /* created in the main module */
extern LATTICE_SITE **world;
extern GENOME *genome;
extern unsigned long *gi;

/***** end of imported stuff *****/


char display_active = 0;

FILE  *world_file;
char   world_file_name[MAX_SLEN];

int ap_id, aes_screenhandle, v_screenhandle, dw_x, dw_y, dw_w, dw_h;
int work_in[11], work_out[57];

static int world_whandle = 0;
static char title[MAX_SLEN], info[MAX_SLEN];

static int celldisplay_size;

static const volatile unsigned long * const p_hz200timer = (unsigned long *) 0x4ba;
static long stack;


int open_world_file(const char *simname)
{
  sprintf(world_file_name, "%s%swld", simname, ".");
  world_file = fopen(world_file_name, "w");
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


unsigned long hz200timer(void)
{
  unsigned long t;

  stack = Super(0);
  t = *p_hz200timer;
  Super((void *) stack);
  return (t);
}


void redraw_world(int w_handle, int x, int y, int w, int h, long xpos, long ypos)
{
  unsigned long xp, yp;
  int xy[4], x0, y0, w0, h0;

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
          if (plant[world[xp][yp].plant_no].cell[world[xp][yp].cell_no].energy)
          {
            vsf_interior(v_screenhandle, FIS_SOLID);
          }
          else
          {
            vsf_interior(v_screenhandle, FIS_HOLLOW);
          }
          v_bar(v_screenhandle, xy);
        }
      }
    }
  }
}


void world_window_closed(int w_handle)
{
  world_whandle = 0;
  display_active = 0;
}


void display_on(unsigned long generation)
{
  if (world_whandle > 0)
  {
    display_active = 1;
    sprintf(info, "generation %5lu", generation);
    wind_set(world_whandle, WF_INFO, info);
  }
}


void display_off(void)
{
  if (world_whandle > 0)
  {
    display_active = 0;
    w_redraw_all(world_whandle);
  }
}


void poll_user_interface(void)
{
  static unsigned long last_timer = 0;

  int  mbuf[8];
  unsigned int mx, my, ks, ret_keybd;
  int  mb, ret_button, event;

  if (hz200timer() - last_timer > 40)
  {
    do
    {
      event = evnt_multi(MU_KEYBD | MU_MESAG | MU_TIMER,
                         0, 0, 0,
                         0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0,
                         mbuf, (unsigned long) 0,
                         &mx, &my, &mb, &ks, &ret_keybd, &ret_button);
      if (event & MU_KEYBD)
      {
        if ((ret_keybd & 0xff) == 17)
        {
          finish_flag = 1;
        }
      }
      if (event & MU_MESAG)
      {
        process_event(mbuf);
      }
    }
    while (event != MU_TIMER);
  }
}


void draw_cell(unsigned long xp, unsigned long yp)
{
  int x, y, w, h, x0, y0, w0, h0, cell_xy[4], clipxy[4];
  long xpos, ypos;

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
      vswr_mode(v_screenhandle, MD_REPLACE);
      vsl_width(v_screenhandle, 1);
      vsl_type(v_screenhandle, LT_SOLID);
      vsl_color(v_screenhandle, BLACK);
      vsf_color(v_screenhandle, BLACK);
      if (plant[world[xp][yp].plant_no].cell[world[xp][yp].cell_no].energy)
      {
        vsf_interior(v_screenhandle, FIS_SOLID);
      }
      else
      {
        vsf_interior(v_screenhandle, FIS_HOLLOW);
      }
      wind_update(BEG_UPDATE);
      wind_get(world_whandle, WF_FIRSTXYWH, &x, &y, &w, &h);
      while ((w > 0) || (h > 0))
      {
        if ((x <= cell_xy[2]) && (cell_xy[0] < x + w) && (y <= cell_xy[3]) && (cell_xy[1] < y + h))
        {
          if (rc_inter(dw_x, dw_y, dw_w, dw_h, &x, &y, &w, &h))
          {
            clipxy[0] = x;
            clipxy[1] = y;
            clipxy[2] = x + w - 1;
            clipxy[3] = y + h - 1;
            vs_clip(v_screenhandle, 1, clipxy);
            v_bar(v_screenhandle, cell_xy);
          }
        }
        wind_get(world_whandle, WF_NEXTXYWH, &x, &y, &w, &h);        
      }
      wind_update(END_UPDATE);
    }
    poll_user_interface();
  }
}


void init_world_display(const char *simname)
{
  int x, y, w, h, x0, y0, w0, h0;

  gem_init(&dw_x, &dw_y, &dw_w, &dw_h, &aes_screenhandle, &v_screenhandle, work_in, work_out);
  wind_calc(WC_WORK, 0xfff, dw_x, dw_y, dw_w, dw_h, &x, &y, &w, &h);
  celldisplay_size = w / world_width;
  celldisplay_size = (h / world_height < celldisplay_size) ? h / world_height : celldisplay_size;
  celldisplay_size = (celldisplay_size > 1) ? celldisplay_size : 1;
  world_whandle = create_window(0xfff, world_width * celldisplay_size, world_height * celldisplay_size,
                                1, 1, celldisplay_size, v_screenhandle, redraw_world, world_window_closed);
  if (world_whandle > 0)
  {
    wind_get(world_whandle, WF_FULLXYWH, &x, &y, &w, &h);
    wind_calc(WC_WORK, 0xfff, x, y, w, h, &x0, &y0, &w0, &h0);
    sprintf(title, " %s ", simname);
    wind_set(world_whandle, WF_NAME, title);
    info[0] = '\0';
    wind_set(world_whandle, WF_INFO, info);
    set_hslider(world_whandle, 0, w0);
    set_vslider(world_whandle, world_height * celldisplay_size - h0, h0);
    wind_open(world_whandle, x, y, w, h);
    poll_user_interface();
  }
}


void close_world_display(void)
{
  if (world_whandle > 0)
  {
    remove_window(world_whandle);
  }
}

#endif /* __atarist__,  Atari specific user interface */

