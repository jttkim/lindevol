/* $Id: lnddx6.c,v 1.2 2000/01/30 03:10:59 kim Exp $ */
/*
 * $Log: lnddx6.c,v $
 * Revision 1.2  2000/01/30 03:10:59  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

/* This module contains the functions for displaying the LindEvol world.
   At present, there is just a primitive ASCII displaying system. */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif

#include "lnderror.h"
#include "lndvals.h"
#include "lndtypes.h"

#include "lndglobl.h"
#include "lnddispl.h"


static long curr_generation;

FILE  *world_file = NULL;
char   world_file_name[MAX_SLEN];


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


void show_world(long generation, FILE *f)
{
  long x, y, l, p;

  fprintf(f, "world at generation %ld, sized %ld * %ld\n", curr_generation, world_width, world_height);
  for (y = world_height - 1;y >= 0 ; y--)
  {
    for (x = 0; x < world_width; x++)
    {
      if ((p = world[x][y].plant_no) > -1)
      {
        if (plant[p]->cell[world[x][y].cell_no].energy)
        {
          fprintf(f, "%c", (char) ('A' + p % 26));
        }
        else
        {
          fprintf(f, "%c", (char) ('a' + p % 26));
        }
      }
      else
      {
        fprintf(f, " ");
      }
    }
    fprintf(f, "\n");
  }
  fprintf(f, "\n");
  for (l = 0; l < 6; l++)
  {
    for (x = 0; x < l; x++) fprintf(f, " ");
    for (x = l; x < world_width; x+= 6)
    {
      if ((p = world[x][0].plant_no) > -1)
      {
        if (world[x][0].cell_no == 0)
        {
          fprintf(f, "%-5ld ", plant[p]->cellular_energy);
        }
        else
        {
          fprintf(f, "      ");
        }
      }
      else
      {
        fprintf(f, "      ");
      }
    }
    fprintf(f, "\n");
  }
  /* fprintf(f, "\nmax. fitness: %ld, min. fitness: %ld\n", genome[gi[0]].fitness, genome[gi[psize - 1]].fitness); */
}


void display_start(long generation)
{
  curr_generation = generation;
}


void display_done(void)
{
  if (!quietmode)
    show_world(curr_generation, stdout);
  if ((worldmode != 0) && (world_file != NULL))
    show_world(curr_generation, world_file);
}


void display_plant_killed(long plant_no)
{
}


void display_plant_mutated(long plant_no)
{
}


void display_cell(long x, long y)
{
}


void poll_user_interface(void)
{
}


void init_world_display(const char *simname)
{
  if (worldmode)
  {
    open_world_file(simname, "a");
  }
}


void close_world_display(void)
{
  if (worldmode)
  {
    close_world_file();
  }
}

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>

#include <sys/ipc.h>
#include <sys/shm.h>
#include <X11/extensions/XShm.h>

#include "wndhndlr.h"
#include "jklib.h"

#define NORMAL_PLAY   1
#define FAST_FORWARD  2
#define FAST_REWIND   3
#define BACKWARD_PLAY 4
#define SINGLE_STEP   5

Status XShmQueryExtension(Display *display);

char *event_name[] = {
"",
"",
"KeyPress",
"KeyRelease",
"ButtonPress",
"ButtonRelease",
"MotionNotify",
"EnterNotify",
"LeaveNotify",
"FocusIn",
"FocusOut",
"KeymapNotify",
"Expose",
"GraphicsExpose",
"NoExpose",
"VisibilityNotify",
"CreateNotify",
"DestroyNotify",
"UnmapNotify",
"MapNotify",
"MapRequest",
"ReparentNotify",
"ConfigureNotify",
"ConfigureRequest",
"GravityNotify",
"ResizeRequest",
"CirculateNotify",
"CirculateRequest",
"PropertyNotify",
"SelectionClear",
"SelectionRequest",
"SelectionNotify",
"ColormapNotify",
"ClientMessage",
"MappingNotify" };


static int mode = SINGLE_STEP;
static int game_over = 0;
static int screen_no;
static Window window = BadValue;
static Drawable paint_pixmap = BadValue, redraw_pixmap = BadValue;
static GC gc;
static GC lnd_gc[6], lndx_gc[6], background_gc[16], text_gc;
static XFontStruct *fixed_font;
static unsigned long lnd_pixel[6], lndx_pixel[6], back_pixel[16];
static Display *display;
static long frame_count = 0, day_count = 0;
long resume_fpos = 0;
long resume_frame = 0;
long world_width, world_height, generation;
int cell_boxsize = 0, window_width, window_height, pixmap_height;

Status mit_shm = False;


void handle_signal(int sig_id)
{
  switch (sig_id)
  {
  case SIGTERM: case SIGINT: case SIGHUP:
    signal(sig_id, handle_signal);
    fprintf(stderr, "xlmovie will terminate\n");
    game_over = 1;
    break;
  default:
    fprintf(stderr, "Caught signal %d but don't know what to do\n", sig_id);
    break;
  }
}

void delay(float delay)
{
  struct timeval t;

  t.tv_sec = floor(delay);
  t.tv_usec = floor((delay - t.tv_sec) * 1e6);
  select(0, NULL, NULL, NULL, &t);
}

void init_signal_handling(void)
{
  signal(SIGTERM, handle_signal);
  signal(SIGINT, handle_signal);
  signal(SIGHUP, handle_signal);
}


void redraw_text(void)
{
  char ts[256];

  switch (mode)
  {
  case FAST_FORWARD:
    sprintf(ts, "time step: %7ld (day %7ld) <fast forward> ", generation, day_count);
    break;
  case FAST_REWIND:
    sprintf(ts, "time step: %7ld (day %7ld) <fast rewind>  ", generation, day_count);
    break;
  case NORMAL_PLAY:
    sprintf(ts, "time step: %7ld (day %7ld) <forward play> ", generation, day_count);
    break;
  case BACKWARD_PLAY:
    sprintf(ts, "time step: %7ld (day %7ld) <backward play>", generation, day_count);
    break;
  case SINGLE_STEP:
    sprintf(ts, "time step: %7ld (day %7ld) <single step>  ", generation, day_count);
    break;
  default:
    sprintf(ts, "time step: %7ld (day %7ld)", generation, day_count);
    break;
  }
  if (fixed_font)
  {
    XDrawImageString(display, window, text_gc, 0, fixed_font->ascent, ts, strlen(ts));
    XDrawLine(display, window, text_gc, 0, window_height - pixmap_height - 2, window_width, window_height - pixmap_height - 2);
  }
}


void redraw_pixworld(Window win, int x, int y, unsigned int w, unsigned int h, long xpos, long ypos)
{
  if ((window != BadValue) && (redraw_pixmap != BadValue))
    XCopyArea(display, redraw_pixmap, window, background_gc[15], 0, 0, window_width, pixmap_height, 0, window_height - pixmap_height);
  redraw_text();
}


int close_pixwindow(Window win)
{
  if (win == window)
  {
    window = BadValue;
    game_over = 1;
    return (1);
  }
/*
  else
  {
    if (window != BadValue)
    {
      printf("sizeable window still exists -- no shutdown!\n");
      return (0);
    }
    printf("other window was closed\n");
    game_over = 1;
    return (1);
  }
*/
  return (1);
}


void read_pixheader(FILE *f)
{
  char buf[1000];

  fgets(buf, 100, f);
  fgets(buf, 100, f);
  fgets(buf, 100, f);
  world_width = strtol(buf, (char **) NULL, 10);
  fgets(buf, 100, f);
  world_height = strtol(buf, (char **) NULL, 10);
}


int create_lndgcs(GC *lnd_gc, GC *lndx_gc, GC *background_gc, GC *text_gc)
{
  XGCValues gcvalues;
  XColor xcolor;
  Colormap cmap;
  int i, j;

  cmap = DefaultColormap(display, screen_no);
  xcolor.flags = DoRed | DoGreen | DoBlue;
  gcvalues.line_width = 1;
  gcvalues.background = WhitePixel(display, screen_no);
  gcvalues.foreground = BlackPixel(display, screen_no);
  if (fixed_font)
    gcvalues.font = fixed_font->fid;
  *text_gc = XCreateGC(display, window, GCForeground | GCBackground | GCLineWidth | GCFont, &gcvalues);
  gcvalues.foreground = WhitePixel(display, screen_no);
  gcvalues.background = BlackPixel(display, screen_no);
  /* printf("WhitePixel=%lu, BlackPixel=%lu\n", (unsigned long) WhitePixel(display, screen_no), (unsigned long) BlackPixel(display, screen_no)); */
  back_pixel[15] = WhitePixel(display, screen_no);
  background_gc[15] = XCreateGC(display, window, GCForeground | GCBackground | GCLineWidth, &gcvalues);
  gcvalues.foreground = BlackPixel(display, screen_no);
  back_pixel[0] = BlackPixel(display, screen_no);
  background_gc[0] = XCreateGC(display, window, GCForeground | GCBackground | GCLineWidth, &gcvalues);
  xcolor.red = 32767;
  xcolor.green = 32767;
  xcolor.blue = 32767;
  if (!XAllocColor(display, cmap, &xcolor))
    return (-1);
  gcvalues.foreground = xcolor.pixel;
  back_pixel[1] = xcolor.pixel;
  background_gc[1] = XCreateGC(display, window, GCForeground | GCBackground | GCLineWidth, &gcvalues);
  for (i = 1; i < 7; i++)
  {
    if ((i & 1))
      xcolor.red = 65535;
    else
      xcolor.red = 0;
    if ((i & 2))
      xcolor.green = 65535;
    else
      xcolor.green = 0;
    if ((i & 4))
      xcolor.blue = 65535;
    else
      xcolor.blue = 0;
    if (!XAllocColor(display, cmap, &xcolor))
      return (-1);
    /* printf("Got color for %d: r=%d, g=%d, b=%d, pixel=%lu\n", i, xcolor.red, xcolor.green, xcolor.blue, xcolor.pixel); */
    gcvalues.foreground = xcolor.pixel;
    lnd_pixel[i - 1] = xcolor.pixel;
    lnd_gc[i - 1] = XCreateGC(display, window, GCForeground | GCBackground | GCLineWidth, &gcvalues);
    xcolor.red = xcolor.red * 3 / 4;
    xcolor.green = xcolor.green * 3 / 4;
    xcolor.blue = xcolor.blue * 3 / 4;
    if (!XAllocColor(display, cmap, &xcolor))
      return (-1);
    /* printf("Got color for x-%d: r=%d, g=%d, b=%d, pixel=%lu\n", i, xcolor.red, xcolor.green, xcolor.blue, xcolor.pixel); */
    gcvalues.foreground = xcolor.pixel;
    lndx_pixel[i - 1] = xcolor.pixel;
    lndx_gc[i - 1] = XCreateGC(display, window, GCForeground | GCBackground | GCLineWidth, &gcvalues);
  }
/*
  getchar();
  for (i = 0; i < 6; i++)
  {
    for (j = 0; j < window_height / 2; j++)
      XDrawLine(display, window, lnd_gc[i], 0, j, window_width - 1, j);
    for (; j < window_height; j++)
      XDrawLine(display, window, lndx_gc[i], 0, j, window_width - 1, j);
    printf("lndx_gc[%d]\n", i);
    getchar();
  }
  for (j = 0; j < window_height; j++)
    XDrawLine(display, window, background_gc[15], 0, j, window_width - 1, j);
  getchar();
*/
  return (0);
}


void free_lndgcs(void)
{
  int i;

  for (i = 0; i < 6; i++)
  {
    XFreeGC(display, lnd_gc[i]);
    XFreeGC(display, lndx_gc[i]);
  }
  XFreeGC(display, background_gc[0]);
  XFreeGC(display, background_gc[1]);
  XFreeGC(display, background_gc[15]);
}


long get_int32(unsigned char *i32)
{
  if (i32[0] & 0x80)
  {
    i32[0] ^= 0xff;
    i32[1] ^= 0xff;
    i32[2] ^= 0xff;
    i32[3] ^= 0xff;
    return (-i32[0] * 16777216 - i32[1] * 65536 - i32[2] * 256 - i32[3] - 1);
  }
  else
    return (+i32[0] * 16777216 + i32[1] * 65536 + i32[2] * 256 + i32[3]);
}


int rewind_generation(FILE *f)
{
  unsigned char b;
  long g = -1, last_fpos = -1, fpos;

  fpos = ftell(f);
  do
    b = fgetc(f);
  while (!feof(f) && !ferror(f) && (b <= 0xf0));
  if (feof(f) || ferror(f))
  {
    /* printf("rewind_generation: resetting to resume_fpos\n"); */
    fseek(f, resume_fpos, SEEK_SET);
    frame_count = resume_frame;
    return (0);
  }
  while ((b > 0xf0) && (last_fpos == -1))
  {
    switch(b)
    {
    case 0xf8:
      g = read_int32(f);
      break;
    case 0xf9:
      last_fpos = read_int32(f);
      if (last_fpos == -1)
      {
	fprintf(stderr, "Failed to rewind: Preceding position does not exist\n");
	fseek(f, fpos, SEEK_SET);
	return (-1);
      }
      break;
    case 0xfa:
    case 0xfb:
      fprintf(stderr, "Failed to rewind: No preceding position found\n");
      fseek(f, fpos, SEEK_SET);
      return (-1);
    default:
      fprintf(stderr, "Unknown reserved code %02x\n", b);
      break;
    }
    b = fgetc(f);
    if (feof(f) || ferror(f))
    {
      /* printf("rewind_generation: resetting to resume_fpos\n"); */
      fseek(f, resume_fpos, SEEK_SET);
      frame_count = resume_frame;
      return (-2);
    }
  }
  if (last_fpos == -1)
  {
    fprintf(stderr, "Failed to rewind: No preceding position found\n");
    fseek(f, fpos, SEEK_SET);
    return (-1);
  }
  /* printf("rewinding generation %ld\n", generation); */
  fseek(f, last_fpos, SEEK_SET);
  frame_count--;
  resume_fpos = fpos;
  day_count = 0;
  if (g > -1)
  {
    generation = g;
    resume_frame = frame_count;
  }
  /* printf("rewound generation %ld\n", g); */
  return (0);
}


long read_frame(FILE *f, unsigned char *buffer, unsigned char *huff)
{
  long g = -1, count, huff_length, length = 0, fpos;
  unsigned char b;

  fpos = ftell(f);
  do
    b = fgetc(f);
  while (!feof(f) && !ferror(f) && (b <= 0xf0));
  if (feof(f) || ferror(f))
  {
    return (0);
  }
  while ((b > 0xf0) && (length == 0))
  {
    switch(b)
    {
    case 0xf8:
      g = read_int32(f);
      break;
    case 0xf9:
      read_int32(f);
      break;
    case 0xfa:
      length = read_int32(f);
      if (fread(buffer, 1, length, f) < length)
      {
	/* printf("read_frame: resetting to resume_fpos\n"); */
	fseek(f, resume_fpos, SEEK_SET);
	frame_count = resume_frame;
	return (0);
      }
      day_count++;
      if ((g > -1) && (g != generation))
      {
	generation = g;
	day_count = 0;
      }
      break;
    case 0xfb:
      huff_length = read_int32(f);
      if (fread(huff, 1, huff_length, f) < huff_length)
      {
	/* printf("read_frame: resetting to resume_fpos\n"); */
	fseek(f, resume_fpos, SEEK_SET);
	frame_count = resume_frame;
	return (0);
      }
      length = huffman_decode_8bit(huff, buffer);
      day_count++;
      if ((g > -1) && (g != generation))
      {
	generation = g;
	day_count = 0;
      }
      break;
    default:
      fprintf(stderr, "Unknown reserved code %02x\n", b);
      break;
    }
    if (length == 0)
      b = fgetc(f);
    if (feof(f) || ferror(f))
    {
      /* printf("read_frame: resetting to resume_fpos\n"); */
      fseek(f, resume_fpos, SEEK_SET);
      frame_count = resume_frame;
      return (0);
    }
  }
  frame_count++;
  resume_fpos = fpos;
  resume_frame = frame_count;
  return (length);
}


int skip_frame(FILE *f)
{
  unsigned char b;
  long g = -1, data_size = -1, fpos;

  fpos = ftell(f);
  do
    b = fgetc(f);
  while (!feof(f) && !ferror(f) && (b <= 0xf0));
  if (feof(f) || ferror(f))
  {
    /* printf("skip_frame: resetting to resume_fpos\n"); */
    fseek(f, resume_fpos, SEEK_SET);
    frame_count = resume_frame;
    return (-2);
  }
  while ((b > 0xf0) && (data_size == -1))
  {
    switch(b)
    {
    case 0xf8:
      g = read_int32(f);
      break;
    case 0xf9:
      read_int32(f);
      break;
    case 0xfa:
    case 0xfb:
      data_size = read_int32(f);
      break;
    default:
      fprintf(stderr, "Unknown reserved code %02x\n", b);
      break;
    }
    if (data_size == -1)
      b = fgetc(f);
    if (feof(f) || ferror(f))
    {
      /* printf("skip_frame: resetting to resume_fpos\n"); */
      fseek(f, resume_fpos, SEEK_SET);
      frame_count = resume_frame;
      return (-2);
    }
  }
  if (data_size == -1)
  {
    fprintf(stderr, "Failed to skip: No data block size found\n");
    fseek(f, fpos, SEEK_SET);
    return (-1);
  }
  /* printf("skipping generation %ld\n", generation); */
  if (fseek(f, data_size, SEEK_CUR) == 0)
    frame_count++;
  resume_fpos = fpos;
  day_count++;
  if ((g > -1) && (g != generation))
  {
    day_count = 0;
    generation = g;
  }
  resume_frame = frame_count;
  /* printf("skipped generation %ld\n", g); */
  return (0);
}


int draw_generation(FILE *f, Drawable d, unsigned char *buffer, unsigned char *huff)
{
  long x, x0, y, count, p, num_xpoints, length;
  unsigned char b;

  length = read_frame(f, buffer, huff);
  if (length == 0)
    return (-1);
  x = 0;
  p = 0;
  num_xpoints = 0;
  y = world_height - 1;
  do
  {
    count = 1;
    b = buffer[p++];
    if (b < 0xf0)
    {
      if ((0x70 <= b) && (b <= 0x7f))
      {
	if (b == 0x70)
	  gc = background_gc[0];
	else if (b == 0x71)
	  gc = background_gc[1];
	else
	  gc = background_gc[15];
      }
      else if ((b & 0x80))
	gc = lndx_gc[(b & 0x7f) % 6];
      else
	gc = lnd_gc[b % 6];
      /* XDrawPoint(display, d, gc, x, world_height - y - 1); */
      XFillRectangle(display, d, gc, x * cell_boxsize, (world_height - y - 1) * cell_boxsize, cell_boxsize, cell_boxsize);
      x++;
      if (x == world_width)
      {
	x = 0;
	y--;
      }
    }
    else if (b == 0xf0)
    {
      count = get_int32(buffer + p);
      p += 4;
      if (p >= length)
	break;
      b = buffer[p++];
      if ((0x70 <= b) && (b <= 0x7f))
      {
	if (b == 0x70)
	  gc = background_gc[0];
	else if (b == 0x71)
	  gc = background_gc[1];
	else
	  gc = background_gc[15];
      }
      else if ((b & 0x80))
	gc = lndx_gc[(b & 0x7f) % 6];
      else
	gc = lnd_gc[b % 6];
      while (count)
      {
	x0 = x;
	x += count - 1;
	if (x >= world_width)
	  x = world_width - 1;
	count -= x - x0 + 1;
	if (x0 == x)
	  XFillRectangle(display, d, gc, x * cell_boxsize, (world_height - y - 1) * cell_boxsize, cell_boxsize, cell_boxsize);
	  /* XDrawPoint(display, d, gc, x, world_height - y - 1); */
	else
	{
	  /* printf("drawing line (%ld, %ld) -> (%ld, %ld)\n", x, world_height - y - 1, x, world_height - y - 1); */
	  /* XDrawLine(display, d, gc, x0, world_height - y - 1, x, world_height - y - 1); */
	  XFillRectangle(display, d, gc, x0 * cell_boxsize, (world_height - y - 1) * cell_boxsize, (x - x0 + 1) * cell_boxsize, cell_boxsize);
	}
	x++;
	if (x == world_width)
	{
	  x = 0;
	  y--;
	}
      }
    }
    else if (b == 0xf8)
    {
      fprintf(stderr, "Unexpected beginning of generation %02x %02x %02x %02x, x=%ld, y=%ld\n",
	      buffer[p], buffer[p + 1], buffer[p + 2], buffer[p + 3], x, y);
      break;
    }
    else
    {
      fprintf(stderr, "Unknown / unexpected reserved code %02x\n", b);
      break;
    }
    if (p >= length)
      break;
  }
  while (y > -1);
  return (0);
  /* getchar(); */
}


void shm_fillrect(XImage *shm_image, long x0, long y0, long w, long h, unsigned long pixel)
{
  long xmax = x0 + w, ymax = y0 + h, x, y;

  if ((shm_image->format == ZPixmap) && (shm_image->bits_per_pixel == 8))
  {
    for (y = y0; y < ymax; y++)
    {
      for (x = x0; x < xmax; x++)
	shm_image->data[y * shm_image->bytes_per_line + x] = pixel & 0xff;
    }
  }
  else
  {
    for (y = y0; y < ymax; y++)
    {
      for (x = x0; x < xmax; x++)
	XPutPixel(shm_image, x, y, pixel);
    }
  }
}


int shm_draw_generation(FILE *f, XImage *shm_image, unsigned char *buffer, unsigned char *huff)
{
  long x, x0, y, count, p, num_xpoints, length;
  unsigned char b;
  unsigned long pixel;

  length = read_frame(f, buffer, huff);
  if (length == 0)
    return (-1);
  x = 0;
  p = 0;
  num_xpoints = 0;
  y = world_height - 1;
  do
  {
    count = 1;
    b = buffer[p++];
    if (b < 0xf0)
    {
      if ((0x70 <= b) && (b <= 0x7f))
      {
	if (b == 0x70)
	  pixel = back_pixel[0];
	else if (b == 0x71)
	  pixel = back_pixel[1];
	else
	  pixel = back_pixel[15];
      }
      else if ((b & 0x80))
	pixel = lndx_pixel[(b & 0x7f) % 6];
      else
	pixel = lnd_pixel[b % 6];
      /* XDrawPoint(display, d, gc, x, world_height - y - 1); */
      shm_fillrect(shm_image, x * cell_boxsize, (world_height - y - 1) * cell_boxsize, cell_boxsize, cell_boxsize, pixel);
      x++;
      if (x == world_width)
      {
	x = 0;
	y--;
      }
    }
    else if (b == 0xf0)
    {
      count = get_int32(buffer + p);
      p += 4;
      if (p >= length)
	break;
      b = buffer[p++];
      if ((0x70 <= b) && (b <= 0x7f))
      {
	if (b == 0x70)
	  pixel = back_pixel[0];
	else if (b == 0x71)
	  pixel = back_pixel[1];
	else
	  pixel = back_pixel[15];
      }
      else if ((b & 0x80))
	pixel = lndx_pixel[(b & 0x7f) % 6];
      else
	pixel = lnd_pixel[b % 6];
      while (count)
      {
	x0 = x;
	x += count - 1;
	if (x >= world_width)
	  x = world_width - 1;
	count -= x - x0 + 1;
	if (x0 == x)
	  shm_fillrect(shm_image, x * cell_boxsize, (world_height - y - 1) * cell_boxsize, cell_boxsize, cell_boxsize, pixel);
	  /* XDrawPoint(display, d, gc, x, world_height - y - 1); */
	else
	{
	  /* printf("drawing line (%ld, %ld) -> (%ld, %ld)\n", x, world_height - y - 1, x, world_height - y - 1); */
	  /* XDrawLine(display, d, gc, x0, world_height - y - 1, x, world_height - y - 1); */
	  shm_fillrect(shm_image, x0 * cell_boxsize, (world_height - y - 1) * cell_boxsize, (x - x0 + 1) * cell_boxsize, cell_boxsize, pixel);
	}
	x++;
	if (x == world_width)
	{
	  x = 0;
	  y--;
	}
      }
    }
    else if (b == 0xf8)
    {
      fprintf(stderr, "Unexpected beginning of generation %02x %02x %02x %02x, x=%ld, y=%ld\n",
	      buffer[p], buffer[p + 1], buffer[p + 2], buffer[p + 3], x, y);
      break;
    }
    else
    {
      fprintf(stderr, "Unknown / unexpected reserved code %02x\n", b);
      break;
    }
    if (p >= length)
      break;
  }
  while (y > -1);
  /* printf("displayed frame %ld\n", frame_count); */
  return (0);
  /* getchar(); */
}


int main(int argc, char **argv)
{
  int num_queued;
  time_t start_time;
  double t;
  long num_frames, fast_skip = 200, fast_rewind = 200;
  unsigned int width, height;
  XEvent xevent;
  const char *prgname = argv[0];
  const char *classname = "XLndmovie";
  char *pixel_file_name = NULL;
  FILE *pixel_file;
  unsigned char *buffer = NULL, *huff = NULL;
  XImage *shm_image;
  XShmSegmentInfo shminfo;
  int shm_major, shm_minor;
  Bool shm_pixmaps;
  int ShmCompletionType;
  int shmtransfer_completed;
  int i, count, step_key = 0;
  char xkey[32];
  KeySym keysym;
  XComposeStatus compose;
  int oc;
  extern int optind;
  extern char *optarg;


  init_signal_handling();
  display = xwin_init(NULL, prgname, classname, argc, argv, &screen_no, &width, &height);
  fixed_font = XLoadQueryFont(display, "fixed");
  mit_shm = XShmQueryVersion(display, &shm_major, &shm_minor, &shm_pixmaps);
  mit_shm = XShmQueryExtension(display);
  if (mit_shm == True)
  {
    /* printf("Shared memory extension used (found V %d.%d)\n", shm_major, shm_minor); */
    ShmCompletionType = XShmGetEventBase(display) + ShmCompletion;
  }
/*
  else
    printf("Standard Xlib used, no shared memory\n");
*/
  while ((oc = getopt(argc, argv, "c:f:r:h")) != -1)
  {
    switch (oc)
    {
    case 'c':
      cell_boxsize = strtol(optarg, NULL, 10);
      break;
    case 'f':
      fast_skip = strtol(optarg, NULL, 10);
      if (fast_skip < 1)
	fast_skip = 200;
      break;
    case 'r':
      fast_rewind = strtol(optarg, NULL, 10);
      if (fast_rewind < 1)
	fast_rewind = 200;
      break;
    default:
      fprintf(stderr, "Unknown option \'-%c\' -- ignored\n", oc);
      break;
    }
  }
  if (optind < argc)
    pixel_file_name = argv[optind++];
  else
  {
    fprintf(stderr, "No pixel file specified -- exit\n");
    exit (EXIT_FAILURE);
  }
  pixel_file = fopen(pixel_file_name, "rb");
  if (pixel_file == NULL)
  {
    fprintf(stderr, "Failed to open pixel file \"%s\" -- exit\n", pixel_file_name);
    exit (EXIT_FAILURE);
  }
  read_pixheader(pixel_file);
  if ((cell_boxsize <= 0) || ((width / world_width) < cell_boxsize))
    cell_boxsize = width / world_width;
  if ((height / world_height) < cell_boxsize)
    cell_boxsize = height / world_height;
  if (cell_boxsize <= 0)
    cell_boxsize = 1;
  if ((height / world_height) < cell_boxsize)
    cell_boxsize = height / world_height;
  window_width = world_width * cell_boxsize;
  pixmap_height = world_height * cell_boxsize;
  if (fixed_font)
    window_height = pixmap_height + fixed_font->ascent + fixed_font->descent + 5;
  else
    window_height = pixmap_height;
  buffer = (unsigned char *) malloc(world_width * world_height * sizeof(char));
  if (buffer == NULL)
  {
    fprintf(stderr, "Failed to allocate internal buffer\n");
    exit (EXIT_FAILURE);
  }
  huff = (unsigned char *) malloc(world_width * world_height * sizeof(char));
  if (buffer == NULL)
  {
    fprintf(stderr, "Failed to allocate internal buffer\n");
    free(buffer);
    exit (EXIT_FAILURE);
  }
  window = create_window(MOVER | FULLER | CLOSER | RTARROW | LFARROW | UPARROW | DNARROW,
	  window_width, window_height, redraw_pixworld, close_pixwindow, 100, 100, window_width, window_height);
  XSelectInput(display, window, WH_EVENTMASK | WH_SELECTMASK);
  map_window(window);
  if (create_lndgcs(lnd_gc, lndx_gc, background_gc, &text_gc))
  {
    fprintf(stderr, "Failed to create necessary GCs (insufficient colors?)\n");
    free(buffer);
    free(huff);
    exit (EXIT_FAILURE);
  }
  if (mit_shm)
  {
    shm_image = XShmCreateImage(display, DefaultVisual(display, screen_no),
	    DefaultDepth(display, screen_no), ZPixmap, NULL, &shminfo, window_width, pixmap_height);
/*
    printf("XImage Structure:\n");
    printf("width = %d, height = %d, depth = %d\n", shm_image->width, shm_image->height, shm_image->depth);
    printf("format: ");
    switch (shm_image->format)
    {
    case XYBitmap:
      printf("XYBitmap\n");
      break;
    case XYPixmap:
      printf("XYPixmap\n");
      break;
    case ZPixmap:
      printf("ZPixmap\n");
      break;
    default:
      printf("UNKNOWN (%d)\n", shm_image->format);
      break;
    }
    printf("%d bytes per line, %d bits per pixel\n", shm_image->bytes_per_line, shm_image->bits_per_pixel);
    switch (shm_image->byte_order)
    {
    case MSBFirst:
      printf("byte order: MSBFirst\n");
      break;
    case LSBFirst:
      printf("byte order: LSBFirst\n");
      break;
    default:
      printf("byte order UNKNOWN (%d)\n", shm_image->byte_order);
      break;
    }
    switch (shm_image->bitmap_bit_order)
    {
    case MSBFirst:
      printf("bitmap bit order: MSBFirst\n");
      break;
    case LSBFirst:
      printf("bitmap bit order: LSBFirst\n");
      break;
    default:
      printf("bitmap bit order UNKNOWN (%d)\n", shm_image->bitmap_bit_order);
      break;
    }
*/
    shminfo.shmid = shmget(IPC_PRIVATE, shm_image->bytes_per_line * shm_image->height, IPC_CREAT | 0777);
    if (shminfo.shmid == -1)
    {
      fprintf(stderr, "Failed to get shared memory segment, falling back to standard Xlib\n");
      mit_shm = 0;
    }
    else
    {
      shminfo.shmaddr = shmat(shminfo.shmid, 0, 0);
      shm_image->data = shminfo.shmaddr;
      shminfo.readOnly = False;
      if (XShmAttach(display, &shminfo))
	; /* printf("Shared mem successfully attached to server\n"); */
      else
      {
	mit_shm = 0;
	XDestroyImage(shm_image);
	shmdt(shminfo.shmaddr);
	shmctl(shminfo.shmid, IPC_RMID, 0);
	printf("Shared memory released\n");
      }
    }
  }
  if (!mit_shm)
  {
    paint_pixmap = XCreatePixmap(display, window, window_width, pixmap_height, DefaultDepth(display, screen_no));
    XFillRectangle(display, paint_pixmap, background_gc[15], 0, 0, window_width, pixmap_height);
  }
  redraw_pixmap = XCreatePixmap(display, window, window_width, pixmap_height, DefaultDepth(display, screen_no));
  XFillRectangle(display, redraw_pixmap, background_gc[15], 0, 0, window_width, pixmap_height);
/*
  for (i = 0; i < 6; i++)
  {
    XPutPixel(shm_image, 0, 0, lnd_pixel[i]);
    printf("lnd_pixel[%d]=%08lx -> data=%02x\n", i, lnd_pixel[i], shm_image->data[0]);
    XPutPixel(shm_image, 0, 0, lndx_pixel[i]);
    printf("lndx_pixel[%d]=%08lx -> data=%02x\n", i, lndx_pixel[i], shm_image->data[0]);
  }
  XPutPixel(shm_image, 0, 0, back_pixel);
  printf("back_pixel=%08lx -> data=%02x\n", back_pixel, shm_image->data[0]);
*/
  do
  {
    XNextEvent(display, &xevent);
    num_queued = XEventsQueued(display, QueuedAfterFlush);
    process_xevent(&xevent);
  }
  while (!((xevent.type == Expose) && (xevent.xexpose.window == window)));
  shmtransfer_completed = 1;
  start_time = time(NULL);
  num_frames = 0;
  while (!game_over)
  {
    if(feof(pixel_file) || ferror(pixel_file))
    {
      if (resume_fpos != -1)
      {
	mode = SINGLE_STEP;
	/* printf("switched to single step due to feof / ferror\n"); */
	step_key = 0;
	/*
	printf("rewinding to resume pos %ld\n", resume_fpos);
	fseek(pixel_file, resume_fpos, SEEK_SET);
	*/
      }
      else
      {
	fprintf(stderr, "Encountered EOF without finding any entry point -- quitting\n");
	game_over = 1;
      }
    }
    num_queued = XEventsQueued(display, QueuedAfterFlush);
    while ((num_queued) || ((mode == SINGLE_STEP) && !step_key))
    {
      XNextEvent(display, &xevent);
      num_queued = XEventsQueued(display, QueuedAfterFlush);
      /* printf("received event: %s, window: %d\n", event_name[xevent.type], (int) xevent.xany.window); */
      if (process_xevent(&xevent))
	continue;
      if (xevent.type == ShmCompletionType)
      {
        /* printf("shm transfer completed\n"); */
	shmtransfer_completed = 1;
      }
      if (xevent.type == KeyPress)
      {
        count = XLookupString(&(xevent.xkey), xkey, 32, &keysym, &compose);
	if (count)
	{
	  switch(xkey[0])
	  {
	  case 'f':
	    /* printf("fast forward\n"); */
	    /* printf("rewind position: %ld\n", resume_fpos); */
	    mode = FAST_FORWARD;
	    break;
	  case 'r':
	    /* printf("fast rewind\n"); */
	    /* printf("rewind position: %ld\n", resume_fpos); */
	    mode = FAST_REWIND;
	    break;
	  case 'b':
	    /* printf("backward play\n"); */
	    /* printf("rewind position: %ld\n", resume_fpos); */
	    mode = BACKWARD_PLAY;
	    rewind_generation(pixel_file);
	    break;
	  case 's':
	    /* printf("single step\n"); */
	    mode = SINGLE_STEP;
	    step_key = 0;
	    break;
	  case 'n':
	    /* printf("normal play\n"); */
	    /* printf("rewind position: %ld\n", resume_fpos); */
	    mode = NORMAL_PLAY;
	    break;
	  case 'q': case 'Q':
	    /* printf("game over, bye bye\n"); */
            mode = NORMAL_PLAY;
            step_key = 1;
	    game_over = 1;
	    break;
	  default:
	    step_key = 1;
	    break;
	  }
	}
      }
    }
    if ((mode == NORMAL_PLAY) || (mode == BACKWARD_PLAY) || ((mode == SINGLE_STEP) && step_key)
            || ((mode == FAST_FORWARD) && ((generation % fast_skip) == (fast_skip - 1)))
            || ((mode == FAST_REWIND) && ((generation % fast_rewind) == 1)))
    {
      if (mit_shm && shmtransfer_completed)
      {
	if (shm_draw_generation(pixel_file, shm_image, buffer, huff))
	{
	  mode = SINGLE_STEP;
	  /* printf("switched to single step after shm_draw_generation\n"); */
	  step_key = 0;
	}
	else
	{
	  if (window != BadValue)
	    XShmPutImage(display, window, background_gc[15], shm_image, 0, 0, 0, window_height - pixmap_height, window_width, pixmap_height, True);
	  if (redraw_pixmap != BadValue)
	  XShmPutImage(display, redraw_pixmap, background_gc[15], shm_image, 0, 0, 0, 0, window_width, pixmap_height, True);
	  shmtransfer_completed = 0;
	  num_frames++;
	}
      }
      else
      {
	if (draw_generation(pixel_file, paint_pixmap, buffer, huff))
	{
	  mode = SINGLE_STEP;
	  /* printf("switched to single step after draw_generation\n"); */
	  step_key = 0;
	}
	else
	{
	  if (window != BadValue)
	    XCopyArea(display, paint_pixmap, window, background_gc[15], 0, 0, window_width, pixmap_height, 0, window_height - pixmap_height);
	  if (redraw_pixmap != BadValue)
	    XCopyArea(display, paint_pixmap, redraw_pixmap, background_gc[15], 0, 0, window_width, pixmap_height, 0, 0);
	  num_frames++;
	}
      }
      redraw_text();
      if (mode == FAST_REWIND)
      {
	if (rewind_generation(pixel_file))
	{
	  mode = SINGLE_STEP;
	  /* printf("switched to single step after rewind_generation\n"); */
	  step_key = 0;
	}
      }
      step_key = 0;
    }
    switch (mode)
    {
    case BACKWARD_PLAY:
      if (rewind_generation(pixel_file))
      {
	mode = SINGLE_STEP;
	/* printf("switched to single step after rewind_generation\n"); */
	step_key = 0;
      }
      /* fall through */
    case FAST_REWIND:
      if (rewind_generation(pixel_file))
      {
	mode = SINGLE_STEP;
	/* printf("switched to single step after rewind_generation\n"); */
	step_key = 0;
      }
      redraw_text();
      break;
    case FAST_FORWARD:
      if (skip_frame(pixel_file))
      {
	mode = SINGLE_STEP;
	/* printf("switched to single step after skip_frame\n"); */
	step_key = 0;
      }
      redraw_text();
      break;
    }
  }
  t = difftime(time(NULL), start_time);
  /* printf("%ld frames in %f seconds (%f frames/sec)\n", num_frames, t, num_frames / t); */
  if (mit_shm)
  {
    XShmDetach(display, &shminfo);
    XDestroyImage(shm_image);
    shmdt(shminfo.shmaddr);
    shmctl(shminfo.shmid, IPC_RMID, 0);
    printf("Shared memory released\n");
  }
  free_lndgcs();
  if (window != BadValue)
    remove_window(window);
  if (redraw_pixmap != BadValue)
    XFreePixmap(display, redraw_pixmap);
  if (paint_pixmap != BadValue)
    XFreePixmap(display, paint_pixmap);
  free(buffer);
  free(huff);
  return (EXIT_SUCCESS);
}

