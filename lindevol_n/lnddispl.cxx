#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include <qapp.h>
#include <qwidget.h>
#include "lndwidgt.h"
#include "ctrlwdgt.h"

#include "lnderror.h"
#include "lndvals.h"
#include "lndtypes.h"

#include "lndglobl.h"
#include "lnddispl.h"


static long curr_generation;
static QApplication *q_app = NULL;
static LndWorldScrollWidget *lnd_widget = NULL;
static ControlWidget *control_widget = NULL;

FILE  *world_file = NULL;
char   world_file_name[MAX_SLEN];


static void show_world(long generation, FILE *f)
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


extern "C" {

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


void display_start(long generation)
{
  curr_generation = generation;
}


void display_done(void)
{
  if ((worldmode != 0) && (world_file != NULL))
    show_world(curr_generation, world_file);
  poll_user_interface();
}


void display_plant_killed(long plant_no)
{
  poll_user_interface();
}


void display_plant_mutated(long plant_no)
{
  poll_user_interface();
}


void display_cell(long x, long y)
{
  if (lnd_widget)
  {
    lnd_widget->world_widget->show_cell(x, y);
    poll_user_interface();
  }
}


void poll_user_interface(void)
{
  if (q_app)
    q_app->processEvents();
}


void init_world_display(const char *simname)
{
  static int dummy_argc;
  static char *dummy_argv[] = { SIMPRGNAME, "-fn", "gem16", "" };

  if (worldmode)
  {
    open_world_file(simname, "a");
  }
  if (quietmode)
    return;
  for (dummy_argc = 0; dummy_argv[dummy_argc][0]; dummy_argc++)
    ;
  q_app = new QApplication(dummy_argc, dummy_argv);
  lnd_widget = new LndWorldScrollWidget();
  q_app->setMainWidget(lnd_widget);
  lnd_widget->show();
  control_widget = new ControlWidget();
  control_widget->resize(200, 100);
  control_widget->show();
  poll_user_interface();
}


void close_world_display(void)
{
  if (worldmode)
  {
    close_world_file();
  }
  if (q_app)
  {
    poll_user_interface();
    delete q_app;
    q_app = NULL;
  }
}

} /* end of extern "C" */

