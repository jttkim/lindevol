/* $Id: timefile.c,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: timefile.c,v $
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

#include "lndvals.h"
#include "lndtypes.h"
#include "lnderror.h"
#include "lndglobl.h"
#include "lndlibin.h"
#include "lndlib.h"


FILE *savetime_file = NULL;
char savetime_file_name[MAX_SLEN] = "";


int open_savetime_file(const char *fname)
{
  strncpy(savetime_file_name, fname, MAX_SLEN);
  if ((savetime_file = fopen(savetime_file_name, "r")) == NULL)
    return (-1);
  return (0);
}


void close_savetime_file(void)
{
  if (savetime_file)
    fclose(savetime_file);
  savetime_file = NULL;
}


long savetime_next(long generation)
{
  char buf[MAX_SLEN], e[MAX_SLEN * 2], *s;
  size_t i;
  long n;

  if (savetime_file == NULL)
    return (-1);
  for (;;)
  {
    do
      fgets(buf, MAX_SLEN, savetime_file);
    while (!feof(savetime_file) && !ferror(savetime_file) && (buf[0] == '\n'));
    if (feof(savetime_file) || ferror(savetime_file))
      break;
    i = strlen(buf);
    if (i && (buf[i - 1] == '\n'))
      buf[i - 1] = '\0';
    n = strtol(buf, &s, 10);
    if (buf == s)
    {
      sprintf(e, "savetime_next: Illegal value \"%s\"", buf);
      do_error(e);
    }
    else if (n <= generation)
    {
      sprintf(e, "savetime_next: Skipping \"%s\", already at generation %ld", buf, generation);
      do_error(e);
    }
    else
      return (n);
  }
  close_savetime_file();
  return (-1);
}

