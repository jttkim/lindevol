/* $Id: lnderror.c,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: lnderror.c,v $
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

#include "lndvals.h"
#include "lndlib.h"

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif


void do_error(char *message)
{
  printf("%s\n", message);
}

