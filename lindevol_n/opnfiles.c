/* $Id: opnfiles.c,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: opnfiles.c,v $
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
#include "lndtypes.h"
#include "lndglobl.h"
#include "lnddispl.h"
#include "lndlib.h"


int open_data_files(const char *mode)
{
  prepare_filenames();
  if (open_pro_file(mode) != 0)
  {
    clear_arrays();
    return (-1);
  }
  if (dmt_savefreq > 0)
    if (open_dmt_file(mode) != 0)
      dmt_savefreq = 0;
  if (phyltest_savefreq > 0)
  {
    if (open_jf_file(mode) != 0)
      phyltest_savefreq = 0;
    if (open_genome_file(mode) != 0)
      phyltest_savefreq = 0;
  }
  if (ddistr_savefreq > 0)
    if (open_dst_file(mode) != 0)
      ddistr_savefreq = 0;
  if (bp_ceiling > 0)
    if (open_bpe_file(mode) != 0)
      bp_ceiling = 0;
  if (worldmode)
    open_world_file(simname, mode);
  if (pixmode)
    if (open_pixfile(mode) != 0)
      pixmode = 0;
  return (0);
}

