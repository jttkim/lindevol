/* $Id: prepfnam.c,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: prepfnam.c,v $
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
#include "lndlibin.h"
#include "lndlib.h"


void prepare_filenames(void)
{
  sprintf(pro_file_name, "%s.pro", simname);
  sprintf(dmt_file_name, "%s.dmt", simname);
  sprintf(jf_file_name, "%s.jf", simname);
  sprintf(genome_file_name, "%s.gen", simname);
  sprintf(dst_file_name, "%s.dst", simname);
  sprintf(bpe_file_name, "%s.bpe", simname);
  sprintf(save_file_name, "%s.dat", simname);
  sprintf(pixel_file_name, "%s.pix", simname);
}

