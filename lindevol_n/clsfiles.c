/* $Id: clsfiles.c,v 1.2 2000/01/30 03:10:59 kim Exp $ */
/*
 * $Log: clsfiles.c,v $
 * Revision 1.2  2000/01/30 03:10:59  kim
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


void close_data_files(void)
{
  if (pro_file != NULL)
    fclose(pro_file);
  if (dmt_file != NULL)
    fclose(dmt_file);
  if (jf_file != NULL)
    fclose(jf_file);
  if (dst_file != NULL)
    fclose(dst_file);
  if (bpe_file != NULL)
    fclose(bpe_file);
  if (pixel_file != NULL)
    fclose(pixel_file);
}

