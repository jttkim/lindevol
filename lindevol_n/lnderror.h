/* $Id: lnderror.h,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: lnderror.h,v $
 * Revision 1.2  2000/01/30 03:11:00  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

#ifndef H_LNDERROR
#define H_LNDERROR

#ifdef __cplusplus
extern "C" {
#endif

extern void do_error(char *message);

#ifdef __cplusplus
}
#endif

#endif


