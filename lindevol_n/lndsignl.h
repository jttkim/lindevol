/* $Id: lndsignl.h,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: lndsignl.h,v $
 * Revision 1.2  2000/01/30 03:11:00  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

#ifndef H_LNDSIGNL
#define H_LNDSIGNL

#ifdef __cplusplus
extern "C" {
#endif

extern void init_signal_handling(void);

#ifdef __cplusplus
}
#endif

#endif

