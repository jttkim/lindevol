/* $Id: pgrow.h,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: pgrow.h,v $
 * Revision 1.2  2000/01/30 03:11:00  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

#ifndef H_PGROW
#define H_PGROW

#ifdef __cplusplus
extern "C" {
#endif

extern long create_cell(long p, long x, long y);
extern long divide(long plant_no, long cell_no, int direction);
extern void cell_activity(long plant_no, long cell_no, long action, unsigned char gene_output);
extern void plant_growth(long plant_no);

#ifdef __cplusplus
}
#endif

#endif /* H_PGROW */

