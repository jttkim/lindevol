/* $Id: lndvals.h,v 1.2 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: lndvals.h,v $
 * Revision 1.2  2000/01/30 03:11:00  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

#ifndef H_LNDVALS
#define H_LNDVALS

#define SIMPRGNAME "LindEvol-6.1"

#define MAX_SLEN 256            /* max. length of strings */

#define NUM_CELL_ACTIONS  9

#define NUM_STATEBITS     16
#define LND_ENERGYBIT     0x8000  /* 1 << (NUM_STATEBITS - 1) */
#define LND_NUTRIENTBIT   0x4000  /* 1 << (NUM_STATEBITS - 2) */
#define LND_INSOILBIT     0x2000  /* 1 << (NUM_STATEBITS - 3) */

#define LND_DIVIDE        1
#define LND_STATEBIT      2
#define LND_BROADCASTBIT  3
#define LND_FLYINGSEED    4
#define LND_LOCALSEED     5
#define LND_MUTMINUS      6
#define LND_MUTPLUS       7
#define LND_TO_EPOOL      8
#define LND_TO_NPOOL      9
#define LND_FROM_EPOOL   10
#define LND_FROM_NPOOL   11


#define PR_BIRTH          0
#define PR_DEATH          1

#define GNM_USG           1
#define GNM_BP            2

#endif

