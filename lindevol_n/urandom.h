/* $Id: urandom.h,v 1.1 2000/01/30 03:11:00 kim Exp $ */
/*
 * $Log: urandom.h,v $
 * Revision 1.1  2000/01/30 03:11:00  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

#ifndef URANDOM_H
#define URANDOM_H

extern void ulong_srandom( unsigned int x);
extern unsigned long *ulong_initstate(unsigned int seed, unsigned long *arg_state);
extern unsigned long *ulong_setstate(unsigned long *arg_state);
extern long int ulong_random(void);
extern unsigned long urandom_long(unsigned long range);
extern double urandom_double(void);
extern double urandom_gauss(void);
extern int write_urandom_state(FILE *f);
extern int read_urandom_state(FILE *f);
extern void *urandom_shuffle(size_t num, size_t s, void *a);

#endif

