/* $Id: geneact6.c,v 1.2 2000/01/30 03:10:59 kim Exp $ */
/*
 * $Log: geneact6.c,v $
 * Revision 1.2  2000/01/30 03:10:59  kim
 * Added cvs tags
 * Switched to urandom dependent lndrandm (this should be moved to a lib)
 * Added nutrient flux: free nutrient may diffuse out of the world and is
 *     generated at random locations. New control parameters:
 *     * nutrient_per_timestep
 *     * organic_nutrient_diffusion
 *
 */

#include <limits.h>
#include <stddef.h>
#include <stdlib.h>

#include "lndvals.h"
#include "lndtypes.h"
#include "lnderror.h"
#include "lndlibin.h"
#include "lndlib.h"
#include "lnd6.h"
#include "lndglobl.h"


int l4_promoter(unsigned char c)
{
  return ((c & 0x80) == 0x80);
}


int l4_promotermask(unsigned char c)
{
  return (c == 0x80);
}


int l4_terminator(unsigned char c)
{
  return ((c & 0x40) == 0x40);
}


unsigned long l4_num_genes(const GENOME *genome)
{
  unsigned long i, n = 0;

  for (i = 0; i < genome->length; i++)
  {
    if (l4_promoter(genome->g[i]))
      n++;
  }
  return (n);
}


unsigned long l4_gene_number(const GENOME *genome, unsigned long pos)
{
  unsigned long i, n = 0;

  if (pos >= genome->length)
    return (genome->num_genes);
  for (i = 0; i <= pos; i++)
  {
    if (l4_promoter(genome->g[i]))
      n++;
  }
  if (n == 0)
    return (ULONG_MAX);
  return (n - 1);
}


unsigned long l4_gene_startposition(const GENOME *genome, unsigned long g_num)
{
  unsigned long i, n = 0;

  for (i = 0; i < genome->length; i++)
  {
    if (l4_promoter(genome->g[i]))
    {
      if (n == g_num)
        return (i);
      n++;
    }
  }
  return (ULONG_MAX);
}


int calculate_activity_codes(GSYS_PARAMETERS *gp)
{
  long sum = gp->num_divide + gp->num_flyingseed + gp->num_localseed
             + gp->num_mutminus + gp->num_mutplus + gp->num_statebit
	     + gp->num_to_epool + gp->num_to_npool
	     + gp->num_from_epool + gp->num_from_npool;

  if (sum != 64L)
  {
    fprintf(stderr, "sum of numbers of code bytes does not amount to 64\n");
    gp->num_divide = gp->num_divide * 64 / sum;
    gp->num_flyingseed = gp->num_flyingseed * 64 / sum;
    gp->num_localseed = gp->num_localseed * 64 / sum;
    gp->num_mutminus = gp->num_mutminus * 64 / sum;
    gp->num_mutplus = gp->num_mutplus * 64 / sum;
    gp->num_to_epool = gp->num_to_epool * 64 / sum;
    gp->num_to_npool = gp->num_to_npool * 64 / sum;
    gp->num_from_epool = gp->num_from_epool * 64 / sum;
    gp->num_from_npool = gp->num_from_npool * 64 / sum;
    gp->num_statebit = 64 - gp->num_divide - gp->num_flyingseed - gp->num_localseed
            - gp->num_mutminus - gp->num_mutplus - gp->num_to_epool - gp->num_to_npool
	    - gp->num_from_epool - gp->num_from_npool;
  }
  gp->divide_code = 0;
  gp->flyingseed_code = gp->divide_code + gp->num_divide;
  gp->localseed_code = gp->flyingseed_code + gp->num_flyingseed;
  gp->mutminus_code = gp->localseed_code + gp->num_localseed;
  gp->mutplus_code = gp->mutminus_code + gp->num_mutminus;
  gp->to_epool_code = gp->mutplus_code + gp->num_mutplus;
  gp->to_npool_code = gp->to_epool_code + gp->num_to_epool;
  gp->from_epool_code = gp->to_npool_code + gp->num_to_npool;
  gp->from_npool_code = gp->from_epool_code + gp->num_from_epool;
  gp->statebit_code = gp->from_npool_code + gp->num_from_npool;
  return (0);
}


/*
 * return the gene activity corresponding to the output passed
 * in argument.
 */

long gene_activity(unsigned char gene_output)
{
  gene_output &= 0x3f;
  if (gsys_parameters.num_divide && (gsys_parameters.divide_code <= gene_output) && (gene_output < gsys_parameters.flyingseed_code))
    return (LND_DIVIDE);
  if (gsys_parameters.num_flyingseed && (gsys_parameters.flyingseed_code <= gene_output) && (gene_output < gsys_parameters.localseed_code))
    return (LND_FLYINGSEED);
  if (gsys_parameters.num_localseed && (gsys_parameters.localseed_code <= gene_output) && (gene_output < gsys_parameters.mutminus_code))
    return (LND_LOCALSEED);
  if (gsys_parameters.num_mutminus && (gsys_parameters.mutminus_code <= gene_output) && (gene_output < gsys_parameters.mutplus_code))
    return (LND_MUTMINUS);
  if (gsys_parameters.num_mutplus && (gsys_parameters.mutplus_code <= gene_output) && (gene_output < gsys_parameters.to_epool_code))
    return (LND_MUTPLUS);
  if (gsys_parameters.num_to_epool && (gsys_parameters.to_epool_code <= gene_output) && (gene_output < gsys_parameters.to_npool_code))
    return (LND_TO_EPOOL);
  if (gsys_parameters.num_to_npool && (gsys_parameters.to_npool_code <= gene_output) && (gene_output < gsys_parameters.from_epool_code))
    return (LND_TO_NPOOL);
  if (gsys_parameters.num_from_epool && (gsys_parameters.from_epool_code <= gene_output) && (gene_output < gsys_parameters.from_npool_code))
    return (LND_FROM_EPOOL);
  if (gsys_parameters.num_from_npool && (gsys_parameters.from_npool_code <= gene_output) && (gene_output < gsys_parameters.statebit_code))
    return (LND_FROM_NPOOL);
  /* if (gsys_parameters.statebit_code <= gene_output) */
  return (LND_STATEBIT);
}

