/* $Id: bpe_file.c,v 1.2 2000/01/30 03:10:59 kim Exp $ */
/*
 * $Log: bpe_file.c,v $
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
#include <stdlib.h>

#include "jklib.h"
#include "lndvals.h"
#include "lndtypes.h"
#include "lndglobl.h"
#include "lnderror.h"
#include "lndlibin.h"
#include "lndlib.h"


FILE *bpe_file = NULL;


int open_bpe_file(const char *mode)
{
  sprintf(buf, "%sb", mode);
  bpe_file = fopen(bpe_file_name, buf);
  if (bpe_file == NULL)
  {
    perror("open_bpe_file failed");
    return (-1);
  }
  return (0);
}


/*
PROCEDURE bedau_and_packard
  '
  ' **************************************************************************
  ' * determines the Bedau and Packard evolutionary activity in a generation *
  ' * and saves it.                                                          *
  ' **************************************************************************
*/

void write_bpe(void)
{
  long i, j;
  long compress_len;
  long bp_max = 0;

  short *bp, *compress;

  bp = (short *) malloc((bp_ceiling + 1) * sizeof(short));
  if (bp == NULL)
  {
    do_error("write_bpe failed: couldn't allocate memory");
  }
  else
  {
    compress = (short *) malloc((bp_ceiling + 1) * sizeof(short));
    if (compress == NULL)
    {
      do_error("write_bpe failed: couldn't allocate memory");
    }
    else
    {
      /* printf("write_bpe: starting\n"); */
      j = 0;
      for (i = 0; i < world_width; i++)
      {
        if (plant[i] != NULL)
        {
          tmp_index[j++] = i;
        }
      }
      for (i = 0; i <= bp_ceiling; i++)
      {
        bp[i] = 0;
      }
      for (i = 0; i < psize; i++)
      {
        /* printf("write_bpe: processing genome %ld\n", i); */
        /* printf("write_bpe: number of genes in genome: %ld\n", genome[i].num_genes); */
        for (j = 0; j < plant[tmp_index[i]]->genome.num_genes; j++)
        {
          bp_max = (plant[tmp_index[i]]->genome.bp_count[j] > bp_max) ? plant[tmp_index[i]]->genome.bp_count[j] : bp_max;
          bp[plant[tmp_index[i]]->genome.bp_count[j]]++;
        }
      }
      /* printf("write_bpe: loop done\n"); */
      compress_len = compress_histogram(bp_max + 1, bp, compress);
      /* printf("write_bpe: histogram compressed\n"); */
      fprintf(bpe_file, "%ld\r\n", compress_len);
      fwrite_int16array(compress, (size_t) compress_len, bpe_file);
      fflush(bpe_file);
      /* printf("write_bpe: %ld elements (%ld bytes) written\n", compress_len, compress_len * sizeof(compress[0])); */
      free((void *) compress);
    }
    free((void *) bp);
  }
}

