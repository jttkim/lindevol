#include <stdio.h>
#include <stdlib.h>

#include "jklib.h"
#include "lndvals.h"
#include "lndtypes.h"
#include "lndglobl.h"
#include "lnderror.h"
#include "lndlibin.h"
#include "lndlib.h"


FILE *dst_file = NULL;


int open_dst_file(const char *mode)
{
  sprintf(buf, "%sb", mode);
  dst_file = fopen(dst_file_name, buf);
  if (dst_file == NULL)
  {
    perror("open_dst_file failed");
    return (-1);
  }
  if (*mode == 'w')
  {
    fprintf(dst_file, "ProPlot-1.0 INT32\r\n");
    fprintf(dst_file, "i%1ld\r\n", ddistr_savefreq);
  }
  return (0);
}


/*
PROCEDURE distance_distribution
  ' *******************************************************************
  ' * Saves the distribution of distance values in the matrix of edit *
  ' * distances between all genomes in the population.                *
  ' *******************************************************************
  write_it controls whether the distance distribution is actually written
  out or not. If no distance distribution is written, the side effect of
  computing distance entropies is still needed.
*/

void write_dst(char write_it, long sample_size, long *sample_index)
{
  long i, j;
  long compress_len;
  long d;
  long max_distance = 0;
  long max_len, max_genome_length;
  long dd_rel[100];

  long *dd, *compress;

  max_genome_length = 0;
  for (i = 0; i < sample_size; i++)
  {
    if (plant[sample_index[i]] && (plant[sample_index[i]]->genome.length > max_genome_length))
      max_genome_length = plant[sample_index[i]]->genome.length;
  }
  dd = (long *) malloc((max_genome_length + 1) * sizeof(long));
  if (dd == NULL)
  {
    do_error("write_dst failed: couldn't allocate memory");
  }
  else
  {
    compress = (long *) malloc((max_genome_length + 1) * sizeof(long));
    if (compress == NULL)
    {
      do_error("write_dst failed: couldn't allocate memory");
    }
    else
    {
      /* printf("write_dst: starting\n"); */
      for (i = 0; i <= max_genome_length; i++)
        dd[i] = 0;
      /* printf("write_dst: dd[] initialized to zero\n"); */
      for (i = 0; i < 100; i++)
        dd_rel[i] = 0;
      for (i = 0; i < sample_size - 1; i++)
      {
        for (j = i + 1; j < sample_size; j++)
        {
          d = (long) edit_distance(plant[sample_index[i]]->genome.length, (char *) plant[sample_index[i]]->genome.g,
                                    plant[sample_index[j]]->genome.length, (char *) plant[sample_index[j]]->genome.g);
          /* printf("d(%ld, %ld) = %ld\n", i, j, d); */
          if (d >= 0)
          {
            dd[d]++;
            max_distance = (d > max_distance) ? d : max_distance;
            max_len = (plant[sample_index[i]]->genome.length > plant[sample_index[j]]->genome.length ) ?
                       plant[sample_index[i]]->genome.length : plant[sample_index[j]]->genome.length;
            if (max_len > 0)
              d = d * 100 / max_len;
            else
              d = 0;
            dd_rel[(d < 100) ? d : 99]++;
          }
        }
      }
      /* printf("write_dst: dd[] filled, max_distance=%ld\n", max_distance); */
      if (write_it)
      {
        compress_len = compress_histogram_long(max_genome_length + 1, dd, compress);
        /* printf("write_dst: dd[] compressed, length now: %ld\n", compress_len); */
        fprintf(dst_file, "%ld\r\n", compress_len);
        fwrite_int32array(compress, compress_len, dst_file);
        fflush(dst_file);
      }
      distance_entropy = shannon_long(max_distance, dd);
      distance_entropy_rel = shannon_long(100, dd_rel);
      /* printf("s = %lf, rel. s = %lf\n", distance_entropy, distance_entropy_rel); */
      free((void *) compress);
    }
    free((void *) dd);
  }
}

