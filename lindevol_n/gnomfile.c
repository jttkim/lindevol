#include <stdio.h>

#include "genomelib.h"

#include "lndvals.h"
#include "lndtypes.h"
#include "lndglobl.h"
#include "lndlibin.h"
#include "lndlib.h"
#include "gnlib.h"


FILE *genome_file = NULL;


int open_genome_file(const char *mode)
{
  genome_file = fopen(genome_file_name, mode);
  if (genome_file == NULL)
  {
    perror("open_genome_file failed");
    return (-1);
  }
  return (0);
}


void write_genomes(long num_samples, const long *sample_index)
{
  long plant_no, i;

  fprintf(genome_file, "g %ld\n", generation);
  fprintf(genome_file, "%ld\n", num_samples);
  for (i = 0; i < num_samples; i++)
  {
    plant_no = sample_index[i];
    if (plant[plant_no])
    {
      gn_save_id(&(plant[plant_no]->node), genome_file);
      write_genome(genome_file, &(plant[plant_no]->genome), 0);
/*
      write_pirseq(genome_file, plant_no);
      for (j = 0; j < plant->genome.length; j++)
      {
        fprintf(genome_file, "%02x", plant->genome.g[j]);
      }
      fprintf(genome_file, "\n");
*/
    }
    else
      fprintf(stderr, "write_genomes: Erroneous sample_index[%ld] = %ld", i, sample_index[i]);
  }
}

