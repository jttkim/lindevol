#include <stdio.h>

#include "lndvals.h"
#include "lndtypes.h"
#include "lndglobl.h"
#include "lndlibin.h"
#include "lndlib.h"


void prepare_filenames(void)
{
  sprintf(pro_file_name, "%s.pro", simname);
  sprintf(dmt_file_name, "%s.dmt", simname);
  sprintf(jf_file_name, "%s.jf", simname);
  sprintf(genome_file_name, "%s.gen", simname);
  sprintf(dst_file_name, "%s.dst", simname);
  sprintf(bpe_file_name, "%s.bpe", simname);
  sprintf(save_file_name, "%s.dat", simname);
  sprintf(pixel_file_name, "%s.pix", simname);
}

