#include <stdio.h>

#include "lndvals.h"
#include "lndtypes.h"
#include "lndglobl.h"
#include "lndlibin.h"
#include "lndlib.h"


void close_data_files(void)
{
  if (pro_file != NULL)
    fclose(pro_file);
  if (dmt_file != NULL)
    fclose(dmt_file);
  if (jf_file != NULL)
    fclose(jf_file);
  if (dst_file != NULL)
    fclose(dst_file);
  if (bpe_file != NULL)
    fclose(bpe_file);
  if (pixel_file != NULL)
    fclose(pixel_file);
}

