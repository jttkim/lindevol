#include <stdio.h>

#include "lndvals.h"
#include "lndlib.h"

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif


void do_error(char *message)
{
  printf("%s\n", message);
}

