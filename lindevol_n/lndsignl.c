#include <stdio.h>
#include <signal.h>

#include "lndvals.h"
#include "lndglobl.h"
#include "lndlibin.h"

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif


void handle_signal(int sig_id)
{
  switch (sig_id)
  {
  case SIGTERM: case SIGINT:
    signal(sig_id, handle_signal);
    fprintf(stderr, "%s will terminate after finishing generation #%1ld\n", SIMPRGNAME, generation);
    finish_flag = 1;
    break;
  case SIGHUP:
    signal(SIGHUP, handle_signal);
    fprintf(stderr, "%s will save after finishing generation #%1ld\n", SIMPRGNAME, generation);
    save_data = 1;
    break;
  case SIGUSR1:
    signal(SIGUSR1, handle_signal);
    quietmode = !quietmode;
    if (quietmode)
    {
      fprintf(stderr, "online output off\n");
    }
    else
    {
      fprintf(stderr, "online output on\n");
    }
    break;
  default:
    fprintf(stderr, "Caught signal %d but don't know what to do\n", sig_id);
    break;
  }
}


void init_signal_handling(void)
{
  signal(SIGTERM, handle_signal);
  signal(SIGINT, handle_signal);
  signal(SIGHUP, handle_signal);
  signal(SIGUSR1, handle_signal);
}

