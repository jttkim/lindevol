/* a module for porting stuff from GFA basic to C. It implements
   some GFA-like functions. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TMPSTRLEN 256
static char tmpstr[TMPSTRLEN];

/* bas_int emulates the typical basic INT() function */

long bas_int(double x)
{
  double i;

  (void) modf(x, &i);
  return ((long) i);
}


/* boolean requester. Writes a requester string to stdout and
   reads input from stdin. If the first character of the input
   is 'Y' or 'y', TRUE is returned, else FALSE is returned. */

short input_boolean(const char *txt)
{
  fprintf(stdout, "%s", txt);
  fprintf(stdout, " (y/n)? ");
  fgets(tmpstr, TMPSTRLEN, stdin);
  return (((tmpstr[0] == 'Y') || (tmpstr[0] == 'y')));
}

/* input_long() emulates the basic input command for signed longs */

signed long input_long(const char *txt)
{
  fprintf(stdout, "%s", txt);
  fgets(tmpstr, TMPSTRLEN, stdin);
  return (strtol(tmpstr, (char **) NULL, 10));
}

/* input_ul() emulates the basic input command for unsigned longs */

unsigned long input_ul(const char *txt)
{
  fprintf(stdout, "%s", txt);
  fgets(tmpstr, TMPSTRLEN, stdin);
  return (atol(tmpstr));
}

/* input_dbl() emulates the basic input command for doubles */

double input_dbl(const char *txt)
{
  fprintf(stdout, "%s", txt);
  fgets(tmpstr, TMPSTRLEN, stdin);
  return (strtod(tmpstr, (char **) NULL));
}

/* input_str() requests a string from stdin. A null-terminated
   string without any trailing \r's or \n's is written into buf. */

char *input_str(const char *txt, char *buf, unsigned long max_len)
{
  unsigned long i;

  fprintf(stdout, "%s", txt);
  fgets(buf, max_len, stdin);

  for (i = strlen(buf) - 1; (buf[i] == '\r') || (buf[i] == '\n'); i--)
  {
    buf[i] = '\0';
    if (i == 0)
    {
      break;
    }
  }
  return (buf);
}

double rnd(void)
{
  return(((double) (random() & 0x7fffff)) / ((double) 0x800000));
}

