#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>


#define BUFLEN 256


void do_mean(char fname[], int s)
{
  int p;
  FILE *f;
  char buf[BUFLEN];
  char *tmp, *final;
  double x;
  auto double sum = 0.0;
  auto long n = 0;

  if ((f = fopen(fname, "r")) != (FILE *) NULL)
  {
    while (!feof(f))
    {
      fgets(buf, BUFLEN, f);
      if (buf[0] != '#')
      {
        tmp = buf;
        for (p = 0; p < s; p++)
        {
          (void) strtod(tmp, &final);
        tmp = final;
        }
        errno = 0;
        x = strtod(tmp, &final);
        if ((!errno) && (tmp < final))
        {
          sum += x;
          n++;
        }
        else
        {
          fprintf(stderr, "mean: no numerical value: %s\n", buf);
        }
      }
    }
    fclose(f);
    if (n > 0)
    {
      fprintf(stdout, "%s: %f\n", fname, sum / n);
    }
/*
    sprintf(buf, "%s.mean", fname);
    if ((f = fopen(buf, "w")) != (FILE *) NULL)
    {
      fprintf(f, "%lf\n", sum / n);
      fclose(f);
    }
    else
    {
      fprintf(stderr, "mean: could not create %s\n", buf);
    }
*/
  }
  else
  {
    fprintf(stderr, "mean: failed to open %s\n", fname);
  }
}


int main(int argc, char *argv[])
{
  int s = 0;
  int a;

  for (a = 1; a < argc; a++)
  {
    if (argv[a][0] == '-')
    {
      if (isdigit((int) (argv[a][1])))
      {
        s = strtol(&(argv[a][1]), (char **) NULL, 10) - 1;
      }
    }
    else
    {
      do_mean(argv[a], s);
    }
  }
  return (0);
}

