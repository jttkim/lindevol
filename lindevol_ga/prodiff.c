#include <stdio.h>
#include <stdlib.h>
#include <time.h>


char   buf[256];
char   datname[50][256];

long count_numentries(FILE *f)
{
  long i, j, x;
  long n = 0;

  fgets(buf, 256, f);
  fgets(buf, 256, f);
  fgets(buf, 256, f);
  fgets(buf, 256, f);
  x = strtol(buf, (char **) NULL, 10);
  for (i = 0; i < x; i++)
    fgets(buf, 256, f);
  fgets(buf, 256, f);
  x = strtol(buf, (char **) NULL, 10);
  for (i = 0; i < x; i++)
  {
    fgets(buf, 256, f);
    fgets(datname[n], 256, f);
    for (j = strlen(datname[n]) - 1; (datname[n][j] == '\n') || (datname[n][j] == '\r'); j--)
      datname[n][j] = '\0';
    if (buf[0] == 'n')
    {
      n++;
    }
    else
    {
      fgets(buf, 256, f);
    }
  }
  return (n);
}

void write_header(FILE *f, long n, const char *fn1, const char *fn2)
{
  long i;
  time_t t;

  t = time(NULL);
  fprintf(f, "Protocol Diff V. 0.0\n");
  fprintf(f, "%s", ctime(&t));
  fprintf(f, "%s", ctime(&t));
  fprintf(f, "1\n");
  fprintf(f, "Diff between %s and %s\n", fn1, fn2);
  fprintf(f, "%ld\n", n);
  for (i = 0; i < n; i++)
  {
    fprintf(f, "n\n");
    fprintf(f, "Diffs of %s\n", datname[i]);
  }
}

int main(int argc, char *argv[])
{
  FILE *i1, *i2, *o;
  long n1, n2;
  double x1, x2;

  if (argc != 4)
  {
    printf("Usage: prodiff <infile1> <infile2> <outfile>\n");
    return (-1);
  }
  if ((i1 = fopen(argv[1], "r")) == NULL)
  {
    printf("failed to open %s\n", argv[1]);
    return (-2);
  }
  if ((i2 = fopen(argv[2], "r")) == NULL)
  {
    printf("failed to open %s\n", argv[2]);
    return (-2);
  }
  if ((o = fopen(argv[3], "w")) == NULL)
  {
    printf("failed to open %s\n", argv[3]);
    return (-2);
  }
  n1 = count_numentries(i1);
  n2 = count_numentries(i2);
  if (n1 != n2)
  {
    printf("%s (n=%ld) and %s (n=%ld): incompatible protocol files\n", argv[1], n1, argv[2], n2);
    return (-3);
  }
  write_header(o, n1, argv[1], argv[2]);
  while ((!feof(i1)) && (!feof(i2)))
  {
    fgets(buf, 256, i1);
    x1 = strtod(buf, (char **) NULL);
    fgets(buf, 256, i2);
    x2 = strtod(buf, (char **) NULL);
    if ((!feof(i1)) && (!feof(i2)))
      fprintf(o, "%f\n", x1 - x2);
  }
  fclose(i1);
  fclose(i2);
  fclose(o);
  return (0);
}

