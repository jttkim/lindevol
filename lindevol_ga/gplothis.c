#include <stdio.h>
#include <stdlib.h>

short c_histogram[2000], histogram[2000];
char  buf[256], h_fname[256], o_fname[256];

int expand_histogram(int l, const short *c, short *h)
{
  auto int i = 0;
  auto int y = 0;
  int j;

  while (i < l)
  {
    if (c[i] < 0)
    {
      for (j = 0; j < -c[i]; j++)
        h[y++] = c[i + 1];
      i += 2;
    }
    else
      h[y++] = c[i++];
  }
  return(y);
}


int main(int argc, char *argv[])
{
  long g = 0;
  long max_hlen = 0;
  int  y;
  int xstep = 5;
  int ystep = 2;
  int  i, j, l;
  FILE *h, *o;

  if (argc > 1)
    strcpy(h_fname, argv[1]);
  else
  {
    printf("name of histogram file: ");
    fgets(h_fname, 256, stdin);
    for (j = strlen(h_fname) - 1; h_fname[j] == '\r' || h_fname[j] == '\n'; j--)
      h_fname[j] = '\0';
  }
  if (argc > 2)
    strcpy(o_fname, argv[2]);
  else
  {
    printf("name of output file: ");
    fgets(o_fname, 256, stdin);
    for (j = strlen(o_fname) - 1; o_fname[j] == '\r' || o_fname[j] == '\n'; j--)
      o_fname[j] = '\0';
  }
  if (argc > 3)
    xstep = atoi(argv[3]);
  if (argc > 4)
    ystep = atoi(argv[4]);
  if ((h = fopen(h_fname, "rb")) == (FILE *) NULL)
  {
    printf("unable to open histogram file \"%s\"\n", h_fname);
    return (-1);
  }
  if ((o = fopen(o_fname, "w")) == (FILE *) NULL)
  {
    printf("unable to open output file \"%s\"\n", o_fname);
    return (-1);
  }
  while (fgets(buf, 256, h) != NULL)
  {
    l = atoi(buf);
    fread(c_histogram, sizeof(short), l, h);
    if ((g % xstep) == 0)
    {
      l = expand_histogram(l, c_histogram, histogram);
      max_hlen = (l > max_hlen) ? l : max_hlen;
    }
  }
  rewind(h);
  g = 0;
  while (fgets(buf, 256, h) != NULL)
  {
    l = atoi(buf);
    fread(c_histogram, sizeof(short), l, h);
    if ((g % xstep) == 0)
    {
      l = expand_histogram(l, c_histogram, histogram);
      for (y = 0; y < l; y += ystep)
        fprintf(o, "%ld %d %d\n", g, y, histogram[y]);
      for (y = l; y < max_hlen; y += ystep)
        fprintf(o, "%ld %d 0\n", g, y);
      fprintf(o, "\n");
    }
    g++;
  }
  fclose(h);
  fclose(o);
  return (0);
}

