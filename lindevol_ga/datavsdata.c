#include <stdio.h>
#include <stdlib.h>

char buf[256];
char pro_fname[256], out_fname[256], ofn[256];
char *datname[256];
FILE *pro_file, *out_file;

int get_num_numdata(FILE *f)
{
  int num_numdata = 0;
  int i, j, n;

  fgets(buf, 256, f);
  fgets(buf, 256, f);
  fgets(buf, 256, f);
  fgets(buf, 256, f);
  n = atoi(buf);
  for (i = 0; i < n; i++)
    fgets(buf, 256, f);
  fgets(buf, 256, f);
  n = atoi(buf);
  for (i = 0; i < n; i++)
  {
    fgets(buf, 256, f);
    if (buf[0] == 'n')
    {
      fgets(buf, 256, f);
      for (j = strlen(buf) - 1; buf[j] == '\r' || buf[j] == '\n'; j--)
        buf[j] = '\0';
      datname[num_numdata] = (char *) malloc(strlen(buf) + 1);
      strcpy(datname[num_numdata], buf);
      num_numdata++;
    }
    else
    {
      fgets(buf, 256, f);
      fgets(buf, 256, f);
    }
  }
  return(num_numdata);
}


void make_outfile(int num_numdata, int nx, int ny, FILE *pf, FILE *f)
{
  auto long g = 0;
  long i, j;
  char xval[256], yval[256];

  fprintf(f, "# plot of %s vs. %s\n", datname[nx], datname[ny]);
  for (i = 0; fgets(buf, 256, pf) != NULL; i++)
  {
    if ((i % num_numdata) == nx)
    {
      for (j = strlen(buf) - 1; buf[j] == '\r' || buf[j] == '\n'; j--)
        buf[j] = '\0';
      strcpy(xval, buf);
    }
    if ((i % num_numdata) == ny)
    {
      for (j = strlen(buf) - 1; buf[j] == '\r' || buf[j] == '\n'; j--)
        buf[j] = '\0';
      strcpy(yval, buf);
    }
    if (i % num_numdata == num_numdata - 1)
    {
      fprintf(f, "%s %s\n", xval, yval);
      g++;
    }
  }
}


int main(int argc, char *argv[])
{
  int  num_numdata;
  int  i, j;
  int  nx = -1, ny = -1;

  if (argc > 1)
    strcpy(pro_fname, argv[1]);
  else
  {
    printf("name of protocol file: ");
    fgets(pro_fname, 256, stdin);
    for (j = strlen(pro_fname) - 1; pro_fname[j] == '\r' || pro_fname[j] == '\n'; j--)
      pro_fname[j] = '\0';
  }
  if (argc > 2)
    strcpy(out_fname, argv[2]);
  else
  {
    printf("name of output files: ");
    fgets(out_fname, 256, stdin);
    for (j = strlen(out_fname) - 1; out_fname[j] == '\r' || out_fname[j] == '\n'; j--)
      out_fname[j] = '\0';
  }
  if (argc > 4)
  {
    nx = strtol(argv[3], (char **) NULL, 10);
    ny = strtol(argv[4], (char **) NULL, 10);
  }
  if ((pro_file = fopen(pro_fname, "r")) == (FILE *) NULL)
  {
    printf("unable to open protocol file \"%s\"\n", pro_fname);
    return (-1);
  }
  num_numdata = get_num_numdata(pro_file);
  if (nx == -1)
  {
    for (i = 0; i < num_numdata; i++)
    {
      printf("%3d: %s\n", i, datname[i]);
    }
    printf("\nData to be plotted on X axis: ");
    fgets(buf, 256, stdin);
    nx = strtol(buf, (char **) NULL, 10);
    printf("Data to be plotted on X axis: ");
    fgets(buf, 256, stdin);
    ny = strtol(buf, (char **) NULL, 10);
  }
  sprintf(ofn, "%s.%03d.%03d", out_fname, nx, ny);
  if ((out_file = fopen(ofn, "w")) == (FILE *) NULL)
  {
    printf("unable to open output file \"%s\"\n", ofn);
    return (-2);
  }
  make_outfile(num_numdata, nx, ny, pro_file, out_file);
  fclose(out_file);
  fclose(pro_file);
  return(0);
}

