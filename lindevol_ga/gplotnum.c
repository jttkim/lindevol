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


void make_outfile(int num_numdata, int n, FILE *pf, FILE *f)
{
  auto long g = 0;
  long i, j;

  fprintf(f, "# %s\n", datname[n]);
  for (i = 0; fgets(buf, 256, pf) != NULL; i++)
  {
    if ((i % num_numdata) == n)
    {
      for (j = strlen(buf) - 1; buf[j] == '\r' || buf[j] == '\n'; j--)
        buf[j] = '\0';
      fprintf(f, "%d %s\n", g++, buf);
    }
  }
}


void make_viewfile(int num_numdata)
{
  int i;
  FILE *f;

  if ((f = fopen("viewme.gnp", "w")) != (FILE *) NULL)
  {
    for (i = 0; i < num_numdata; i++)
    {
      fprintf(f, "plot \'%s.%03d\' t \'%s\' with lines\n", out_fname, i, datname[i]);
      fprintf(f, "pause -1 \"Bash the return key to continue\"\n");
    }
    fclose(f);
  }
}


int main(int argc, char *argv[])
{
  long datstart;
  int  num_numdata;
  int  i, j;

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
  if ((pro_file = fopen(pro_fname, "r")) == (FILE *) NULL)
  {
    printf("unable to open protocol file \"%s\"\n", pro_fname);
    return (-1);
  }
  num_numdata = get_num_numdata(pro_file);
  datstart = ftell(pro_file);
  for (i = 0; i < num_numdata; i++)
  {
    sprintf(ofn, "%s.%03d", out_fname, i);
    fseek(pro_file, datstart, 0);
    if ((out_file = fopen(ofn, "w")) == (FILE *) NULL)
    {
      printf("unable to open output file \"%s\"\n", ofn);
      return (-2);
    }
    make_outfile(num_numdata, i, pro_file, out_file);
    fclose(out_file);
  }
  fclose(pro_file);
  make_viewfile(num_numdata);
  return(0);
}

