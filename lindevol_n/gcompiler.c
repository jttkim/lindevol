#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>


#define GENE_LENGTH   3
#define MAX_SLEN    256
#define skipspaces(s) while (isspace(*(s))) (s)++


unsigned long line_no = 0;


int process_gene_line(char *gene_hex, const char *gene_source,
        long lnd_divide, long lnd_statebit, long lnd_broadcastbit,
        long lnd_flyingseed, long lnd_localseed,
        long lnd_mutminus, long lnd_mutplus)
{
  const char *src_pos;
  unsigned char state, mask;
  int site, i, n = 0;

  if ((src_pos = strchr(gene_source, ':')) == NULL)
  {
    src_pos = gene_source;
  }
  else
  {
    src_pos += 1;
  }
  skipspaces(src_pos);
  state = 0;
  mask = 0xff;
  for (i = 7; i >= 0; i--)
  {
    if (*src_pos == '*')
    {
      mask ^= 1 << i;
    }
    else if (*src_pos == '1')
    {
      state |= 1 << i;
    }
    else if (*src_pos != '0')
    {
      fprintf(stderr, "line %lu: error: unexpected character \'%c\' in input part\n", line_no, *src_pos);
      return (-1);
    }
    src_pos++;
  }
  while (isspace(*src_pos) || (*src_pos == '-') || (*src_pos == '>'))
  {
    src_pos++;
  }
  if (strstr(src_pos, "divide") == src_pos)
  {
    if (lnd_divide < 0)
    {
      fprintf(stderr, "line %lu: error: lnd_divide not set\n", line_no);
      return (-1);
    }
    src_pos += strlen("divide");
    skipspaces(src_pos);
    if (!isdigit(*src_pos))
    {
      fprintf(stderr, "line %lu: error: divide position argument is not numeric\n", line_no);
      return (-1);
    }
    site = strtol(src_pos, (char **) NULL, 10);
    if ((site < 0) || (site > 7))
    {
      fprintf(stderr, "line %lu: error: divide position argument must be in [0, 7]\n", line_no);
      return (-1);
    }
    sprintf(gene_hex, "%02x%02x%02x", state, mask, lnd_divide + site);
    return (0);
  }
  else if (strstr(src_pos, "flying seed") == src_pos)
  {
    if (lnd_flyingseed < 0)
    {
      fprintf(stderr, "line %lu: error: lnd_flyingseed not set\n", line_no);
      return (-1);
    }
    sprintf(gene_hex, "%02x%02x%02x", state, mask, lnd_flyingseed);
    return (0);
  }
  else if (strstr(src_pos, "local seed") == src_pos)
  {
    if (lnd_localseed < 0)
    {
      fprintf(stderr, "line %lu: error: lnd_localseed not set\n", line_no);
      return (-1);
    }
    sprintf(gene_hex, "%02x%02x%02x", state, mask, lnd_localseed);
    return (0);
  }
  else if (strstr(src_pos, "mut-") == src_pos)
  {
    if (lnd_mutminus < 0)
    {
      fprintf(stderr, "line %lu: error: lnd_mutminus not set\n", line_no);
      return (-1);
    }
    sprintf(gene_hex, "%02x%02x%02x", state, mask, lnd_mutminus);
    return (0);
  }
  else if (strstr(src_pos, "mut+") == src_pos)
  {
    if (lnd_mutplus < 0)
    {
      fprintf(stderr, "line %lu: error: lnd_mutplus not set\n", line_no);
      return (-1);
    }
    sprintf(gene_hex, "%02x%02x%02x", state, mask, lnd_mutplus);
    return (0);
  }
  else if ((*src_pos  == '.') || (*src_pos == '1'))
  {
    if (lnd_statebit < 0)
    {
      fprintf(stderr, "line %lu: error: lnd_statebit not set\n", line_no);
      return (-1);
    }
    n = -1;
    for (i = 0; i < 8; i++)
    {
      if (src_pos[i] == '1')
      {
        if (n == -1)
        {
            n = 7 - i;
        }
        else
        {
          fprintf(stderr, "line %lu: error: multiple 1's in state bit specification\n", line_no);
          return (-1);
        }
      }
      else if (src_pos[i] != '.')
      {
        fprintf(stderr, "line %lu: error: unknown character \'%c\' in state bit specification\n",
                line_no, src_pos[i]);
        return (-1);
      }
    }
    if ((src_pos[8] == '.') || (src_pos[8] == '1'))
    {
      fprintf(stderr, "line %lu: warning: extraneous characters in state bit specification\n", line_no);
    }
    if (n == -1)
    {
      fprintf(stderr, "line %lu: error: no \'1\' in state bit specification\n", line_no);
      return (-1);
    }
    sprintf(gene_hex, "%02x%02x%02x", state, mask, lnd_statebit + n);
    return (0);
  }
  else if ((*src_pos  == '\'') || (*src_pos == 'b'))
  {
    if (lnd_broadcastbit < 0)
    {
      fprintf(stderr, "line %lu: error: lnd_broadcastbit not set\n", line_no);
      return (-1);
    }
    n = -1;
    for (i = 0; i < 8; i++)
    {
      if (src_pos[i] == 'b')
      {
        if (n == -1)
        {
            n = 7 - i;
        }
        else
        {
          fprintf(stderr, "line %lu: error: multiple b's in broadcast bit specification\n", line_no);
          return (-1);
        }
      }
      else if (src_pos[i] != '\'')
      {
        fprintf(stderr, "line %lu: error: unknown character \'%c\' in broadcast bit specification\n",
                line_no, src_pos[i]);
        return (-1);
      }
    }
    if ((src_pos[8] == '\'') || (src_pos[8] == 'b'))
    {
      fprintf(stderr, "line %lu: warning: extraneous characters in broadcast bit specification\n", line_no);
    }
    if (n == -1)
    {
      fprintf(stderr, "line %lu: error: no \'b\' in broadcast bit specification\n", line_no);
      return (-1);
    }
    sprintf(gene_hex, "%02x%02x%02x", state, mask, lnd_broadcastbit + n);
    return (0);
  }
  else
  {
    fprintf(stderr, "line %lu: error: unknown output part\'%s\'\n", line_no, src_pos);
    return (-1);
  }
}


int process_lnd_setting(char *line, const char *name, long *lndpar, int *errcode)
{
  if (strncmp(line, name, strlen(name)))
  {
    return (0);
  }
  line += strlen(name);
  while ((isspace(*line)) || (*line == '='))
  {
    line++;
  }
  if (*lndpar >= 0)
  {
    fprintf(stderr, "line %lu: error: %s is already set to %ld\n", line_no, name, *lndpar);
    *errcode = -1;
  }
  if (!isdigit(*line))
  {
    fprintf(stderr, "line %lu: error: invalid value specified for %s\n", line_no, name);
    *errcode = -1;
  }
  *lndpar = strtol(line, (char **) NULL, 10);
  return (1);
}


int process_file(FILE *outfile, FILE *infile)
{
  long lnd_divide = -1, lnd_statebit = -1, lnd_broadcastbit = -1,
       lnd_flyingseed = -1, lnd_localseed = -1,
       lnd_mutminus = -1, lnd_mutplus = -1;
  long genome_length = 0, n, i;
  char *genome_hex = NULL, *gh;
  int errcode = 0;
  char line[MAX_SLEN], gene_hex[GENE_LENGTH * 2 + 1];
  char *line_pos;

  while (!feof(infile) && !errcode)
  {
    if (fgets(line, MAX_SLEN, infile) == NULL)
    {
      break;
    }
    line_no++;
    if (*line == '#')
    {
      continue;
    }
    if (*line == ':')
    {
      line_pos = line + 1;
      if (!strncmp(line_pos, "write", 5))
      {
        if (genome_length == 0)
        {
          fprintf(stderr, "line %lu: error: no genome to write\n", line_no);
          continue;
        }
        line_pos += 5;
        while ((isspace(*line_pos)) || (*line_pos == '='))
        {
          line_pos++;
        }
        n = strtol(line_pos, (char **) NULL, 10);
        if (n <= 0)
        {
          fprintf(stderr, "line %lu: error: invalid argument to :write\n", line_no);
        }
        for (i = 0; i < n; i++)
        {
          fprintf(outfile, "%ld\n", genome_length);
          fprintf(outfile, "%s\n", genome_hex);
        }
        free(genome_hex);
        genome_hex = NULL;
        genome_length = 0;
        continue;
      }
      if (process_lnd_setting(line_pos, "lnd_divide", &lnd_divide, &errcode)) continue;
      if (process_lnd_setting(line_pos, "lnd_statebit", &lnd_statebit, &errcode)) continue;
      if (process_lnd_setting(line_pos, "lnd_broadcastbit", &lnd_broadcastbit, &errcode)) continue;
      if (process_lnd_setting(line_pos, "lnd_flyingseed", &lnd_flyingseed, &errcode)) continue;
      if (process_lnd_setting(line_pos, "lnd_localseed", &lnd_localseed, &errcode)) continue;
      if (process_lnd_setting(line_pos, "lnd_mutminus", &lnd_mutminus, &errcode)) continue;
      if (process_lnd_setting(line_pos, "lnd_mutplus", &lnd_mutplus, &errcode)) continue;
      fprintf(stderr, "line %lu: error: unknown keyword\n", line_no);
    }
    else
    {
      if (process_gene_line(gene_hex, line, lnd_divide, lnd_statebit, lnd_broadcastbit,
                  lnd_flyingseed, lnd_localseed, lnd_mutminus, lnd_mutplus) < 0)
      {
        errcode = -1;
        continue;
      }
      else
      {
        if (genome_length  == 0)
        {
          genome_hex = (char *) malloc(GENE_LENGTH * 2 + 1);
          if (genome_hex == NULL)
          {
            fprintf(stderr, "line %lu: out of memory error\n", line_no);
            errcode = -2;
            break;
          }
        }
        else
        {
          gh = (char *) realloc((void *) genome_hex, (genome_length + 1) * GENE_LENGTH * 2 + 1);
          if (gh)
          {
            genome_hex = gh;
          }
          else
          {
            fprintf(stderr, "line %lu: out of memory error\n", line_no);
            errcode = -2;
            break;
          }
        }
      }
      sprintf(genome_hex + genome_length * GENE_LENGTH * 2, "%s", gene_hex);
      genome_length++;
    }
  }
  if (genome_length)
  {
    free(genome_hex);
  }
  return (errcode);
}


int main(int argc, char **argv)
{
  FILE *infile, *outfile;
  int errcode;

  if (argc > 1)
  {
    if ((infile = fopen(argv[1], "r")) == NULL)
    {
      fprintf(stderr, "failed to open \"%s\" for input -- exit\n", argv[1]);
      return (EXIT_FAILURE);
    }
  }
  else
  {
    infile = stdin;
  }
  if (argc > 2)
  {
    if ((outfile = fopen(argv[2], "w")) == NULL)
    {
      fprintf(stderr, "failed to open \"%s\" for output -- exit\n", argv[2]);
      return (EXIT_FAILURE);
    }
  }
  else
  {
    outfile = stdout;
  }
  line_no = 0;
  errcode = process_file(outfile, infile);
  if (infile != stdin)
  {
    fclose(infile);
  }
  if (outfile != stdout)
  {
    fclose(outfile);
  }
  if (errcode < 0)
  {
    return (EXIT_FAILURE);
  }
  return (EXIT_SUCCESS);
}

