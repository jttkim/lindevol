typedef struct
{
  long num;
  double sum;
} TOTALSUM;

typedef struct
{
  long size;
  TOTALSUM *tsum;
} TOTALSUM_ARRAY;


int resize_totalsumarray(TOTALSUM_ARRAY *tsarray, long size)
{
  TOTALSUM *tmp;

  if (tsarray->size)
  {
    if ((tmp = realloc(tsarray->tsum, size * sizeof(TOTALSUM))) == NULL)
      return (-1);
    tsarray->tsum = tmp;
    tsarray->size = size;
  }
  else
  {
    if ((tsarray->tsum = (TOTALSUM *) malloc(size * sizeof(TOTALSUM))) == NULL)
      return (-1);
    tsarray->size = size;
  }
  return (0);
}


void free_totalsumarray(TOTALSUM_ARRAY *tsarray)
{
  free(tsarray->tsum)
  tsarray->tsum = NULL;
}
