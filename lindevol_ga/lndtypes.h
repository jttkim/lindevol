#ifndef H_LNDTYPES
#define H_LNDTYPES

#include "lndvals.h"

#include "genome.h"

#include "gntypes.h"


typedef struct
{
  long divide_code;
  long statebit_code;
  long broadcastbit_code;
  long flyingseed_code;
  long localseed_code;
  long mutminus_code;
  long mutplus_code;
  long num_divide;
  long num_statebit;
  long num_broadcastbit;
  long num_flyingseed;
  long num_localseed;
  long num_mutminus;
  long num_mutplus;
} GSYS_PARAMETERS;

typedef struct
{
  GENOME genome;
  unsigned long fitness;                  /* fitness value of the genome */
  GN_NODE_ID    node_id;                  /* node ID for genealogy package */
} LND_GENOME;

typedef struct
{
  long plant_no;          /* number of the plant that occupies that site, -1 for empty */
  unsigned long cell_no;  /* number of the cell in the plant */
} LATTICE_SITE;

typedef struct
{
  long x;               /* coordinates of */
  long y;               /* the cell */
  char          energy; /* state of the cell (energy-rich or -less) */
} PL_CELL;

typedef struct
{
  unsigned long num_cells;     /* number of cells of the plant */
  PL_CELL cell[MAX_NUM_CELLS]; /* the cells of the plant */
} PLANT;

typedef struct
{
  unsigned long num;   /* number of genomes of that species */
  unsigned long index; /* index of an individual of that species */
} SPECIES_D;           /* species descriptor */

typedef struct
{
  unsigned char valid_bits, state_bits;
} STATE_SPEC;

#endif


