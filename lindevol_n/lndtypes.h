#ifndef H_LNDTYPES
#define H_LNDTYPES

#include "lndvals.h"

#include "genome.h"
#include "gntypes.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  long divide_code;
  long flyingseed_code;
  long localseed_code;
  long mutminus_code;
  long mutplus_code;
  long to_epool_code;
  long to_npool_code;
  long from_epool_code;
  long from_npool_code;
  long statebit_code;
  long num_divide;
  long num_flyingseed;
  long num_localseed;
  long num_mutminus;
  long num_mutplus;
  long num_to_epool;
  long num_to_npool;
  long num_from_epool;
  long num_from_npool;
  long num_statebit;
} GSYS_PARAMETERS;

typedef struct
{
  long plant_no;          /* number of the plant that occupies that site, -1 for empty */
  long cell_no;           /* number of the cell in the plant */
  long organic_nutrient;
  unsigned char nutrient;
} LATTICE_SITE;

typedef struct
{
  long x;                        /* coordinates of */
  long y;                        /* the cell */
  unsigned char energy;          /* state of the cell (energy-rich or -less) */
  unsigned char nutrient;        /* nutrient state of cell */
  unsigned long state, next_state;
} PL_CELL;

typedef struct
{
  GENOME      genome;            /* genome of the plant */
  GN_NODE_ID  node;              /* node ID in genealogy tree */
  long        num_cells;         /* number of cells of the plant */
  PL_CELL    *cell;              /* the cells of the plant */
  long        cellular_energy;   /* number of energyrich cells */
  long        cellular_nutrient; /* number of nutrient containing cells */
  long        energy_pool;
  long        nutrient_pool;
  long        age;               /* age of plant in generations */
} PLANT;

typedef struct
{
  long num;   /* number of genomes of that species */
  long index; /* index of an individual of that species */
} SPECIES_D;  /* species descriptor */

typedef struct
{
  unsigned long valid_bits, state_bits;
} STATE_SPEC;

typedef struct
{
  STATE_SPEC statespec;
  int        output;
} GENE_SPEC;

#ifdef __cplusplus
}
#endif

#endif

