#include "lndvals.h"
#include "lndtypes.h"
#include "lndglobl.h"
#include "lnderror.h"
#include "lndlib.h"
#include "lnd6.h"


/*
 * indicate whether an action code encodes an energy dependent
 * action
 */

int is_energy_dependent(long action_code)
{
  return ((action_code == LND_DIVIDE) || (action_code == LND_FLYINGSEED)
          || (action_code == LND_LOCALSEED) || (action_code == LND_MUTMINUS)
          || (action_code == LND_MUTPLUS));
}


/*
 * Determine the actual state of a cell
 */

unsigned char cell_state(long plant_no, long cell_no)
{
  unsigned char cs = 0;
  long x, y;
  int i;

  for (i = 0; i < 8; i++)
  {
    x = (plant[plant_no]->cell[cell_no].x + x_offset[i] + world_width) % world_width;
    y = plant[plant_no]->cell[cell_no].y + y_offset[i];
    if ((y >= 0) && (y < world_height))
    {
      if (world[x][y].plant_no == plant_no)
      {
        cs |= (1 << i);
      }
    }
  }
  return (cs);
}


/*
 * Compute the state specification encoded by the regulatory
 * part of a gene
 */

long compute_statespec(const GENOME *genome, long gpos, STATE_SPEC *statespec, long *num_bitchecks)
{
  long p = gpos;
  unsigned long regulator_bit;
  int regulator_bitno;

  statespec->valid_bits = 0;
  statespec->state_bits = 0;
  while ((p < genome->length) && !l4_promoter(genome->g[p]))
  {
    if (l4_terminator(genome->g[p]))
      break;
    regulator_bitno = genome->g[p] & 0xf;
    regulator_bit = 1UL << regulator_bitno;
    statespec->valid_bits |= regulator_bit;
    if (genome->g[p] & 0x20)
      statespec->state_bits |= regulator_bit;
    if (num_bitchecks)
    {
      num_bitchecks[regulator_bitno]++;
      /* printf("bit #%lu was checked\n", regulator_bitno); */
    }
    p++;
  }
  return (p);
}


/*
 * The state of a cell is determined by the pattern of cells
 * of the same plant in the moore neighborhood of the cell
 * in question.
 */

long process_cell(long plant_no, long cell_no, unsigned char *p_gene_output, const GENE_SPEC *gene_spec)
{
  unsigned char output;
  unsigned long cs;
  unsigned long gpos = 0, gene_no = 0;
  STATE_SPEC statespec;
  long ga, action = -1;
  GENOME *genome = &(plant[plant_no]->genome);

  cs = cell_state(plant_no, cell_no);
  /* printf("cell state: %2x\n", cs); */
  cs |= plant[plant_no]->cell[cell_no].state;
  if (plant[plant_no]->cell[cell_no].energy)
    cs |= LND_ENERGYBIT;
  if (plant[plant_no]->cell[cell_no].nutrient)
    cs |= LND_NUTRIENTBIT;
  if (plant[plant_no]->cell[cell_no].y <= world_soil)
    cs |= LND_INSOILBIT;
  plant[plant_no]->cell[cell_no].next_state = 0;
  if (gene_spec)
  {
    for (gene_no = 0; gene_no < plant[plant_no]->genome.num_genes; gene_no++)
    {
      if (((cs ^ gene_spec[gene_no].statespec.state_bits) & gene_spec[gene_no].statespec.valid_bits) == 0)
      {
        if (gene_spec[gene_no].output >= 0)
        {
          output = gene_spec[gene_no].output;
          ga = gene_activity(output);
          if (is_energy_dependent(ga))
          {
            if (action == -1)
            {
              action = ga;
              *p_gene_output = output;
              if (genome->flags & GNM_USG)
                genome->usg_count[gene_no]++;
              if ((genome->flags & GNM_BP) && (genome->bp_count[gene_no] < bp_ceiling))
                genome->bp_count[gene_no]++;
            }
          }
          else
          {
            switch (ga)
            {
            case LND_STATEBIT:
              plant[plant_no]->cell[cell_no].next_state |= 1 << (output & 0xf);
              num_bitsets[output & 0xf]++;
              if (genome->flags & GNM_USG)
                genome->usg_count[gene_no]++;
              if ((genome->flags & GNM_BP) && (genome->bp_count[gene_no] < bp_ceiling))
                genome->bp_count[gene_no]++;
              break;
            case LND_TO_EPOOL:
              if (plant[plant_no]->cell[cell_no].energy)
              {
                plant[plant_no]->cell[cell_no].energy = 0;
                plant[plant_no]->cellular_energy--;
                plant[plant_no]->energy_pool++;
              }
              num_to_epool++;
              break;
            case LND_TO_NPOOL:
              if ((plant[plant_no]->cell[cell_no].nutrient) && (plant[plant_no]->nutrient_pool < plant[plant_no]->num_cells))
              {
                plant[plant_no]->cell[cell_no].nutrient = 0;
                plant[plant_no]->cellular_nutrient--;
                plant[plant_no]->nutrient_pool++;
              }
              num_to_npool++;
              break;
            case LND_FROM_EPOOL:
              if (plant[plant_no]->energy_pool > 0)
              {
                plant[plant_no]->energy_pool--;
                if (!plant[plant_no]->cell[cell_no].energy)
                {
                  plant[plant_no]->cell[cell_no].energy = 1;
                  plant[plant_no]->cellular_energy++;
                }
              }
              num_from_epool++;
              break;
            case LND_FROM_NPOOL:
              if ((plant[plant_no]->nutrient_pool > 0) && (!plant[plant_no]->cell[cell_no].nutrient))
              {
                plant[plant_no]->nutrient_pool--;
                plant[plant_no]->cell[cell_no].nutrient = 1;
                plant[plant_no]->cellular_nutrient++;
              }
              num_from_npool++;
              break;
            default:
              fprintf(stderr, "process_cell: unknown action code %ld\n", ga);
              break;
            }
          }
        }
      }
    }
  }
  else
  {
    while (gpos < genome->length)
    {
      while ((gpos < genome->length) && !l4_promoter(genome->g[gpos]))
        gpos++;
      gpos++;
      if (gpos >= genome->length)
        break;
      if (cell_no)  /* HACK: do num_bitcheck statistic only once per plant, at cell #0 */
        gpos = compute_statespec(genome, gpos, &statespec, NULL);
      else
        gpos = compute_statespec(genome, gpos, &statespec, num_bitchecks);
      if (gpos >= genome->length)
        break;
      if (l4_terminator(genome->g[gpos]))
      {
        if (((cs ^ statespec.state_bits) & statespec.valid_bits) == 0)
        {
          output = genome->g[gpos];
          ga = gene_activity(output);
          if (is_energy_dependent(ga))
          {
            if (action == -1)
            {
              action = ga;
              *p_gene_output = output;
              if (genome->flags & GNM_USG)
                genome->usg_count[gene_no]++;
              if ((genome->flags & GNM_BP) && (genome->bp_count[gene_no] < bp_ceiling))
                genome->bp_count[gene_no]++;
            }
          }
          else
          {
            switch (ga)
            {
            case LND_STATEBIT:
              plant[plant_no]->cell[cell_no].next_state |= 1 << (output & 0xf);
              num_bitsets[output & 0xf]++;
              if (genome->flags & GNM_USG)
                genome->usg_count[gene_no]++;
              if ((genome->flags & GNM_BP) && (genome->bp_count[gene_no] < bp_ceiling))
                genome->bp_count[gene_no]++;
              break;
            case LND_TO_EPOOL:
              if (plant[plant_no]->cell[cell_no].energy)
              {
                plant[plant_no]->cell[cell_no].energy = 0;
                plant[plant_no]->cellular_energy--;
                plant[plant_no]->energy_pool++;
              }
              num_to_epool++;
              break;
            case LND_TO_NPOOL:
              if ((plant[plant_no]->cell[cell_no].nutrient) && (plant[plant_no]->nutrient_pool < plant[plant_no]->num_cells))
              {
                plant[plant_no]->cell[cell_no].nutrient = 0;
                plant[plant_no]->cellular_nutrient--;
                plant[plant_no]->nutrient_pool++;
              }
              num_to_npool++;
              break;
            case LND_FROM_EPOOL:
              if (plant[plant_no]->energy_pool > 0)
              {
                plant[plant_no]->energy_pool--;
                plant[plant_no]->cell[cell_no].energy = 1;
                plant[plant_no]->cellular_energy++;
              }
              num_from_epool++;
              break;
            case LND_FROM_NPOOL:
              if ((plant[plant_no]->nutrient_pool > 0) && (!plant[plant_no]->cell[cell_no].nutrient))
              {
                plant[plant_no]->nutrient_pool--;
                plant[plant_no]->cell[cell_no].nutrient = 1;
                plant[plant_no]->cellular_nutrient++;
              }
              num_from_npool++;
              break;
            default:
              fprintf(stderr, "process_cell: unknown action code %ld\n", ga);
              break;
            }
          }
        }
      }
      gene_no++;
    }
  }
  return (action);
}

