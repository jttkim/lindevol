#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

#include "genomelib.h"
#include "lndtypes.h"
#include "lndlib.h"
#include "lnderror.h"
#include "lnd6.h"


#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif


static void delete_counter(GENOME *genome, unsigned long gene_no)
{
  long i;

  if (genome->flags & GNM_USG)
  {
    for (i = gene_no + 1; i < genome->num_genes; i++)
      genome->usg_count[i - 1] = genome->usg_count[i];
    if (gene_no > 0)
      genome->usg_count[gene_no - 1] = 0;
  }
  if (genome->flags & GNM_BP)
  {
    for (i = gene_no; i < genome->num_genes; i++)
      genome->bp_count[i - 1] = genome->bp_count[i];
    if (gene_no > 0)
      genome->bp_count[gene_no - 1] = 0;
  }
}


static void insert_counter(GENOME *genome, unsigned long gene_no)
{
  long i, n;

#ifdef MEMDEBUG
  verify_Watchdogs();
#endif
  if (gene_no == ULONG_MAX)
    n = 0;
  else
    n = gene_no + 1;
  if (genome->flags & GNM_USG)
  {
    for (i = genome->num_genes - 1; n < i; i--)
      genome->usg_count[i] = genome->usg_count[i - 1];
    if (n < genome->num_genes)
      genome->usg_count[n] = 0;
    if (gene_no != ULONG_MAX)
      genome->usg_count[gene_no] = 0;
  }
  if (genome->flags & GNM_BP)
  {
    for (i = genome->num_genes - 1; n < i; i--)
      genome->bp_count[i] = genome->bp_count[i - 1];
    if (n < genome->num_genes)
      genome->bp_count[n] = 0;
    if (gene_no != ULONG_MAX)
      genome->bp_count[gene_no] = 0;
  }
#ifdef MEMDEBUG
  verify_Watchdogs();
#endif
}


/*
 * mutate genome at location passed in first argument with the mutationrates
 * specified in the following arguments.
 * Return: -1: error
 *          0: successful
 */

int mutate_genome(GENOME *genome, double m_replacement, double m_insertion, double m_deletion, double m_duplication, double m_factor)
{
  long j, p, gene_no, num_genes;
  unsigned long gene_startpos, gpos, target_pos, gene_length;
  unsigned char flipmask;
  double mf = 1.0;

  if (m_factor && genome->mut_flag)
  {
    errno = 0;
    mf = pow(m_factor, genome->mut_flag);
    if (errno)
    {
      fprintf(stderr, "mutate_genome: error in pow(%f, %ld): ", m_factor, genome->mut_flag);
      perror("");
    }
    else
    {
      m_replacement *= mf;
      m_insertion *= mf;
      m_deletion *= mf;
      m_duplication *= mf;
      if (m_replacement > 1.0)
        m_replacement = 1.0;
      if (m_insertion > 1.0)
        m_insertion = 1.0;
      if (m_deletion > 1.0)
        m_deletion = 1.0;
      if (m_duplication > 1.0)
        m_duplication = 1.0;
    }
  }
#ifdef MEMDEBUG
  verify_Watchdogs();
#endif
  genome->num_mutations++;
  if (m_replacement > 0.0)
  {
    p = next_mutpos(m_replacement, 0);
    while (p < genome->length)
    {
      gene_no = l4_gene_number(genome, p);
      flipmask = 1 << lnd_random(8);
      genome->g[p] ^= flipmask;
      if (l4_promotermask(flipmask))
      {
        if (l4_promoter(genome->g[p]))
        {
          if (resize_genome(genome, genome->length, genome->num_genes + 1) < 0)
            do_error("Error resizing genome while creating new promoter");
          else
            insert_counter(genome, gene_no);
#ifdef MEMDEBUG
          verify_Watchdogs();
#endif
        }
        else
        {
          delete_counter(genome, gene_no);
#ifdef MEMDEBUG
          verify_Watchdogs();
#endif
          if (resize_genome(genome, genome->length, genome->num_genes - 1) < 0)
            do_error("Error resizing genome while destroying promoter");
        }
      }
      else if (gene_no != ULONG_MAX)
      {
        if (genome->flags & GNM_USG)
          genome->usg_count[gene_no] = 0;
        if (genome->flags & GNM_BP)
          genome->bp_count[gene_no] = 0;
      }
      p = next_mutpos(m_replacement, p) + 1;
    }
  }
#ifdef MUTATE_CONSITENCY
  if (l4_num_genes(genome) != genome->num_genes)
    fprintf(stderr, "Inconsistent genome after replacement\n");
#endif
#ifdef MEMDEBUG
  verify_Watchdogs();
#endif
  if (m_insertion > 0.0)
  {
    p = next_mutpos(m_insertion, 0);
    while (p <= genome->length)
    {
      gene_no = l4_gene_number(genome, p);
#ifdef MEMDEBUG
      verify_Watchdogs();
#endif
      if (resize_genome(genome, genome->length + 1, genome->num_genes) < 0)
      {
        do_error("resizing failure during insertion");
        break;
      }
      for (j = genome->length - 1; p < j; j--)
        genome->g[j] = genome->g[j - 1];
      genome->g[p] = (unsigned char) lnd_random(256);
      if (l4_promoter(genome->g[p]))
      {
        if (resize_genome(genome, genome->length, genome->num_genes + 1) < 0)
        {
          do_error("resizing failure during insertion of counter");
          break;
        }
        insert_counter(genome, gene_no);
      }
#ifdef MEMDEBUG
      verify_Watchdogs();
#endif
      p = next_mutpos(m_insertion, p) + 2;
    }
  }
#ifdef MUTATE_CONSITENCY
  if (l4_num_genes(genome) != genome->num_genes)
    fprintf(stderr, "Inconsistent genome after inesertion\n");
#endif
  /* printf("  deleting at random, length = %lu\n", genome->length); */
#ifdef MEMDEBUG
  verify_Watchdogs();
#endif
  if ((genome->length > 0) && (m_deletion > 0.0))
  {
    p = next_mutpos(m_deletion, 0);
    while ((genome->length > 0) && (p < (genome->length - 1)))
    {
      gene_no = l4_gene_number(genome, p);
      if (l4_promoter(genome->g[p]))
      {
        delete_counter(genome, gene_no);
        if (resize_genome(genome, genome->length, genome->num_genes - 1) < 0)
          do_error("resizing failure during deletion of counter");
      }
      for (j = p + 1; j < genome->length; j++)
        genome->g[j - 1] = genome->g[j];
      if (resize_genome(genome, genome->length - 1, genome->num_genes))
      {
        do_error("resizing failure during deletion");
        break;
      }
      p = next_mutpos(m_deletion, p);
    }
  }
#ifdef MUTATE_CONSITENCY
  if (l4_num_genes(genome) != genome->num_genes)
    fprintf(stderr, "Inconsistent genome after deletion\n");
#endif
#ifdef MEMDEBUG
  verify_Watchdogs();
#endif
  if (m_duplication > 0.0)
  {
    num_genes = genome->num_genes;
    p = next_mutpos(m_duplication, 0);
    while (p < num_genes)
    {
      gene_startpos = l4_gene_startposition(genome, p);
      if (gene_startpos == ULONG_MAX)
      {
        do_error("failed to map gene to character position during gene duplication");
        break;
      }
      if (gene_startpos == genome->length - 1)
        gpos = gene_startpos;
      else
      {
        gpos = gene_startpos + 1;
        while ((gpos < genome->length) && !l4_terminator(genome->g[gpos]) && !l4_promoter(genome->g[gpos]))
          gpos++;
        if (gpos == genome->length)
          gpos--;
        if (l4_promoter(genome->g[gpos]))
          gpos--;
      }
      gene_length = gpos - gene_startpos + 1;
      target_pos = genome->length;
      if (resize_genome(genome, genome->length + gene_length, genome->num_genes + 1))
      {
        do_error("resizing failure during gene duplication");
        break;
      }
      for (j = 0; j < gene_length; j++)
        genome->g[target_pos + j] = genome->g[gene_startpos + j];
      insert_counter(genome, genome->num_genes - 1);
      p = next_mutpos(m_duplication, p) + 1;
    }
  }
#ifdef MUTATE_CONSITENCY
  if (l4_num_genes(genome) != genome->num_genes)
    fprintf(stderr, "Inconsistent genome after duplication\n");
#endif
#ifdef MEMDEBUG
  verify_Watchdogs();
#endif
  return (0);
}

