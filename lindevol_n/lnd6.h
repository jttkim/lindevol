#ifndef H_LND4
#define H_LND4

#include "lndtypes.h"

#define free0(p) free(p); (p) = NULL

#ifdef __cplusplus
extern "C" {
#endif

extern int mutate_genome(GENOME *genome, double m_replacement, double m_insertion, double m_deletion, double m_duplication, double m_factor);

extern unsigned char cell_state(long plant_no, long cell_no);
extern long compute_statespec(const GENOME *genome, long gpos, STATE_SPEC *statespec, long *num_bitchecks);
extern long process_cell(long plant_no, long cell_no, unsigned char *p_gene_output, const GENE_SPEC *gene_spec);
extern long gene_activity(unsigned char gene_output);
extern int write_named_savefile(const char *save_file_name);
extern int write_savefile(void);
extern int load_named_savefile(const char *save_file_name);
extern int load_savefile(void);

extern int l4_promoter(unsigned char c);
extern int l4_promotermask(unsigned char c);
extern int l4_terminator(unsigned char c);
extern unsigned long l4_num_genes(const GENOME *genome);
extern unsigned long l4_gene_number(const GENOME *genome, unsigned long pos);
extern unsigned long l4_gene_startposition(const GENOME *genome, unsigned long g_num);
extern int calculate_activity_codes(GSYS_PARAMETERS *gp);
extern void cell_activity(long plant_no, long cell_no, long action, unsigned char gene_output);
extern void plant_growth(long plant_no);

#ifdef __cplusplus
}
#endif

#include "pgrow.h"

#endif /* H_LND4 */

