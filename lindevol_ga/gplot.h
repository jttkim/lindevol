/* header file for lnd2 gplot functions */

#ifndef H_GPLOT
#define H_GPLOT

extern long max_numgenes(void);
extern long num_usedgenes(const LND_GENOME *lnd_genome);
extern long max_num_usedgenes(void);
extern long gene2statespec(const GENOME *genome, long pos, STATE_SPEC *statespec);
extern int calculate_statespecs(const unsigned char *gene, STATE_SPEC *statespec);
extern int statespec_match(const STATE_SPEC *current_statespec, int n, const STATE_SPEC *statespec);
extern void list_genome(FILE *f, long plant_no, int used_only);
extern int postscript_genebox(FILE *f, const GENOME *genome, long gene_no, double x0, double y0, double width, double height);
extern int postscript_connectgraph(FILE *f, const LND_GENOME *lnd_genome, int used_only,
        double x0, double y0, double width, double height, double gbox_width, double gbox_height, unsigned long counter_flags);

#endif /* H_GPLOT */

