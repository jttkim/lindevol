/* header file for lnd2 gplot functions */

#ifndef H_GPLOT
#define H_GPLOT

#ifdef __cplusplus
extern "C" {
#endif

extern int num_setbits(unsigned long x);
extern long max_numgenes(void);
extern long num_usedgenes(const PLANT *plant);
extern long max_num_usedgenes(void);
extern int calculate_statespecs(const unsigned char *gene, STATE_SPEC *statespec);
extern int statespec_match(const STATE_SPEC *current_statespec, int n, const STATE_SPEC *statespec);
extern void list_genome(FILE *f, long plant_no, int used_only);
extern int postscript_counterbox(FILE *f, const GENOME *genome, long gene_no, double x0, double y0, double width, double height, double fontheight, unsigned long flags);
extern int postscript_genebox(FILE *f, const GENOME *genome, long gene_no, double x0, double y0, double width, double height);
extern int postscript_connectgraph(FILE *f, const PLANT *plant, int used_only,
        double x0, double y0, double width, double height, double gbox_width, double gbox_height, unsigned long counter_flags);

#ifdef __cplusplus
}
#endif

#endif /* H_GPLOT */

