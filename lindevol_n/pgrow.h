#ifndef H_PGROW
#define H_PGROW

#ifdef __cplusplus
extern "C" {
#endif

extern long create_cell(long p, long x, long y);
extern long divide(long plant_no, long cell_no, int direction);
extern void cell_activity(long plant_no, long cell_no, long action, unsigned char gene_output);
extern void plant_growth(long plant_no);

#ifdef __cplusplus
}
#endif

#endif /* H_PGROW */

