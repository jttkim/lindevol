#ifndef H_LNDIO
#define H_LNDIO

extern void prepare_filenames(void);
extern int calculate_activity_codes(GSYS_PARAMETERS *gp);
extern long gene_activity(unsigned char gene_output);
extern int open_pro_file(const char *mode);
extern int open_dmt_file(const char *mode);
extern int open_stb_file(const char *mode);
extern int open_genome_file(const char *mode);
extern int open_dst_file(const char *mode);
extern int open_bpe_file(const char *mode);
extern void close_data_files(void);
extern void write_pro(void);
extern void write_stb(void);
extern void write_genomes(void);
extern void write_dmt(void);
extern void write_dst(void);
extern void write_bpe(void);
extern int write_named_savefile(const char *savefile_name);
extern int load_named_savefile(const char *savefile_name);
extern int write_savefile(void);
extern int load_savefile(void);
extern int open_savetime_file(const char *fname);
extern void close_savetime_file(void);
extern long savetime_next(long generation);

#endif



