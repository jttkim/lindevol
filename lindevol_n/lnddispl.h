#ifndef H_LNDDISPL
#define H_LNDDISPL

#ifdef __cplusplus
extern "C" {
#endif

extern int open_world_file(const char *simname, const char *mode);
extern void close_world_file(void);
extern void display_start(long generation);
extern void display_done(void);
extern void display_cell(long x, long y);
extern void display_plant_killed(long plant_no);
extern void display_plant_mutated(long plant_no);
extern void poll_user_interface(void);
extern void init_world_display(const char *simname);
extern void close_world_display(void);

#ifdef __cplusplus
}
#endif

#endif

