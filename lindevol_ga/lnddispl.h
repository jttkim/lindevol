#ifndef H_LNDDISPL
#define H_LNDDISPL

extern int open_world_file(const char *simname);
extern void close_world_file(void);
extern void display_on(unsigned long generation);
extern void display_off(void);
extern void draw_cell(unsigned long x, unsigned long y);
extern void poll_user_interface(void);
extern void init_world_display(const char *simname);
extern void close_world_display(void);

#endif

