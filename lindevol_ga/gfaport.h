#ifndef H_GFAPORT
#define H_GFAPORT

extern long bas_int(double x);
extern short input_boolean(const char *txt);
extern signed long input_long(const char *txt);
extern unsigned long input_ul(const char *txt);
extern double input_dbl(const char *txt);
extern char *input_str(const char *txt, char *buf, unsigned long max_len);

extern double rnd(void);

#endif



