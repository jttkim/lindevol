#ifndef URANDOM_H
#define URANDOM_H

extern void ulong_srandom( unsigned int x);
extern unsigned long *ulong_initstate(unsigned int seed, unsigned long *arg_state);
extern unsigned long *ulong_setstate(unsigned long *arg_state);
extern long int ulong_random(void);
extern unsigned long urandom_long(unsigned long range);
extern double urandom_double(void);
extern double urandom_gauss(void);
extern int write_urandom_state(FILE *f);
extern int read_urandom_state(FILE *f);
extern void *urandom_shuffle(size_t num, size_t s, void *a);
extern int main(int argc, char *argv[]);

#endif

