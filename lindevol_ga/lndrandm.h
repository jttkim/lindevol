#ifndef H_LNDRANDM
#define H_LNDRANDM

extern unsigned long lnd_random(unsigned long range);
extern double lnd_rnd(void);
extern int write_rndgenerator_state(FILE *f);
extern int read_rndgenerator_state(FILE *f);
extern void seed_lnd_random(int seed);
extern void random_shuffle(long num, long *a);

#endif

