/* include file for all modules in lndlib */

#ifndef H_LNDLIBIN
#define H_LNDLIBIN

#include <stdio.h>

#include "lndvals.h"
#include "lndtypes.h"

#include "pgrow.h"

#ifdef MEMDEBUG
#  ifdef __atarist__
#    include <memdebug.h>
#  else
#    include "memdebug.h"
#  endif
#endif


#ifdef __cplusplus
extern "C" {
#endif

extern FILE *pro_file;     /* the protocol file */
extern FILE *dmt_file;     /* the distance matrices file */
extern FILE *dst_file;     /* the distance distribution file */
extern FILE *jf_file;      /* the jf tree file */
extern FILE *genome_file;  /* the genome file */
extern FILE *bpe_file;     /* the Bedau and Packard evolutionary activity waves file */
extern FILE *pixel_file;

extern char pro_file_name[MAX_SLEN];
extern char dmt_file_name[MAX_SLEN];
extern char dst_file_name[MAX_SLEN];
extern char jf_file_name[MAX_SLEN];
extern char genome_file_name[MAX_SLEN];
extern char bpe_file_name[MAX_SLEN];
extern char save_file_name[MAX_SLEN];
extern char savetime_file_name[MAX_SLEN];
extern char pixel_file_name[MAX_SLEN];

#ifdef __cplusplus
}
#endif

#include "lndlib.h"

#endif

