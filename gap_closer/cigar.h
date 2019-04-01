/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-08-27 18:02:32
  *Edit History: 
***********************************************************/

#ifndef XDK_CIGAR_H
#define XDK_CIGAR_H

#include <stdio.h>
#include <stdint.h>

#include "str.h"

struct cigar_s;
typedef struct cigar_s cigar_t;

struct cigar_s {
  uint32_t * c;
  int32_t n, m;
};

cigar_t * cigar_init (void);
void cigar_init2 (cigar_t * c);

void cigar_clear (cigar_t * c);

void cigar_free (cigar_t * c);
void cigar_free2 (cigar_t * c);

void cigar_add (cigar_t * c, uint32_t cigar_elem);

void cigar_resize (cigar_t * c, int32_t cnt);

void cigar_copy (cigar_t * dst, cigar_t * src);

void cigar_reverse (cigar_t * c);

void cigar_dump (FILE * fp, cigar_t * c);
int cigar_write (FILE * fp, cigar_t * c);
int cigar_read (FILE * fp, cigar_t * c);

void cigar2str (cigar_t * c, str_t * s);

int32_t cigar2ref_len (cigar_t * c);
int32_t cigar2qry_len (cigar_t * c);

int cigar_cleanup (cigar_t * in, cigar_t * out);

int cigar_unclip (cigar_t * in, cigar_t * out);

int cigar_has_zero_size_element (cigar_t * c);

/*****************************************************************************/
/*-------------------------- Cigar Related Macros ---------------------------*/
/*****************************************************************************/

/*
#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6
#define BAM_CEQUAL      7
#define BAM_CDIFF       8
#define BAM_CBACK       9

#define BAM_CIGAR_STR   "MIDNSHP=XB"
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  0xf
#define BAM_CIGAR_TYPE  0x3C1A7

#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
#define bam_cigar_opchr(c) (BAM_CIGAR_STR[bam_cigar_op(c)])
#define bam_cigar_gen(l, o) ((l)<<BAM_CIGAR_SHIFT|(o))
 
#define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3)
*/

#endif
