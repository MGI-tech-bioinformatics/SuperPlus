/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-08-27 18:12:11
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "str.h"
#include "cigar.h"

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

#define INIT_CIGAR_MAX 4

cigar_t *
cigar_init (void)
{
  cigar_t * c;

  c = (cigar_t *) ckalloc (1, sizeof(cigar_t));
  cigar_init2 (c);

  return c;
}

void
cigar_init2 (cigar_t * c)
{
  c->n = 0;
  c->m = INIT_CIGAR_MAX;
  c->c = (uint32_t *) ckalloc (INIT_CIGAR_MAX, sizeof(uint32_t));
}

void
cigar_clear (cigar_t * c)
{
  c->n = 0;
}

void
cigar_free (cigar_t * c)
{
  free (c->c);
  free (c);
}

void
cigar_free2 (cigar_t * c)
{
  free (c->c);
}

void
cigar_add (cigar_t * c, uint32_t cigar_elem)
{
  cigar_resize (c, c->n+1);
  c->c[c->n++] = cigar_elem;
}

void
cigar_resize (cigar_t * c, int32_t cnt)
{
  if (cnt <= c->m)
    return;

  if (c->m <= 0)
    c->m = INIT_CIGAR_MAX;

  while (cnt > c->m)
    c->m <<= 1;
  c->c = (uint32_t *) ckrealloc (c->c, c->m*sizeof(uint32_t));
}

void
cigar_copy (cigar_t * dst, cigar_t * src)
{
  cigar_resize (dst, src->n);
  dst->n = src->n;
  memcpy (dst->c, src->c, src->n<<2);
}

void
cigar_reverse (cigar_t * c)
{
  int32_t i, h;
  uint32_t t;

  if (c->n <= 0)
    warn_mesg ("empty_cigar");

  h = c->n >> 1;
  for (i=0; i<h; ++i) {
    t = c->c[i];
    c->c[i] = c->c[c->n-1-i];
    c->c[c->n-1-i] = t;
  }
}

void
cigar_dump (FILE * fp, cigar_t * c)
{
  int32_t i;
  uint32_t opr;
  uint32_t len;

  for (i=0; i<c->n; ++i) {
    opr = bam_cigar_op (c->c[i]);
    len = bam_cigar_oplen (c->c[i]);
    fprintf (fp, "%d%c", len, BAM_CIGAR_STR[opr]);
  }
}

int
cigar_write (FILE * fp, cigar_t * c)
{
  fwrite (&c->n, 4, 1, fp);
  fwrite (c->c, 4, c->n, fp);
}

int
cigar_read (FILE * fp, cigar_t * c)
{
  fread (&c->n, 4, 1, fp);
  cigar_resize (c, c->n);
  fread (c->c, 4, c->n, fp);
}

void
cigar2str (cigar_t * c, str_t * s)
{
  char * ch;
  int32_t i;
  int32_t old_len;
  uint32_t opr;
  uint32_t len;

  if (c->n <= 0) {
    str_resize (s, 1);

    s->s[0] = '*';
    s->s[1] = '\0';
    s->l = 1;

    return;
  }

  s->l = 0;
  str_resize (s, c->n<<2);
  for (i=0; i<c->n; ++i) {
    opr = bam_cigar_op (c->c[i]);
    len = bam_cigar_oplen (c->c[i]);

    old_len = s->l;
    s->l += int2deci_nbits(len) + 1;
    str_resize (s, s->l);
    ch = s->s + old_len;
    sprintf (ch, "%d%c", len, BAM_CIGAR_STR[opr]);
  }
}

int32_t
cigar2ref_len (cigar_t * c)
{
  int32_t i;
  int32_t len;

  for (i=0,len=0; i<c->n; ++i)
    if (bam_cigar_type(bam_cigar_op(c->c[i])) & 2)
      len += bam_cigar_oplen (c->c[i]);

  return len;
}

int32_t
cigar2qry_len (cigar_t * c)
{
  int32_t i;
  int32_t len;

  for (i=0,len=0; i<c->n; ++i)
    if (bam_cigar_type(bam_cigar_op(c->c[i])) & 1)
      len += bam_cigar_oplen (c->c[i]);

  return len;
}

int
cigar_cleanup (cigar_t * in, cigar_t * out)
{
  int32_t i;
  uint32_t opr;
  uint32_t len;

  cigar_clear (out);
  for (i=0; i<in->n; ++i) {
    opr = bam_cigar_op (in->c[i]);
    len = bam_cigar_oplen (in->c[i]);
    if (len!=0 && (out->n!=0 || opr!=BAM_CDEL)) {
      cigar_add (out, in->c[i]);
    }
  }

  return 0;
}

int
cigar_unclip (cigar_t * in, cigar_t * out)
{
  int32_t i;
  uint32_t opr;

  cigar_clear (out);
  for (i=0; i<in->n; ++i) {
    opr = bam_cigar_op (in->c[i]);
    if (opr != BAM_CSOFT_CLIP
        && opr != BAM_CHARD_CLIP
        && opr != BAM_CPAD)
      cigar_add (out, in->c[i]);
  }

  return 0;
}

int
cigar_has_zero_size_element (cigar_t * c)
{
  int32_t i;
  uint32_t len;

  for (i=0; i<c->n; ++i) {
    len = bam_cigar_oplen (c->c[i]);
    if (len == 0)
      return 1;
  }

  return 0;
}
