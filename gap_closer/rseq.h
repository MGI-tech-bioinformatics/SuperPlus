/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-12-10 17:31:16
  *Edit History: 
***********************************************************/

#ifndef XDK_FQ_H
#define XDK_FQ_H

#include <stdint.h>

#include "mp.h"
#include "str.h"

/*****************************************************************************/
/*-------------------------- read seq definition ----------------------------*/
/*****************************************************************************/

struct rseq_s;
typedef struct rseq_s rseq_t;

struct rseq_s {
  char * b, * q;
  int32_t l, m;
};

rseq_t * rseq_init (void);
void rseq_init2 (rseq_t * r, void * data);
void rseq_resize (rseq_t * r, int32_t new_size);
void rseq_clear (rseq_t * r);
void rseq_free (rseq_t * r);
void rseq_free2 (rseq_t * r);
void rseq_assign (rseq_t * r, const char * name, const char * anno,
    const char * bases, const char * quals);
void rseq_reverse_complement (rseq_t * r);

MP_DEF (rs, rseq_t);

/*****************************************************************************/
/*---------------------- pair end read seq definition -----------------------*/
/*****************************************************************************/

struct pe_seq_s;
typedef struct pe_seq_s pe_seq_t;

struct pe_seq_s {
  str_t * n;
	char * b[2];
	char * q[2];
	int32_t l[2];
	int32_t m[2];
};

pe_seq_t * pe_seq_init (void);
void pe_seq_init2 (pe_seq_t * p, void * data);
void pe_seq_resize (pe_seq_t * p, int32_t new_size1, int32_t new_size2);
void pe_seq_clear (pe_seq_t * p);
void pe_seq_free (pe_seq_t * p);
void pe_seq_free2 (pe_seq_t * p);
void pe_seq_assign (pe_seq_t * p, int end, const char * name, const char * anno,
    const char * bases, const char * quals);

MP_DEF (pe, pe_seq_t);

/*****************************************************************************/
/*-------------------------- block of read seqs -----------------------------*/
/*****************************************************************************/

struct rseq_blk_s;
typedef struct rseq_blk_s rseq_blk_t;

struct rseq_blk_s {
  rseq_t * seqs;
  int32_t n, m;
};

void rseq_blk_init2 (rseq_blk_t * blk);

void rseq_blk_clear (rseq_blk_t * blk);

void rseq_blk_free2 (rseq_blk_t * blk);

/*****************************************************************************/
/*------------------------- single end fastq reader -------------------------*/
/*****************************************************************************/

/*
struct sefq_reader_s;
typedef struct sefq_reader_s sefq_reader_t;

struct sefq_reader_s {
  int nth;
  char * file;
  pthread_t pid;
  rseq_blk_t ** blks;
};

sefq_reader_t * sefq_reader_init (const char * fq_file, int nth);

void sefq_reader_destroy (sefq_reader_t * reader);
*/

mp_t(rs) * sefq_load (const char * fq_file);

/*****************************************************************************/
/*------------------------- paired end fastq reader -------------------------*/
/*****************************************************************************/

/*
struct pefq_reader_s;
typedef struct pefq_reader_s pefq_reader_t;

struct pefq_reader_s {

};

pefq_reader_t * pefq_reader_init (const char * fq1_file, const char * fq2_file, int nth);
*/

mp_t(pe) * pefq_load (const char * fq1_file, const char * fq2_file);

#endif
