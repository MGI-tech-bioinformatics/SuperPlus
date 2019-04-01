/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-12-10 16:59:20
  *Edit History: 
***********************************************************/

#ifndef ONT_H
#define ONT_H

#include "mp.h"
#include "rseq.h"
#include "kmer.h"
#include "def.h"

static inline void
okseq_init2 (okseq_t * seq, void * data)
{
  seq->segs = mp_init (oseg, NULL, NULL);
  seq->okmers = mp_init (okmer, NULL, NULL);
  //seq->ctg_ids = arr_init (itg);
}

static inline void
okseq_clear (okseq_t * seq)
{
  mp_clear (oseg, seq->segs, NULL);
	mp_clear (okmer, seq->okmers, NULL);
}

static inline void
okseq_free2 (okseq_t * seq)
{
  mp_free (oseg, seq->segs, NULL);
  mp_free (okmer, seq->okmers, NULL);
}

static inline void
okseq_copy (okseq_t * dst, okseq_t * src)
{
  dst->ont_id = src->ont_id;
  dst->flag = src->flag;
  dst->seq = src->seq;
  mp_copy (okmer, dst->okmers, src->okmers, NULL);
  //arr_copy (itg, dst->ctg_ids, src->ctg_ids);
}

int search_kmers_on_ont_reads (mp_t(rs) * raw_seqs,
    mp_t(ctg) * ctg_seqs, xh_t ** ctg_khashs,
    mp_t(okseq) * okseqs, xh_set_t(kmer) ** anchored_ksets,
    const char * prefix, int n_thread, int kmer_len);

void okseq_init2 (okseq_t * okseq, void * data);

int ont_kseqs_init (mp_t(rs) * raw_seqs, mp_t(okseq) * okseqs);

int ont_kseqs_dump (mp_t(okseq) * okseqs);

#endif
