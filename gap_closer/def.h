/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-01-02 21:53:33
  *Edit History: 
***********************************************************/

#ifndef GC_DEF_H
#define GC_DEF_H

#include <stdint.h>

#include "mp.h"
#include "rseq.h"
#include "hash.h"
#include "kseq1.h"
#include "array.h"

struct kmer_s;
typedef struct kmer_s kmer_t;

//struct kpos_s;
//typedef struct kpos_s kpos_t;

struct ctg_s;
typedef struct ctg_s ctg_t;

struct scaf_s;
typedef struct scaf_s scaf_t;

struct ont_seg_s;
typedef struct ont_seg_s ont_seg_t;

struct ont_kmer_s;
typedef struct ont_kmer_s ont_kmer_t;

struct okseq_s;
typedef struct okseq_s okseq_t;

/*****************************************************************************/
/*---------------------------- Kmer Definitions -----------------------------*/
/*****************************************************************************/

//struct kpos_s {
//  int64_t tid;
//  int32_t pos;
//  uint32_t flag;
//};

//MP_DEF (kpos, kpos_t);

#define KMER_REV 0x1
#define KMER_ONT 0x2 // found on ONT reads
#define KMER_NGS 0x4 // found on NGS reads

struct kmer_s {
	kseq1_t kseq;
  int32_t hs_id; // hash set id
	int32_t tid;
	int32_t pos;
  uint16_t flag;
  int16_t kmer_len;
  //mp_t(kpos) * ont_pos;
};

HASH_SET_DEF (kmer, kmer_t);

/*****************************************************************************/
/*--------------------------- Contig Definitions ----------------------------*/
/*****************************************************************************/

ARR_DEF (itg, int32_t);
ARR_DEF (dbl, double);

struct ctg_s {
	str_t * seq;
  int32_t id;
  int32_t n_kmer;
  int32_t m_kmer;
  int32_t node_id;
  int32_t link_id;
  int32_t l_pre_gap;
  struct kmer_s * kmers;
  arr_t(itg) * ont_ids;
};

MP_DEF (ctg, ctg_t);

struct scaf_s {
  arr_t(itg) * ctg_ids;
  int32_t id;
};

static inline void scaf_init2 (scaf_t * sf_link, void * data)
{
  sf_link->ctg_ids = arr_init (itg);
}

static inline void scaf_free2 (scaf_t * sf_link)
{
  arr_free (itg, sf_link->ctg_ids);
}

MP_DEF (sf, scaf_t);

/*****************************************************************************/
/*----------------------------- ONT Definitions -----------------------------*/
/*****************************************************************************/

#define ONT_KMER_REV 0x1
#define ONT_SCAF_REV 0x2
#define ONT_DEL      0x4
#define ONT_MISLEAD  0x8
#define ONT_NO_ANK   0x10
#define ONT_NO_LINK  0x20
#define ONT_NEW_LINK 0x40

#define ONT_KMER_REV_MASK 0xfffe
#define ONT_SCAD_REV_MASK 0xfffd
#define ONT_DEL_MASK      0xfffb

struct ont_kmer_s {
  kmer_t * kmer;
  int32_t ont_pos;
  int16_t hs_id; // hash set id
  uint16_t flag;
};

MP_DEF (okmer, ont_kmer_t);

struct ont_seg_s {
  int32_t beg;
  int32_t end;
};

MP_DEF (oseg, ont_seg_t);

struct okseq_s {
  int ont_id;
  uint64_t flag;
	rseq_t * seq;
  mp_t(okmer) * okmers;
  mp_t(oseg) * segs;
  //arr_t(itg) * ctg_ids;
};

MP_DEF (okseq, okseq_t);

#endif
