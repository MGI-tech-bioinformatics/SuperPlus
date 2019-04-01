#ifndef GC_LFR_H
#define GC_LFR_H

#include <stdint.h>

#include "mp.h"
#include "contig.h"
#include "hash.h"
#include "rseq.h"

struct lfr_s;
typedef struct lfr_s lfr_t;

struct lfr_bc_s;
typedef struct lfr_bc_s lfr_bc_t;

struct lfr_set_s;
typedef struct lfr_set_s lfr_set_t;

struct lfr_s {
  pe_seq_t * r;
  uint64_t bc;
  int64_t bc_id;
  int16_t hs_id;
};

MP_DEF (lfr, lfr_t);

struct lfr_bc_s {
  uint64_t bc;
  int64_t bc_id;
  mp_t(lfr) * reads;
};

MP_DEF (lbc, lfr_bc_t);

struct lfr_set_s {
  mp_t(lbc) * bcs;
  mp_t(pe) * pes;
  mp_t(lfr) * lfrs;
  int64_t * pth_idx_arr;
  int64_t n_thread;
};

int ankor_lfr_reads (mp_t(ctg) * ctgs, xh_t ** ctg_khashs,
    lfr_set_t * lfr_seqs, int n_thread, int kmer_len,
    const char * fq1_file, const char * fq2_file);

void lfr_free (lfr_set_t * set);

#endif
