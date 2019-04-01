/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-12-18 15:35:11
  *Edit History: 
***********************************************************/

#ifndef ONT_GAP_CLOSER_H
#define ONT_GAP_CLOSER_H

#include "lfr.h"

int gap_closing (mp_t(ctg) * ctg_seqs, mp_t(rs) * ont_seqs, lfr_set_t * lfr_seqs,
    mp_t(okseq) * okseqs, mp_t(sf) * scafs, ctg_graph_t * ctg_graph,
    xh_set_t(kmer) ** anchored_ksets, int n_thread);

#endif
