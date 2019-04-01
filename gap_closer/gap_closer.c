/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-12-18 15:48:30
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mp.h"
#include "sw.h"
#include "lfr.h"
#include "ont.h"
#include "hash.h"
#include "kmer.h"
#include "utils.h"
#include "contig.h"
#include "gc_graph.h"
#include "gap_closer.h"
#include "ctg_graph.h"

int
gap_closing (mp_t(ctg) * ctg_seqs, mp_t(rs) * ont_seqs,
    lfr_set_t * lfr_seqs, mp_t(okseq) * okseqs, mp_t(sf) * scafs,
    ctg_graph_t * ctg_graph, xh_set_t(kmer) ** anchored_ksets, int n_thread)
{

  // map ont to contigs
  map_ont2contigs (okseqs, ctg_graph, ctg_seqs);

  // fix errors only using ont
  fix_ont1 (okseqs, ctg_graph, ctg_seqs, scafs);

  return 0;
}
