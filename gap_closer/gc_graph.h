/*************************************************
 * File Name: gc_graph.h
 * Description: 
 * Author: Chen Xi
 * Mail: chenxi1@genomics.cn
 * Created Time: 2018/12/20
 * Edit History: 
 *************************************************/

#ifndef ONT_GC_GRAPH_H
#define ONT_GC_GRAPH_H

#include <stdint.h>

#include "mp.h"
#include "rseq.h"
#include "ont.h"
#include "sw.h"
#include "str.h"
#include "hash.h"
#include "array.h"
#include "digraph.h"
#include "contig.h"
#include "gap_closer.h"

#define GCV_ANC 0x1 // seq for anchor
#define GCV_SEQ 0x2 // normal seq
#define GCV_GAP 0x4 // gap seq
#define GCV_SCF 0x8 // scaffold sequence
#define GCV_ONT 0x10 // ont sequence
#define GCV_LFR 0x20 // lfr sequence

struct gc_vinfo_s; // gap closer vertex info
typedef struct gc_vinfo_s gc_vinfo_t;

struct gc_graph_s;
typedef struct gc_graph_s gc_graph_t;

struct gc_vinfo_s {
	str_t * s;
  int32_t ctg_id;
  int32_t ctg_pos;
	uint64_t flag;
};
BMP_DEF (gcv, gc_vinfo_t);

struct gc_graph_s {
	digraph_t * g;
	bmp_t(gcv) * vinfos;
  xh_t * k2v_map; // kmer to kmer-vertex map

  sw_t * aligner;
	str_t * qry;
	str_t * tgt;
};

gc_graph_t * gc_graph_init (void);
void gc_graph_clear (gc_graph_t * gcg);
void gc_graph_free (gc_graph_t * gcg);

void gc_vinfo_dump (FILE * fp, void * info);

#endif
