/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-01-02 20:46:37
  *Edit History: 
***********************************************************/

#ifndef CTG_ORDER_H
#define CTG_ORDER_H

#include <stdint.h>

#include "bmp.h"
#include "def.h"
#include "digraph.h"

#define CTG_NODE_DEL 1
#define CTG_GAP      2
#define CTG_N_GAP    4
#define CTG_NEGA_GAP 8

#define CTG_EDGE_SCAF_LINK 1

struct ctg_vinfo_s;
typedef struct ctg_vinfo_s ctg_vinfo_t;

struct ont2gap_s;
typedef struct ont2gap_s ont2gap_t;

struct ctg_einfo_s;
typedef struct ctg_einfo_s ctg_einfo_t;

struct ctg_graph_s;
typedef struct ctg_graph_s ctg_graph_t;

struct ont_node_s;
typedef struct ont_node_s ont_node_t;

struct ctg_vinfo_s {
  int32_t tid;
  int32_t nid;
  int32_t n_okmers;
	//int16_t direct;
	uint32_t flag;

  int64_t l_gap;
  str_t * gap_seq;
};
BMP_DEF (ctg_vinfo, ctg_vinfo_t);

struct ont2gap_s {
  int32_t ont_id;
  int32_t left_okid;
  int32_t right_okid;
};
MP_DEF (o2g, ont2gap_t);

struct ctg_einfo_s {
  uint32_t flag;
  mp_t(o2g) * ont_infos;
};
BMP_DEF (ctg_einfo, ctg_einfo_t);

struct ont_node_s {
	int32_t n_okmers;
	int32_t tid;
	int32_t nid;
	int16_t direct;
	uint16_t flag;

	int64_t oid;
	int64_t beg_okid;;
	int64_t end_okid;
};

struct ctg_graph_s {
  digraph_t * g;
  bmp_t(ctg_vinfo) * vinfos;
  bmp_t(ctg_einfo) * einfos;
};

ctg_graph_t * ctg_graph_init (void);
void ctg_graph_clear (ctg_graph_t * g);
void ctg_graph_free (ctg_graph_t * g);

int map_ont2contigs (mp_t(okseq) * okseqs, ctg_graph_t * g, mp_t(ctg) * ctg_seqs);

int fix_ont1 (mp_t(okseq) * okseqs, ctg_graph_t * g,
    mp_t(ctg) * ctg_seqs, mp_t(sf) * scafs);

#endif
