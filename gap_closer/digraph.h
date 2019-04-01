/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-12-20 14:17:33
  *Edit History: 
***********************************************************/

#ifndef XDK_GRAPH_H
#define XDK_GRAPH_H

#include <stdint.h>

#include "mp.h"
#include "bmp.h"

typedef void (*NodeInfoDump) (FILE * fp, void * info);

struct dg_node_s;
typedef struct dg_node_s dg_node_t;

struct dg_edge_s;
typedef struct dg_edge_s dg_edge_t;

struct dg_in_s;
typedef struct dg_in_s dg_in_t;

struct dg_out_s;
typedef struct dg_out_s dg_out_t;

struct digraph_s;
typedef struct digraph_s digraph_t;

struct dg_node_s {
  void * info;
	int32_t id;
  int32_t n_in;
  int32_t n_out;
	int32_t in;
	int32_t out;
};

struct dg_edge_s {
  void * info;
	int32_t id;
	int32_t from_node;
	int32_t to_node;
};

struct dg_in_s {
	int32_t id;
	int32_t edge;
	int32_t from_node;
	int32_t next;
};

struct dg_out_s {
	int32_t id;
	int32_t edge;
	int32_t to_node;
	int32_t next;
};

BMP_DEF (dg_node, dg_node_t);
MP_DEF (dg_edge, dg_edge_t);
MP_DEF (dg_in, dg_in_t);
MP_DEF (dg_out, dg_out_t);

struct digraph_s {
  bmp_t(dg_node) * nodes;
  mp_t(dg_edge) * edges;
  mp_t(dg_in) * ins;
  mp_t(dg_out) * outs;
};

digraph_t * digraph_init (void);
void digraph_clear (digraph_t * dg);
void digraph_free (digraph_t * dg);

dg_node_t * digraph_add1node (digraph_t * dg, void * info);
dg_node_t * digraph_node_at (digraph_t * dg, int32_t node_id);

dg_edge_t * digraph_add1edge (digraph_t * dg, int32_t from_node, int32_t to_node);
dg_edge_t * digraph_edge_between (digraph_t * dg, int32_t from_node, int32_t to_node);

int digraph_dump (digraph_t * dg, const char * file, NodeInfoDump info_dump_func);

#endif
