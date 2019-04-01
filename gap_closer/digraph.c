/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-12-20 14:36:47
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "digraph.h"

digraph_t *
digraph_init (void)
{
  digraph_t * g;

	g = (digraph_t *) ckmalloc (sizeof(digraph_t));
	g->nodes = bmp_init (dg_node, 4096, NULL, NULL);
	g->edges = mp_init (dg_edge, NULL, NULL);
	g->ins = mp_init (dg_in, NULL, NULL);
	g->outs = mp_init (dg_out, NULL, NULL);

	return g;
}

void
digraph_clear (digraph_t * dg)
{
	bmp_clear (dg_node, dg->nodes, NULL);
	mp_clear (dg_edge, dg->edges, NULL);
	mp_clear (dg_in, dg->ins, NULL);
	mp_clear (dg_out, dg->outs, NULL);
}

void
digraph_free (digraph_t * dg)
{
	bmp_free (dg_node, dg->nodes, NULL);
	mp_free (dg_edge, dg->edges, NULL);
	mp_free (dg_in, dg->ins, NULL);
	mp_free (dg_out, dg->outs, NULL);
  free (dg);
}

dg_node_t *
digraph_add1node (digraph_t * dg, void * info)
{
  dg_node_t * node;

  node = bmp_alloc (dg_node, dg->nodes);
  node->info = info;
  node->id = bmp_cnt(dg->nodes) - 1;
  node->n_in = 0;
  node->n_out = 0;
  node->in = -1;
  node->out = -1;

  return node;
}

dg_node_t *
digraph_node_at (digraph_t * dg, int32_t node_id)
{
  return bmp_at (dg_node, dg->nodes, node_id);
}

dg_edge_t *
digraph_add1edge (digraph_t * dg, int32_t from_node, int32_t to_node)
{
  dg_in_t * in;
  dg_out_t * out;
  dg_node_t * node;
	dg_edge_t * edge;

  edge = mp_alloc (dg_edge, dg->edges);
  edge->id = mp_cnt(dg->edges) - 1;
  edge->from_node = from_node;
  edge->to_node = to_node;

  node = bmp_at (dg_node, dg->nodes, from_node);
  out = mp_alloc (dg_out, dg->outs);
  out->id = mp_cnt(dg->outs) - 1;
  out->edge = edge->id;
  out->to_node = to_node;
  out->next = node->out;
  node->out = out->id;
  ++node->n_out;

  node = bmp_at (dg_node, dg->nodes, to_node);
  in = mp_alloc (dg_in, dg->ins);
  in->id = mp_cnt(dg->ins) - 1;
  in->edge = edge->id;
  in->from_node = from_node;
  in->next = node->in;
  node->in = in->id;
  ++node->n_in;

  return edge;
}

dg_edge_t *
digraph_edge_between (digraph_t * dg, int32_t from_node, int32_t to_node)
{
  int32_t out;
  dg_out_t * link;
  dg_node_t * node;

  node = bmp_at (dg_node, dg->nodes, from_node);
  out = node->out;
  while (out != -1) {
    link = mp_at (dg_out, dg->outs, out);
    if (link->to_node == to_node)
      return mp_at (dg_edge, dg->edges, link->edge);
    out = link->next;
  }

  return NULL;
}

int
digraph_dump (digraph_t * dg, const char * file, NodeInfoDump info_dump_func)
{
  int idx;
  FILE * fp;
  dg_node_t * node;
  dg_edge_t * edge;
  dg_out_t * out;

  fp = ckopen (file, "w");

  fprintf (fp, "digraph G {\n");
  fprintf (fp, "\tsize = \"512,512\";\n");

  bmp_iter_init (dg_node, dg->nodes);
  while (NULL != (node=bmp_iter_next(dg_node,dg->nodes))) {
    fprintf (fp, "\tV%d [lable=\"", node->id);
    info_dump_func (fp, node->info);
    fprintf (fp, "\"];\n");

    idx = node->out;
    while (idx != -1) {
      out = mp_at (dg_out, dg->outs, idx);
      edge = mp_at (dg_edge, dg->edges, out->edge);

      if (edge->from_node != node->id)
        err_mesg ("edge->from_node != node->id");
      if (edge->to_node != out->to_node)
        err_mesg ("edge->to_node != out->to_node");

      fprintf (fp, "\tV%d -> V%d;\n", edge->from_node, edge->to_node);

      idx = out->next;
    }
  }

  fprintf (fp, "}\n");
}
