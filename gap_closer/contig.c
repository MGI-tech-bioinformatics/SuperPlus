/*************************************************
 * File Name: contig.c
 * Description: 
 * Author: Chen Xi
 * Mail: chenxi1@genomics.cn
 * Created Time: 2018/12/09
 * Edit History: 
 *************************************************/

#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "mp.h"
#include "def.h"
#include "utils.h"
#include "contig.h"
#include "kmer.h"
#include "ctg_graph.h"

static void
contig_init2 (ctg_t * seq, void * data)
{
	seq->seq = str_init ();
  seq->ont_ids = arr_init (itg);
}

static void
contig_free2 (ctg_t * seq)
{
	int32_t i;

  str_free (seq->seq);
  arr_free (itg, seq->ont_ids);

	if (seq->kmers) {
		for (i=0; i<seq->n_kmer; ++i)
			kmer_free2 (seq->kmers+i);
		free (seq->kmers);
	}
}

static dg_node_t *
contig_node_create (ctg_graph_t * ctg_graph, ctg_t * seq)
{
	ctg_vinfo_t * vinfo;
  dg_node_t * node;

	vinfo = bmp_alloc (ctg_vinfo, ctg_graph->vinfos);
	vinfo->tid = seq->id;

	node = digraph_add1node (ctg_graph->g, vinfo);
  seq->node_id = node->id;

  return node;
}

static int
contig_edge_create (ctg_graph_t * ctg_graph, int32_t from_node_id, int32_t to_node_id)
{
	ctg_einfo_t * einfo;
	dg_edge_t * edge;

	einfo = bmp_alloc (ctg_einfo, ctg_graph->einfos);
	einfo->flag = CTG_EDGE_SCAF_LINK;

	edge = digraph_add1edge (ctg_graph->g, from_node_id, to_node_id);
	edge->info = einfo;

	return 0;
}

static void
add1scaf (mp_t(ctg) * ctgs, mp_t(sf) * scafs, int32_t beg_ctg_id, int32_t end_ctg_id)
{
  int32_t i, j;
  ctg_t * ctg;
  scaf_t * sf;

  sf = mp_alloc (sf, scafs);
  sf->id = mp_cnt(scafs) - 1;

  sf->ctg_ids->n = end_ctg_id - beg_ctg_id;
  arr_resize (itg, sf->ctg_ids, sf->ctg_ids->n);
  for (i=0,j=beg_ctg_id; j<end_ctg_id; ++i,++j) {
    sf->ctg_ids->arr[i] = j;

    ctg = mp_at (ctg, ctgs, j);
    ctg->link_id = sf->id;
  }
}

mp_t(ctg) *
contig_seqs_load (const char * seq_file, ctg_graph_t * ctg_graph, mp_t(sf) * scafs)
{
	char * line;
  int pre_st;
  int cur_st;
  time_t time_beg;
	int32_t i, j;
  int32_t beg;
  int32_t n_sf_seqs;
  int32_t pre_node_id;
  int32_t beg_ctg_id;
  int32_t end_ctg_id;
  int32_t gap_size;
	FILE * fp;
	ctg_t * seq;
  str_t * sf_seq;
	mp_t(ctg) * set;
  str_set_t * sf_seqs;
  dg_node_t * node;

  time (&time_beg);

	set = mp_init (ctg, contig_init2, NULL);
  sf_seqs = str_set_init ();

	line = ALLOC_LINE;
	fp = ckopen (seq_file, "r");
	while (fgets(line, LINE_MAX, fp)) {
		if (*line == '>') {
      sf_seq = str_set_alloc (sf_seqs);
			continue;
		}
		chomp (line);
		str_append (sf_seq, line, strlen(line));
	}
	fclose (fp);
	free (line);

  n_sf_seqs = str_set_cnt (sf_seqs);
  for (i=0; i<n_sf_seqs; ++i) {
    sf_seq = str_set_at (sf_seqs, i);

    pre_st = NOT_N;
    gap_size = -1;
    pre_node_id = -1;
    beg_ctg_id = mp_cnt (set);
    for (j=0,beg=0; j<sf_seq->l; ++j) {
      cur_st = is_N (sf_seq->s[j]);
      switch (pre_st) {
      case NOT_N:
        switch (cur_st) {
        case IS_N:
			    seq = mp_alloc (ctg, set);
          seq->id = mp_cnt(set) - 1;
          seq->seq->l = j - beg;
          str_resize (seq->seq, seq->seq->l);
          memcpy (seq->seq->s, sf_seq->s+beg, seq->seq->l);
          seq->seq->s[seq->seq->l] = '\0';
          seq->l_pre_gap = gap_size;
          seq->link_id = -1;
          gap_size = 1;

          // add contig to digraph
          node = contig_node_create (ctg_graph, seq);
          if (pre_node_id != -1)
            contig_edge_create (ctg_graph, pre_node_id, node->id);
          pre_node_id = node->id;

          break;
        }
        break;

      case IS_N:
        switch (cur_st) {
        case NOT_N:
          beg = j;
          break;
        case IS_N:
          ++gap_size;
          break;
        }
        break;
      }
      pre_st = cur_st;
    }

    if (pre_st == NOT_N) {
		  seq = mp_alloc (ctg, set);
      seq->id = mp_cnt(set) - 1;
      seq->seq->l = j - beg;
      str_resize (seq->seq, seq->seq->l);
      memcpy (seq->seq->s, sf_seq->s+beg, seq->seq->l);
      seq->seq->s[seq->seq->l] = '\0';
      seq->l_pre_gap = gap_size;
      seq->link_id = -1;

      // add contig to digraph
      node = contig_node_create (ctg_graph, seq);
      if (pre_node_id != -1)
        contig_edge_create (ctg_graph, pre_node_id, node->id);
      pre_node_id = node->id;
    }
    end_ctg_id = mp_cnt (set);

    add1scaf (set, scafs, beg_ctg_id, end_ctg_id);
  }

  for (i=0; i<mp_cnt(set); ++i) {
    seq = mp_at (ctg, set, i);
    seq->n_kmer = 0;
    seq->m_kmer = seq->seq->l;
    seq->kmers = (kmer_t *) ckalloc (seq->m_kmer, sizeof(kmer_t));
  }

  str_set_free (sf_seqs);

	printf ("Total Scaffolds: %d, Total Contigs: %ld, cost: %lds\n", n_sf_seqs, mp_cnt(set), time(NULL)-time_beg);

	return set;
}

void
contig_seqs_info_clear (mp_t(ctg) * seqs)
{
	int32_t i;
	int32_t n_ctgs;
	ctg_t * ctg;

	n_ctgs = mp_cnt (seqs);
	for (i=0; i<n_ctgs; ++i) {
		ctg = mp_at (ctg, seqs, i);
		ctg->n_kmer = 0;
		arr_clear (itg, ctg->ont_ids);
	}
}

void
contig_seqs_free (mp_t(ctg) * seqs)
{
  int64_t i;
  ctg_t * seq;

  for (i=0; i<mp_cnt(seqs); ++i) {
    seq = mp_at (ctg, seqs, i);
    contig_free2 (seq);
  }
  //mp_free (ctg, seqs, contig_free2);
}
