/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-01-02 20:47:29
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "def.h"
#include "ctg_graph.h"

static void
ctg_vinfo_init2 (ctg_vinfo_t * vinfo, void * data)
{
  vinfo->gap_seq = str_init ();
}

static void
ctg_vinfo_clear (ctg_vinfo_t * vinfo)
{
  str_clear (vinfo->gap_seq);
}

static void
ctg_vinfo_free2 (ctg_vinfo_t * vinfo)
{
  str_free (vinfo->gap_seq);
}

static void
ctg_einfo_init2 (ctg_einfo_t * einfo, void * data)
{
  einfo->ont_infos = mp_init (o2g, NULL, NULL);
}

static void
ctg_einfo_clear (ctg_einfo_t * einfo)
{
  mp_clear (o2g, einfo->ont_infos, NULL);
}

static void
ctg_einfo_free2 (ctg_einfo_t * einfo)
{
  mp_free (o2g, einfo->ont_infos, NULL);
}

ctg_graph_t *
ctg_graph_init (void)
{
  ctg_graph_t * g;

  g = (ctg_graph_t *) ckalloc (1, sizeof(ctg_graph_t));
  g->g = digraph_init ();
  g->vinfos = bmp_init (ctg_vinfo, 4096, ctg_vinfo_init2, NULL);
  g->einfos = bmp_init (ctg_einfo, 4096, ctg_einfo_init2, NULL);

  return g;
}

void
ctg_graph_clear (ctg_graph_t * g)
{
  digraph_clear (g->g);
  bmp_clear (ctg_vinfo, g->vinfos, ctg_vinfo_clear);
  bmp_clear (ctg_einfo, g->einfos, ctg_einfo_clear);
}

void
ctg_graph_free (ctg_graph_t * g)
{
  digraph_free (g->g);
  bmp_free (ctg_vinfo, g->vinfos, ctg_vinfo_free2);
  bmp_free (ctg_einfo, g->einfos, ctg_einfo_free2);
}

/*****************************************************************************/
/*--------------------------- Map ONT to Contigs ----------------------------*/
/*****************************************************************************/

#define FORW 0
#define BACK 1

#define MIN_ANKOR_KMERS 100

MP_DEF (ont_node, ont_node_t);

static int
ont_node_init (ont_node_t * ont_node, int32_t node_id, int32_t * direct, int32_t ctg_id,
		int32_t ont_id, int32_t beg_okid, int32_t end_okid, okseq_t * okseq)
{
  int32_t i;
  ont_kmer_t * okmer;

  if (direct[FORW] > 20*direct[BACK]) {
//    if (direct[FORW] < MIN_ANKOR_KMERS) {
//      ont_node->flag = CTG_NODE_DEL;
//      return 0;
//    }

    ont_node->direct = FORW;
    ont_node->n_okmers = direct[FORW];
    ont_node->tid = ctg_id;
    ont_node->nid = node_id;
    ont_node->flag = 0;
		ont_node->oid = ont_id;

    if (direct[BACK] > 0) {
      for (i=beg_okid; i<=end_okid; ++i) {
        okmer = mp_at (okmer, okseq->okmers, i);
        if (NULL == okmer->kmer)
          continue;
        if (!(okmer->flag & ONT_SCAF_REV))
          break;
      }
      assert (i <= end_okid);
		  ont_node->beg_okid = okmer->ont_pos;

      for (i=end_okid; i>=beg_okid; --i) {
        okmer = mp_at (okmer, okseq->okmers, i);
        if (NULL == okmer->kmer)
          continue;
        if (!(okmer->flag & ONT_SCAF_REV))
          break;
      }
      assert (i >= beg_okid);
		  ont_node->end_okid = okmer->ont_pos;
    } else {
      ont_node->beg_okid = beg_okid;
      ont_node->end_okid = end_okid;
    }
  } else if (direct[BACK] > 20*direct[FORW]) {
//    if (direct[BACK] < MIN_ANKOR_KMERS) {
//      ont_node->flag = CTG_NODE_DEL;
//      return 0;
//    }

    ont_node->direct = BACK;
    ont_node->n_okmers = direct[BACK];
    ont_node->tid = ctg_id;
    ont_node->nid = node_id;
    ont_node->flag = 0;
		ont_node->oid = ont_id;

    if (direct[FORW] > 0) {
      for (i=beg_okid; i<=end_okid; ++i) {
        okmer = mp_at (okmer, okseq->okmers, i);
        if (NULL == okmer->kmer)
          continue;
        if (okmer->flag & ONT_SCAF_REV)
          break;
      }
      assert (i <= end_okid);
		  ont_node->beg_okid = okmer->ont_pos;

      for (i=end_okid; i>=beg_okid; --i) {
        okmer = mp_at (okmer, okseq->okmers, i);
        if (NULL == okmer->kmer)
          continue;
        if (okmer->flag & ONT_SCAF_REV)
          break;
      }
      assert (i >= beg_okid);
		  ont_node->end_okid = okmer->ont_pos;
    } else {
      ont_node->beg_okid = beg_okid;
      ont_node->end_okid = end_okid;
    }
  } else {
    ont_node->flag = CTG_NODE_DEL;
  }

  return 0;
}

static void
ont_node_update (digraph_t * g, ont_node_t * ont_node)
{
  dg_node_t * node;
  ctg_vinfo_t * vinfo;

  node = digraph_node_at (g, ont_node->nid);
  vinfo = (ctg_vinfo_t *) node->info;
  vinfo->n_okmers += ont_node->n_okmers;
}

static void
ont_edge_create (ctg_graph_t * g, ont_node_t * from_ont_node, ont_node_t * to_ont_node)
{
  dg_edge_t * edge;
  ctg_einfo_t * einfo;
	ont_node_t * tmp;
	ont2gap_t * o2g;

  ont_node_update (g->g, from_ont_node);
  ont_node_update (g->g, to_ont_node);

  if (from_ont_node->direct == BACK) {
		tmp = from_ont_node;
		from_ont_node = to_ont_node;
		to_ont_node = tmp;
	}

	if ((edge = digraph_edge_between(g->g,from_ont_node->nid,to_ont_node->nid)) == NULL) {
    edge = digraph_add1edge (g->g, to_ont_node->nid, from_ont_node->nid);
		einfo = bmp_alloc (ctg_einfo, g->einfos);
		edge->info = einfo;
	} else {
		einfo = edge->info;
	}

	o2g = mp_alloc (o2g, einfo->ont_infos);
	o2g->ont_id = from_ont_node->oid;

	if (from_ont_node->direct == FORW) {
		o2g->left_okid = from_ont_node->end_okid;
		o2g->right_okid = to_ont_node->beg_okid;
	} else {
		o2g->left_okid = from_ont_node->beg_okid;
		o2g->right_okid = to_ont_node->end_okid;
	}
}

static int
add_ont_link2graph (ctg_graph_t * g, mp_t(ont_node) * ont_nodes, FILE * fp, FILE * raw_fp)
{
  int link_found;
  int direct;
  int32_t i;
  int32_t max_idx;
  int32_t max_cvg;
  int32_t n_ont_nodes;
  ont_node_t * ont_node;
  ont_node_t * pre_ont_node;

  n_ont_nodes = mp_cnt (ont_nodes);
  if (n_ont_nodes <= 1)
    return ONT_NO_LINK;

  // find the node with most okmers
  ont_node = mp_at (ont_node, ont_nodes, 0);
  fprintf (raw_fp, "%d\t%d\t%c\n",
        ont_node->tid, ont_node->n_okmers,
        ont_node->direct==FORW?'F':'B');
  max_idx = 0;
  max_cvg = ont_node->n_okmers;
  for (i=1; i<n_ont_nodes; ++i) {
    ont_node = mp_at (ont_node, ont_nodes, i);
    fprintf (raw_fp, "%d\t%d\t%c\n",
        ont_node->tid, ont_node->n_okmers,
        ont_node->direct==FORW?'F':'B');
    if (ont_node->n_okmers > max_cvg) {
      max_idx = i;
      max_cvg = ont_node->n_okmers;
    }
  }

  // extend and mark from the node with most okmers
  // backward
  pre_ont_node = mp_at (ont_node, ont_nodes, max_idx);
  direct = pre_ont_node->direct;
  for (i=max_idx-1; i>=0; --i) {
    ont_node = mp_at (ont_node, ont_nodes, i);
    if ((ont_node->direct != direct)
        || (direct==FORW && ont_node->nid!=pre_ont_node->nid-1 && ont_node->nid!=pre_ont_node->nid)
        || (direct==BACK && ont_node->nid!=pre_ont_node->nid+1 && ont_node->nid!=pre_ont_node->nid)) {
      if (ont_node->n_okmers < MIN_ANKOR_KMERS)
        ont_node->flag |= CTG_NODE_DEL;
      else
        return ONT_MISLEAD;
    } else
      pre_ont_node = ont_node;
  }
  // forward
  pre_ont_node = mp_at (ont_node, ont_nodes, max_idx);
  direct = pre_ont_node->direct;
  for (i=max_idx+1; i<n_ont_nodes; ++i) {
    ont_node = mp_at (ont_node, ont_nodes, i);
    if ((ont_node->direct != direct)
        || (direct==FORW && ont_node->nid!=pre_ont_node->nid+1 && ont_node->nid!=pre_ont_node->nid)
        || (direct==BACK && ont_node->nid!=pre_ont_node->nid-1 && ont_node->nid!=pre_ont_node->nid)) {
      if (ont_node->n_okmers < MIN_ANKOR_KMERS)
        ont_node->flag |= CTG_NODE_DEL;
      else
        return ONT_MISLEAD;
    } else
      pre_ont_node = ont_node;
  }

  /*
  // check ont anker direction
  for (i=0; i<mp_cnt(ont_nodes); ++i) {
    ont_node = mp_at (ont_node, ont_nodes, i);
    if (!(ont_node->flag & CTG_NODE_DEL))
      break;
  }
  if (i == mp_cnt(ont_nodes))
    return ONT_NO_ANK;
  pre_direct = ont_node->direct;
  for (++i; i<mp_cnt(ont_nodes); ++i) {
    ont_node = mp_at (ont_node, ont_nodes, i);
    if (ont_node->flag & CTG_NODE_DEL)
      continue;
    if (ont_node->direct != pre_direct)
      return ONT_MISLEAD;
  }
  */

  // merge ont_nodes belonging to the same contig
  for (i=0; i<mp_cnt(ont_nodes); ++i) {
    ont_node = mp_at (ont_node, ont_nodes, i);
    if (!(ont_node->flag & CTG_NODE_DEL))
      break;
  }
  pre_ont_node = ont_node;
  for (++i; i<mp_cnt(ont_nodes); ++i) {
    ont_node = mp_at (ont_node, ont_nodes, i);
    if (ont_node->flag & CTG_NODE_DEL)
      continue;

    if (ont_node->nid == pre_ont_node->nid) {
      ont_node->n_okmers += pre_ont_node->n_okmers;
			ont_node->beg_okid = pre_ont_node->beg_okid;
      pre_ont_node->flag |= CTG_NODE_DEL;
    } else {
      if ((ont_node->direct==FORW && ont_node->nid!=pre_ont_node->nid+1)
          || (ont_node->direct==BACK && ont_node->nid!=pre_ont_node->nid-1)) {
        return ONT_NEW_LINK;
      }
    }
    pre_ont_node = ont_node;
  }

  // add links to graph
  link_found = 0;
  for (i=0; i<mp_cnt(ont_nodes); ++i) {
    ont_node = mp_at (ont_node, ont_nodes, i);
    if (!(ont_node->flag & CTG_NODE_DEL))
      break;
  }
  pre_ont_node = ont_node;
  for (++i; i<mp_cnt(ont_nodes); ++i) {
    ont_node = mp_at (ont_node, ont_nodes, i);
    if (ont_node->flag & CTG_NODE_DEL)
      continue;
    link_found = 1;
    fprintf (fp, "%d:%c:%d\n", pre_ont_node->tid, pre_ont_node->direct==FORW?'F':'B', pre_ont_node->n_okmers);
    ont_edge_create (g, pre_ont_node, ont_node);
    pre_ont_node = ont_node;
  }
  fprintf (fp, "%d:%c:%d\n", pre_ont_node->tid, pre_ont_node->direct==FORW?'F':'B', pre_ont_node->n_okmers);

  if (!link_found)
    return ONT_NO_LINK;

  return 0;
}

static int
stat_ont_link (FILE * fp, ctg_graph_t * g, mp_t(ctg) * ctg_seqs,
    okseq_t * okseq, mp_t(ont_node) * ont_nodes)
{
  int i;
  int n_valid_nodes;
  int n_ont_nodes;
  ctg_t * ctg;
  ont_kmer_t * beg_okmer;
  ont_kmer_t * end_okmer;
  ont_node_t * ont_node;
  ont_node_t * pre_ont_node;

  n_ont_nodes = mp_cnt (ont_nodes);
  for (i=0; i<n_ont_nodes; ++i) {
    ont_node = mp_at (ont_node, ont_nodes, i);
    if (!(ont_node->flag & CTG_NODE_DEL))
      break;
  }
  pre_ont_node = ont_node;
  for (++i; i<n_ont_nodes; ++i) {
    ont_node = mp_at (ont_node, ont_nodes, i);
    if (ont_node->flag & CTG_NODE_DEL)
      continue;

    if (ont_node->nid == pre_ont_node->nid) {
      ont_node->n_okmers += pre_ont_node->n_okmers;
      ont_node->beg_okid = pre_ont_node->beg_okid;
      pre_ont_node->flag |= CTG_NODE_DEL;
    }
    pre_ont_node = ont_node;
  }

  for (i=0,n_valid_nodes=0; i<n_ont_nodes; ++i) {
    ont_node = mp_at (ont_node, ont_nodes, i);
    if (!(ont_node->flag & CTG_NODE_DEL))
      ++n_valid_nodes;
  }

  fprintf (fp, ">%ld\t%d\t0-%d\n", ont_node->oid, n_valid_nodes, okseq->seq->l);
  for (i=0; i<n_ont_nodes; ++i) {
    ont_node = mp_at (ont_node, ont_nodes, i);
    if (ont_node->flag & CTG_NODE_DEL)
      continue;
    ctg = mp_at (ctg, ctg_seqs, ont_node->tid);
    fprintf (fp, "%d:0-%d\t%d\t%c",
        ont_node->tid, ctg->seq->l, ont_node->n_okmers,
        ont_node->direct==FORW?'F':'B');

    beg_okmer = mp_at (okmer, okseq->okmers, ont_node->beg_okid);
    end_okmer = mp_at (okmer, okseq->okmers, ont_node->end_okid);

    fprintf (fp, "\t%d - %d", beg_okmer->ont_pos, end_okmer->ont_pos);
    fprintf (fp, "\t%d:%d - %d:%d\n", beg_okmer->kmer->tid, beg_okmer->kmer->pos,
        end_okmer->kmer->tid, end_okmer->kmer->pos);
  }
}

/*****************************************************************************/
/*----------------------------- Fix ONT Step 1 ------------------------------*/
/*****************************************************************************/

typedef struct {
  int32_t info_idx;
  int32_t len;
} ont_gap_t;
MP_DEF (ont_gap, ont_gap_t);

static int
cmp_ont_gap_len (const void * a, const void * b)
{
  ont_gap_t * pa = (ont_gap_t *) a;
  ont_gap_t * pb = (ont_gap_t *) b;

  if (pa->len < pb->len)
    return -1;
  else
    return 1;
}

static ctg_vinfo_t *
create_origin_N_gap_info (ctg_graph_t * g, int32_t ori_N_len)
{
  ctg_vinfo_t * vinfo;

  vinfo = bmp_alloc (ctg_vinfo, g->vinfos);
  vinfo->flag = CTG_GAP | CTG_N_GAP;
  vinfo->l_gap = ori_N_len;

  return vinfo;
}

static ctg_vinfo_t *
combine_ont_info (ctg_graph_t * g, ctg_einfo_t * einfo,
    mp_t(okseq) * okseqs, mp_t(ont_gap) * ont_gaps,
    ctg_t * from_ctg, ctg_t * to_ctg)
{
  char * src;
  int32_t i;
  int32_t idx;
  int32_t n_nega;
  int32_t n_onts;
  int32_t total_len;
  int32_t final_len;
  okseq_t * okseq;
  ont2gap_t * ont_info;
  ont_kmer_t * left_okmer;
  ont_kmer_t * right_okmer;
  ont_gap_t * ont_gap;
  ctg_vinfo_t * vinfo;

  n_onts = mp_cnt (einfo->ont_infos);
  mp_clear (ont_gap, ont_gaps, NULL);
  for (i=0; i<n_onts; ++i) {
    ont_info = mp_at (o2g, einfo->ont_infos, i);
    okseq = mp_at (okseq, okseqs, ont_info->ont_id);
    left_okmer = mp_at (okmer, okseq->okmers, ont_info->left_okid);
    right_okmer = mp_at (okmer, okseq->okmers, ont_info->right_okid);

    if (ont_info->left_okid < ont_info->right_okid)
      total_len = right_okmer->ont_pos - left_okmer->ont_pos;
    else
      total_len = left_okmer->ont_pos - right_okmer->ont_pos;

    ont_gap = mp_alloc (ont_gap, ont_gaps);
    ont_gap->info_idx = i;
    ont_gap->len = total_len - (from_ctg->seq->l-left_okmer->kmer->pos)
                        - right_okmer->kmer->pos;
  }

  qsort (ont_gaps->pool, n_onts, sizeof(ont_gap_t), cmp_ont_gap_len);

  idx = n_onts >> 1;
  ont_gap = mp_at (ont_gap, ont_gaps, idx);

  ont_info = mp_at (o2g, einfo->ont_infos, ont_gap->info_idx);
  okseq = mp_at (okseq, okseqs, ont_info->ont_id);
  left_okmer = mp_at (okmer, okseq->okmers, ont_info->left_okid);
  right_okmer = mp_at (okmer, okseq->okmers, ont_info->right_okid);

  vinfo = bmp_alloc (ctg_vinfo, g->vinfos);
  vinfo->flag = CTG_GAP;
  vinfo->l_gap = ont_gap->len;

  if (ont_gap->len <= 0) {
    vinfo->flag |= CTG_NEGA_GAP;
    return vinfo;
  }

  str_resize (vinfo->gap_seq, ont_gap->len);
  vinfo->gap_seq->l = ont_gap->len;

  // TODO: kmer length is better to be considered
  if (ont_info->left_okid < ont_info->right_okid) {
    src = okseq->seq->b + left_okmer->ont_pos
      + (from_ctg->seq->l - left_okmer->kmer->pos);
    memcpy (vinfo->gap_seq->s, src, ont_gap->len);
  } else {
    src = okseq->seq->b + left_okmer->ont_pos
      - (from_ctg->seq->l - left_okmer->kmer->pos);
    for (i=0; i<ont_gap->len; ++i,--src)
      vinfo->gap_seq->s[i] = base_rc_tbl[*src];
  }
  vinfo->gap_seq->s[ont_gap->len] = '\0';

  return vinfo;
}

static int
contig_seq_dump (FILE * fp, ctg_t * ctg, ctg_vinfo_t * gap_vinfo, int32_t * acu_len)
{
  int32_t i;
  int32_t len;

  len = *acu_len;

  if (gap_vinfo!=NULL && gap_vinfo->l_gap>0) {
    //printf (">pre_gap_size: %ld\n", gap_vinfo->l_gap);
    if (gap_vinfo->flag & CTG_N_GAP) {
      for (i=0; i<gap_vinfo->l_gap; ++i) {
        fprintf (fp, "N");
        if (++len % 60 == 0)
          fprintf (fp, "\n");
      }
    } else {
      for (i=0; i<gap_vinfo->l_gap; ++i) {
        fprintf (fp, "%c", gap_vinfo->gap_seq->s[i]);
        if (++len % 60 == 0)
          fprintf (fp, "\n");
      }
    }
  }

  //printf (">contig_size: %d\n", ctg->seq->l);
  for (i=0; i<ctg->seq->l; ++i) {
    fprintf (fp, "%c", ctg->seq->s[i]);
    if (++len % 60 == 0)
      fprintf (fp, "\n");
  }

  *acu_len = len;

  return 0;
}

/*****************************************************************************/
/*---------------------------- Extern Functions -----------------------------*/
/*****************************************************************************/

int
map_ont2contigs (mp_t(okseq) * okseqs, ctg_graph_t * g, mp_t(ctg) * ctg_seqs)
{
  int32_t j;
	int32_t beg;
  int32_t n_okmers;
  int32_t pre_ctg_id;
  int32_t direct[2];
  int64_t i;
  int64_t n_okseqs;
  FILE * fp;
  FILE * vfp;
  ctg_t * seq;
  okseq_t * okseq;
  ont_kmer_t * okmer;
  ont_kmer_t * pre_okmer;
  dg_node_t * node;
  ont_node_t * ont_node;
  mp_t(ont_node) * ont_nodes;

  fp = ckopen ("ont_link.txt", "w");
  vfp = ckopen ("valid_ont_link.txt", "w");

  // Load ONT Linkage
  n_okseqs = mp_cnt (okseqs);
  ont_nodes = mp_init (ont_node, NULL, NULL);
  for (i=0; i<n_okseqs; ++i) {
    okseq = mp_at (okseq, okseqs, i);
    n_okmers = mp_cnt (okseq->okmers);
    for (j=0; j<n_okmers; ++j) {
      okmer = mp_at (okmer, okseq->okmers, j);
      if (NULL != okmer->kmer)
        break;
    }
    if (j == n_okmers) {
      okseq->flag |= ONT_NO_ANK;
      continue;
    }

		beg = j;
    mp_clear (ont_node, ont_nodes, NULL);
    direct[FORW] = direct[BACK] = 0;
    pre_ctg_id = okmer->kmer->tid;
    if ((okmer->flag&ONT_KMER_REV) && (okmer->kmer->flag&KMER_REV)
        || !(okmer->flag&ONT_KMER_REV) && !(okmer->kmer->flag&KMER_REV)) {
      ++direct[FORW];
    } else {
      ++direct[BACK];
      okmer->flag |= ONT_SCAF_REV;
    }
    pre_okmer = okmer;

    //printf ("> ONT %d\n", i);
    for (++j; j<n_okmers; ++j) {
      okmer = mp_at (okmer, okseq->okmers, j);
      if (NULL == okmer->kmer)
        continue;
      if (okmer->kmer->tid != pre_ctg_id) {
        seq = mp_at (ctg, ctg_seqs, pre_ctg_id);
        ont_node = mp_alloc (ont_node, ont_nodes);
        ont_node_init (ont_node, seq->node_id, direct, pre_ctg_id, i, beg, pre_okmer->ont_pos, okseq);
        //printf ("%d -> ", pre_ctg_id);

				beg = j;
        pre_ctg_id = okmer->kmer->tid;
        direct[FORW] = direct[BACK] = 0;
      }

      if ((okmer->flag&ONT_KMER_REV) && (okmer->kmer->flag&KMER_REV)
          || !(okmer->flag&ONT_KMER_REV) && !(okmer->kmer->flag&KMER_REV)) {
        ++direct[FORW];
      } else {
        ++direct[BACK];
        okmer->flag |= ONT_SCAF_REV;
      }

      pre_okmer = okmer;
    }

    seq = mp_at (ctg, ctg_seqs, pre_ctg_id);
    ont_node = mp_alloc (ont_node, ont_nodes);
    ont_node_init (ont_node, seq->node_id, direct, pre_ctg_id, i, beg, pre_okmer->ont_pos, okseq);
    //printf ("%d\n", pre_ctg_id);

    //if (mp_cnt(ont_nodes) > 1) {
      fprintf (vfp, "> ONT %ld\n", i);
      fprintf (fp, "> ONT %ld\n", i);
      okseq->flag |= add_ont_link2graph (g,ont_nodes,vfp,fp);
    //}

    //stat_ont_link (fp, g, ctg_seqs, okseq, ont_nodes);
  }
  //printf ("\n");

  fclose (fp);
  fclose (vfp);

  return 0;
}

int
fix_ont1 (mp_t(okseq) * okseqs, ctg_graph_t * g,
    mp_t(ctg) * ctg_seqs, mp_t(sf) * scafs)
{
  char * ch;
	int32_t i, j, k, l;
	int32_t n_ctgs;
  int32_t n_sfs;
  int32_t raw_gap_size;
  int32_t acu_len;
  int32_t n_scaf;
  FILE * fp;
	ctg_t * ctg;
	ctg_t * f_ctg;
	ctg_t * t_ctg;
  scaf_t * sf;
	dg_edge_t * edge;
	ctg_einfo_t * einfo;
	ont2gap_t * ont_info;
	okseq_t * okseq;
	ont_kmer_t * left_okmer;
	ont_kmer_t * right_okmer;
  mp_t(ont_gap) * ont_gaps;
  ctg_vinfo_t * gap_vinfo;

  fp = ckopen ("gc_fix1.fa", "w");

  n_scaf = 0;
  ont_gaps = mp_init (ont_gap, NULL, NULL);
  n_sfs = mp_cnt (scafs);
	for (i=0; i<n_sfs; ++i) {
		sf = mp_at (sf, scafs, i);
		n_ctgs = mp_cnt (sf->ctg_ids);

    fprintf (fp, ">%d\n", n_scaf++);
    acu_len = 0;

    contig_seq_dump (fp, mp_at(ctg,ctg_seqs,sf->ctg_ids->arr[0]), NULL, &acu_len);

		for (j=1; j<n_ctgs; ++j) {
			f_ctg = mp_at (ctg, ctg_seqs, sf->ctg_ids->arr[j-1]);
			t_ctg = mp_at (ctg, ctg_seqs, sf->ctg_ids->arr[j]);

			if ((edge = digraph_edge_between(g->g,f_ctg->node_id,t_ctg->node_id)) == NULL)
				err_mesg ("no edge between %dth contig and %dth contig!", f_ctg->id, t_ctg->id);
			einfo = edge->info;

      if (mp_cnt(einfo->ont_infos) <= 0)
        gap_vinfo = create_origin_N_gap_info (g, t_ctg->l_pre_gap);
      else {
        // combine ONT info
        mp_resize (ont_gap, ont_gaps, mp_cnt(einfo->ont_infos));
        gap_vinfo = combine_ont_info (g, einfo, okseqs, ont_gaps, f_ctg, t_ctg);
      }

      contig_seq_dump (fp, t_ctg, gap_vinfo, &acu_len);
		}

    if (acu_len % 60 != 0)
      fprintf (fp, "\n");
	}

  mp_free (ont_gap, ont_gaps, NULL);

  /*
  n_ctgs = mp_cnt (ctg_seqs);
  for (i=0; i<n_ctgs; ++i) {
    ctg = mp_at (ctg, ctg_seqs, i);
    if (ctg->link_id >= 0)
      continue;
    acu_len = 0;
    fprintf (fp, "> %d\n", n_scaf++);
    contig_seq_dump (fp, ctg, NULL, &acu_len);
    if (acu_len % 60 != 0)
      fprintf (fp, "\n");
  }
  */

  fclose (fp);

  return 0;
}
