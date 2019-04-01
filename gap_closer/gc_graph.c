/*************************************************
 * File Name: gc_graph.c
 * Description: 
 * Author: Chen Xi
 * Mail: chenxi1@genomics.cn
 * Created Time: 2018/12/20
 * Edit History: 
 *************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sw.h"
#include "kmer.h"
#include "cigar.h"
#include "digraph.h"
#include "gc_graph.h"

static inline void
gc_vinfo_init2 (gc_vinfo_t * vinfo, void * data)
{
	vinfo->s = str_init ();
	vinfo->flag = 0;
}

static inline void
gc_vinfo_clear (gc_vinfo_t * vinfo)
{
	str_clear (vinfo->s);
	vinfo->flag = 0;
}

static inline void
gc_vinfo_free2 (gc_vinfo_t * vinfo)
{
	str_free (vinfo->s);
}

static dg_node_t *
gc_node_create (gc_graph_t * g, const char * seq,
    int32_t beg, int32_t end, uint32_t flag, int32_t ctg_id, int32_t ctg_pos)
{
  int32_t l_seq;
  dg_node_t * node;
  gc_vinfo_t * vinfo;

  vinfo = bmp_alloc (gcv, g->vinfos);

  if ((l_seq=end-beg) <= 0)
    err_mesg ("seq len in vinfo <= 0!");

  str_resize (vinfo->s, l_seq);
  memcpy (vinfo->s->s, seq+beg, l_seq);
  vinfo->s->s[l_seq] = '\0';
  vinfo->s->l = l_seq;
  vinfo->flag = flag;

  vinfo->ctg_id = ctg_id;
  vinfo->ctg_pos = ctg_pos;

  node = digraph_add1node(g->g, vinfo);

  return node;
}

/*
#define SW_SCORE_M 1
#define SW_SCORE_X (-3)
#define SW_PENALTY_O 5
#define SW_PENALTY_E 2
*/

#define SW_SCORE_M 1
#define SW_SCORE_X (-5)
#define SW_PENALTY_O 2
#define SW_PENALTY_E 1

/*
#define SW_SCORE_M 50
#define SW_SCORE_X (-25)
#define SW_PENALTY_O 110
#define SW_PENALTY_E 6
*/

static sw_t *
gc_graph_sw_init (void)
{
  sw_t * sw;
  int i, j;
  int32_t mat[25];

  for (i=0; i<5; ++i) {
    for (j=0; j<5; ++j) {
      if (i == j)
        mat[5*i+j] = SW_SCORE_M;
      else
        mat[5*i+j] = SW_SCORE_X;
    }
  }

  sw = sw_init ();
  sw_set_parameter (sw, 5, mat, SW_PENALTY_O, SW_PENALTY_E, SW_PENALTY_O, SW_PENALTY_E, SWOS_SOFTCLIP);

  return sw;
}

/*****************************************************************************/
/*---------------------------- Extern Functions -----------------------------*/
/*****************************************************************************/

gc_graph_t *
gc_graph_init (void)
{
	gc_graph_t * g;

	g = (gc_graph_t *) ckmalloc (sizeof(gc_graph_t));
	g->g = digraph_init ();
	g->vinfos = bmp_init (gcv, 4096, gc_vinfo_init2, NULL);
  g->k2v_map = _xh_init (4096, 0.75, kmer_hash_func, kmer_is_equal);

  g->aligner = gc_graph_sw_init ();
	g->qry = str_init ();
	g->tgt = str_init ();

	return g;
}

void
gc_graph_clear (gc_graph_t * g)
{
	digraph_clear (g->g);
	bmp_clear (gcv, g->vinfos, gc_vinfo_clear);
  _xh_clear (g->k2v_map);
}

void
gc_graph_free (gc_graph_t * g)
{
	digraph_free (g->g);
	bmp_free (gcv, g->vinfos, gc_vinfo_free2);
  _xh_free (g->k2v_map);
  free (g);
}

void
gc_vinfo_dump (FILE * fp, void * info)
{
  gc_vinfo_t * vinfo;

  vinfo = (gc_vinfo_t *) info;

  fprintf (fp, "%d:%d\t%d:", vinfo->ctg_id, vinfo->ctg_pos, vinfo->s->l);

  if (vinfo->flag & GCV_ANC)
    fprintf (fp, "ANC|");
  else if (vinfo->flag & GCV_SEQ)
    fprintf (fp, "SEQ|");
  else if (vinfo->flag & GCV_GAP)
    fprintf (fp, "GAP|");

  if (vinfo->flag & GCV_SCF)
    fprintf (fp, "S");
  if (vinfo->flag & GCV_ONT)
    fprintf (fp, "O");
  if (vinfo->flag & GCV_LFR)
    fprintf (fp, "L");
}
