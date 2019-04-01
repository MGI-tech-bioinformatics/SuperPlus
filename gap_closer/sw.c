/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-05-16 17:21:53
  *Edit History: 
***********************************************************/

// TODO

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "sw.h"
#include "utils.h"
#include "cigar.h"

#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6
#define BAM_CEQUAL      7
#define BAM_CDIFF       8
#define BAM_CBACK       9

#define BAM_CIGAR_STR   "MIDNSHP=XB"
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  0xf
#define BAM_CIGAR_TYPE  0x3C1A7

#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
#define bam_cigar_opchr(c) (BAM_CIGAR_STR[bam_cigar_op(c)])
#define bam_cigar_gen(l, o) ((l)<<BAM_CIGAR_SHIFT|(o))
 
#define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3)

#define SW_M 1
#define SW_I 2
#define SW_D 4

#define INIT_QRY_LEN_MAX 128
#define INIT_TGT_LEN_MAX 512

#define INIT_QRY_NBITS 7

#define SW_UNDO 0
#define SW_DONE 1

/**********************************************************
 ******************* Static Functions *********************
 **********************************************************/

static void
init_matrix_values (sw_t * sw)
{
	int32_t i;
	int32_t cur_score;

	sw->sm[0].is = INT32_MIN / 2;
	sw->sm[0].il = 0;
	sw->sm[0].ds = INT32_MIN / 2;
	sw->sm[0].dl = 0;

	sw->sm[1].is = INT32_MIN / 2;
	sw->sm[1].il = 0;
	sw->sm[1].ds = INT32_MIN / 2;
	sw->sm[1].dl = 0;

	for (i=2; i<sw->m_qry; i++) {
		sw->sm[i].is = INT32_MIN / 2;
		sw->sm[i].il = 0;
		sw->sm[i].ds = INT32_MIN / 2;
		sw->sm[i].dl = 0;
	}

	sw->sm[sw->m_qry].is = INT32_MIN / 2;
	sw->sm[sw->m_qry].il = 0;
	sw->sm[sw->m_qry].ds = INT32_MIN / 2;
	sw->sm[sw->m_qry].dl = 0;
	for (i=2; i<sw->m_tgt; i++) {
		sw->sm[i<<sw->qry_nbits].is = INT32_MIN / 2;
		sw->sm[i<<sw->qry_nbits].il = 0;
		sw->sm[i<<sw->qry_nbits].ds = INT32_MIN / 2;
		sw->sm[i<<sw->qry_nbits].dl = 0;
	}

	if (sw->overhang_strategy == SWOS_SOFTCLIP)
		return;

	sw->sm[1].score = -sw->ins_o;
	cur_score = -sw->ins_o;
	for (i=2; i<sw->m_qry; i++) {
		cur_score -= sw->ins_e;
		sw->sm[i].score = cur_score;
	}

	sw->sm[sw->m_qry].score = -sw->del_o;
	cur_score = -sw->del_o;
	for (i=2; i<sw->m_tgt; i++) {
		cur_score -= sw->del_e;
		sw->sm[i<<sw->qry_nbits].score = cur_score;
	}
}

static int
matrix_dump (sw_t * sw, char * qry, int32_t l_qry, char * tgt, int32_t l_tgt)
{
  sw_cell_t * c;
  sw_cell_t * curr_row;
  int32_t i, j;

  printf ("\n");
  for (i=0; i<l_tgt+1; ++i) {
    curr_row = sw->sm + (i << sw->qry_nbits);
    printf ("%d:\t", i);
    for (j=0; j<l_qry+1; ++j) {
      c = curr_row + j;
      printf ("%d\t", c->score);
    }
    printf ("\n");
  }

  return 0;
}

static int
score_matrix_init (sw_t * sw, int32_t qry_len, int32_t tgt_len)
{
	int is_changed;

	is_changed = 0;

	if (sw->m_qry < qry_len+2) {
		is_changed = 1;
		while (sw->m_qry < qry_len+2) {
			sw->m_qry <<= 1;
			sw->qry_nbits += 1;
		}
	}

	if (sw->m_tgt < tgt_len+2) {
		is_changed = 1;
		while (sw->m_tgt < tgt_len+2)
			sw->m_tgt <<= 1;
	}

	if (is_changed) {
		sw->sm_vol = sw->m_tgt << sw->qry_nbits;
		free (sw->sm);
		sw->sm = (sw_cell_t *) ckalloc (sw->sm_vol, sizeof(sw_cell_t));
		init_matrix_values (sw);
	}

	return 0;
}

static int
align_core (sw_t * sw, int32_t qry_len, char * qry, int32_t tgt_len, char * tgt)
{
	char t, q;
	int type_c;
	int cigar_inited;
	int32_t i, j;
	int32_t ki, kd;
	int32_t m, del, ins;
	int32_t bt_qidx;
	int32_t bt_tidx;
	int32_t cur_score;
	int32_t max_score;
	int32_t pre_gap;
	int32_t del_o;
	int32_t del_e;
	int32_t ins_o;
	int32_t ins_e;
	int32_t * mat;
	uint32_t seg_len;
	uint32_t step_len;
	uint32_t pre_opr;
	uint32_t cur_opr;
	sw_cell_t * c;
	sw_cell_t * uc;
	sw_cell_t * lc;
	sw_cell_t * ulc;
	sw_cell_t * sm;
	sw_cell_t * pre_row;
	sw_cell_t * cur_row;

	type_c = sw->type_c;
	mat    = sw->mat;
	del_o  = sw->del_o;
	del_e  = sw->del_e;
	ins_o  = sw->ins_o;
	ins_e  = sw->ins_e;
	sm     = sw->sm;

	pre_row = sm;
	for (i=1; i<tgt_len+1; i++) {
		t = tgt[i-1];
		cur_row = pre_row + sw->m_qry;
		for (j=1; j<qry_len+1; j++) {
			q = qry[j-1];

			c   = cur_row + j;
			lc  = cur_row + j-1;
			uc  = pre_row + j;
			ulc = pre_row + j-1;

			// match
			c->ms = ulc->score + mat[q*type_c+t];

			// deletion
			pre_gap = uc->score - del_o;
			c->ds = uc->ds - del_e;
			if (pre_gap > c->ds) {
				c->ds = pre_gap;
				c->dl = 1;
			} else
				c->dl = uc->dl + 1;

			// insertion
			pre_gap = lc->score - ins_o;
			c->is = lc->is - ins_e;
			if (pre_gap > c->is) {
				c->is = pre_gap;
				c->il = 1;
			} else
				c->il = lc->il + 1;

			// comparison
			if (c->ms>=c->ds && c->ms>=c->is) {
				c->status = SW_M;
				c->score  = c->ms;
				c->ml     = ulc->ml + 1;
			} else if (c->is > c->ds) {
				c->status = SW_I;
				c->score  = c->is;
				c->ml     = 0;
			} else {
				c->status = SW_D;
				c->score  = c->ds;
				c->ml     = 0;
			}
		}
		pre_row = cur_row;
	}

	// backtrace
	bt_tidx = tgt_len;
	bt_qidx = qry_len;
	seg_len = 0;

	max_score = INT32_MIN;
	for (i=1; i<tgt_len+1; i++) {
		cur_score = sm[(i<<sw->qry_nbits)+qry_len].score;
		if (cur_score >= max_score) {
			bt_tidx = i;
			max_score = cur_score;
		}
	}

	if (sw->overhang_strategy != SWOS_LEADING_INDEL) {
		cur_row = sm + (tgt_len<<sw->qry_nbits);
		for (i=1; i<qry_len+1; i++) {
			cur_score = cur_row[i].score;
			if (cur_score > max_score
					|| (cur_score==max_score && abs(tgt_len-i)<abs(bt_tidx-bt_qidx))) {
				bt_tidx = tgt_len;
				bt_qidx = i;
				max_score = cur_score;
				seg_len = qry_len - i;
			}
		}
	}

	sw->has_softclip = 0;
	if (seg_len>0 && sw->overhang_strategy==SWOS_SOFTCLIP) {
		cigar_add (sw->cigar, (seg_len<<BAM_CIGAR_SHIFT)|BAM_CSOFT_CLIP);
		seg_len = 0;
		sw->has_softclip = 1;
	}

	c = sm + (bt_tidx<<sw->qry_nbits) + bt_qidx;
	sw->score = c->score;

	cigar_inited = 0;
	do {
		if (c->status == SW_M) {
			step_len = c->ml;
			cur_opr  = BAM_CMATCH;
			bt_tidx -= step_len;
			bt_qidx -= step_len;
		} else if (c->status == SW_D) {
			step_len = c->dl;
			cur_opr  = BAM_CDEL;
			bt_tidx -= step_len;
		} else {
			step_len = c->il;
			cur_opr  = BAM_CINS;
			bt_qidx -= step_len;
		}

		if (cigar_inited && cur_opr!=pre_opr) {
			cigar_add (sw->cigar, (seg_len<<BAM_CIGAR_SHIFT)|pre_opr);
			seg_len = 0;
		}

		seg_len += step_len;
		pre_opr  = cur_opr;

		if (!cigar_inited)
			cigar_inited = 1;
	} while (bt_tidx>0 && bt_qidx>0);

	cigar_add (sw->cigar, (seg_len<<BAM_CIGAR_SHIFT)|pre_opr);

	if (sw->overhang_strategy == SWOS_SOFTCLIP) {
		if (bt_qidx > 0)
			cigar_add (sw->cigar, (bt_qidx<<BAM_CIGAR_SHIFT)|BAM_CSOFT_CLIP);
		sw->alignment_offset = bt_tidx;
	} else {
		if (bt_tidx > 0)
			cigar_add (sw->cigar, (bt_tidx<<BAM_CIGAR_SHIFT)|BAM_CDEL);
		if (bt_qidx > 0)
			cigar_add (sw->cigar, (bt_qidx<<BAM_CIGAR_SHIFT)|BAM_CINS);
		sw->alignment_offset = 0;
	}

	cigar_reverse (sw->cigar);

	return 0;
}

/**********************************************************
 ******************* Extern Functions *********************
 **********************************************************/

sw_t *
sw_init (void)
{
	sw_t * sw;

	sw = (sw_t *) ckalloc (1, sizeof(sw_t));
	sw->sm_vol = INIT_TGT_LEN_MAX * INIT_QRY_LEN_MAX;
	sw->sm = (sw_cell_t *) ckalloc (sw->sm_vol, sizeof(sw_cell_t));

	sw->m_qry = INIT_QRY_LEN_MAX;
	sw->m_tgt = INIT_TGT_LEN_MAX;

	sw->cigar = cigar_init ();

	sw->type_c = -1;
	sw->mat = NULL;

	sw->qry_nbits = INIT_QRY_NBITS;

	sw->overhang_strategy = SWOS_SOFTCLIP;

	return sw;
}

void
sw_set_parameter (sw_t * sw, int type_c, int32_t * mat,
		int32_t del_o, int32_t del_e, int32_t ins_o, int32_t ins_e, int overhang_strategy)
{
	int32_t mat_nbytes;

	if (type_c <= 0)
		err_mesg ("[%s] invalid 'type_c'!", __func__);

	if (mat == NULL)
		err_mesg ("[%s] mat==NULL!", __func__);

	mat_nbytes = type_c * type_c * 4;
	if (sw->type_c != type_c) {
		if (sw->type_c > 0)
			free (sw->mat);
		sw->type_c = type_c;
		sw->mat = (int32_t *) ckmalloc (mat_nbytes);
	}
	memcpy (sw->mat, mat, mat_nbytes);

	sw->del_o = del_o;
	sw->del_e = del_e;
	sw->ins_o = ins_o;
	sw->ins_e = ins_e;

	sw->overhang_strategy = overhang_strategy;

	init_matrix_values (sw);
}

int
sw_align (sw_t * sw, int32_t qry_len, char * qry, int32_t tgt_len, char * tgt)
{
	sw->status = SW_UNDO;

	score_matrix_init (sw, qry_len, tgt_len);

	sw->score = -1;
	sw->alignment_offset = -1;
	cigar_clear (sw->cigar);

	align_core (sw, qry_len, qry, tgt_len, tgt);

	sw->status = SW_DONE;

	return 0;
}
