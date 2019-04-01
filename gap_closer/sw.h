/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-05-16 17:08:52
  *Edit History: 
***********************************************************/

#ifndef XDK_SW_H
#define XDK_SW_H

#include <stdint.h>

#include "cigar.h"

#define XSW_STATUS "untested"

// SW alignment Overhang Strategy
#define SWOS_SOFTCLIP      0
#define SWOS_LEADING_INDEL 1
#define SWOS_INDEL         2
#define SWOS_IGNORE        3

struct sw_cell_s;
typedef struct sw_cell_s sw_cell_t;

struct sw_s;
typedef struct sw_s sw_t;

struct sw_cell_s {
	int32_t ms, is, ds;
	int32_t ml, dl, il;
	int32_t score;
	int status;
};

struct sw_s {
	sw_cell_t * sm; // score matrix
	int32_t sm_vol; // score matrix volumn
	int32_t m_qry;
	int32_t m_tgt;

	int32_t score;
	int32_t alignment_offset; // query alignment start on target
	int8_t status;

	int8_t type_c;
	int8_t qry_nbits;
	int8_t overhang_strategy;
	int64_t has_softclip;

	int32_t * mat;
	int32_t del_o;
	int32_t del_e;
	int32_t ins_o;
	int32_t ins_e;

	cigar_t * cigar;
};

#ifdef __cplusplus
extern "C" {
#endif

	sw_t * sw_init (void);

	void sw_set_parameter (sw_t * sw, int type_c, int32_t * mat,
			int32_t del_o, int32_t del_e, int32_t ins_o, int32_t ins_e, int overhang_strategy);

	int sw_align (sw_t * sw, int32_t qry_len, char * qry, int32_t tgt_len, char * tgt);

#ifdef __cplusplus
}
#endif

#endif
