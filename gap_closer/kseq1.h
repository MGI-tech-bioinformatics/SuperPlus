/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-12-19 18:08:56
  *Edit History: 
***********************************************************/

#ifndef KSEQ1_H
#define KSEQ1_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "bio.h"
#include "utils.h"

typedef uint64_t kseq1_t;

static inline void kseq1_plus (kseq1_t * kseq, uint64_t base)
{
	*kseq <<= 2;
	*kseq |= base;
}

static inline void kseq1_new_b (kseq1_t * kseq, char * s, int kmer_len)
{
	int i;

	*kseq = 0;
	for (i=0; i<kmer_len; ++i)
		kseq1_plus (kseq, (uint64_t)base2int(s[i]));
}

static inline void kseq1_fast_reverse_comp (kseq1_t * kseq, kseq1_t * rv_kseq, int kmer_len)
{
	*rv_kseq = *kseq ^ 0xAAAAAAAAAAAAAAAALLU;
	*rv_kseq = ((*rv_kseq & 0x3333333333333333LLU)<<2)  | ((*rv_kseq & 0xCCCCCCCCCCCCCCCCLLU)>>2);
	*rv_kseq = ((*rv_kseq & 0x0F0F0F0F0F0F0F0FLLU)<<4)  | ((*rv_kseq & 0xF0F0F0F0F0F0F0F0LLU)>>4);
	*rv_kseq = ((*rv_kseq & 0x00FF00FF00FF00FFLLU)<<8)  | ((*rv_kseq & 0xFF00FF00FF00FF00LLU)>>8);
	*rv_kseq = ((*rv_kseq & 0x0000FFFF0000FFFFLLU)<<16) | ((*rv_kseq & 0xFFFF0000FFFF0000LLU)>>16);
	*rv_kseq = ((*rv_kseq & 0x00000000FFFFFFFFLLU)<<32) | ((*rv_kseq & 0xFFFFFFFF00000000LLU)>>32);
	*rv_kseq >>= (64-(kmer_len<<1));
}

static inline void kseq1_prev (kseq1_t * kseq, uint64_t base, int kmer_len)
{
	*kseq >>= 2;
	*kseq |= (base << ((kmer_len-1)<<1));
}

static inline void kseq1_next (kseq1_t * kseq, uint64_t base, kseq1_t * kseq_mask)
{
	*kseq <<= 2;
	*kseq &= *kseq_mask;
	*kseq |= base;
}

static inline uint8_t last_char_in_kseq1 (kseq1_t * kseq)
{
	return (uint8_t) (*kseq & 0x03);
}

static inline uint8_t first_char_in_kseq1 (kseq1_t * kseq, int kmer_len)
{
	return (uint8_t) (*kseq >> ((kmer_len-1)<<1) & 0x03);
}

static inline void create_kseq1_mask (kseq1_t * kseq, int kmer_len)
{
	*kseq = (((uint64_t)1) << (kmer_len<<1)) - 1;
}

static inline uint64_t cal_kseq1_os (kseq1_t * kseq, uint64_t size)
{
	return *kseq % size;
}

static inline void kseq12seq (kseq1_t * kseq, char * seq, int kmer_len)
{
	int i;
	uint64_t tmp;

	tmp = *kseq;
	for (i=kmer_len-1; i>=0; --i) {
		seq[i] = "ACTG"[tmp & 0x03];
		tmp >>= 2;
	}
  seq[kmer_len] = '\0';
}

static inline void kseq12seq_i (kseq1_t * kseq, char * seq, int kmer_len)
{
	int i;
	uint64_t tmp;

	tmp = *kseq;
	for (i=kmer_len-1; i>=0; --i) {
		seq[i] = tmp & 0x03;
		tmp >>= 2;
	}
}

static inline void kseq12seq_r (kseq1_t * kseq, char * seq, int kmer_len)
{
	int i;
	uint64_t tmp;

	tmp = *kseq;
	for (i=kmer_len-1; i>=0; --i) {
		seq[kmer_len-1-i] = "ACTG"[tmp & 0x03];
		tmp >>= 2;
	}
  seq[kmer_len] = '\0';
}

static inline void kseq12seq_ri (kseq1_t * kseq, char * seq, int kmer_len)
{
	int i;
	uint64_t tmp;

	tmp = *kseq;
	for (i=kmer_len-1; i>=0; --i) {
		seq[kmer_len-1-i] = tmp & 0x03;
		tmp >>= 2;
	}
}

static inline int kseq1_cmp (kseq1_t * kseq1, kseq1_t * kseq2)
{
	if (*kseq1 < *kseq2)
		return -1;
	else if (*kseq1 > *kseq2)
		return 1;
	else
		return 0;
}

static inline int kseq1_equal (kseq1_t * kseq1, kseq1_t * kseq2)
{
	if (*kseq1 != *kseq2)
		return 0;
	else
		return 1;
}

#endif
