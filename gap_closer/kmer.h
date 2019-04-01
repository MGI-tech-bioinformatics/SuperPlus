/*************************************************
 * File Name: kmer.h
 * Description: 
 * Author: Chen Xi
 * Mail: chenxi1@genomics.cn
 * Created Time: 2018/12/09
 * Edit History: 
 *************************************************/

#ifndef KMER_H
#define KMER_H

#include "hash.h"
#include "kseq1.h"
#include "crc32.h"
#include "contig.h"

#define is_N(ch) (((ch)&0xf)==0xe ? 1 : 0)

static inline void kmer_copy_func (kmer_t * dst, kmer_t * src)
{
  memcpy (dst, src, sizeof(kmer_t));
}

static inline uint64_t kmer_hash_func (const void * key)
{
	kmer_t * k = (kmer_t *) key;

	return k->kseq;
}

static inline int kmer_is_equal (const void * a, const void * b)
{
	kmer_t * ka = (kmer_t *) a;
	kmer_t * kb = (kmer_t *) b;

	return ka->kseq == kb->kseq;
}

static inline void kmer_free2 (kmer_t * kmer)
{
  //if (NULL != kmer->ont_pos)
  //  mp_free (kpos, kmer->ont_pos, NULL);
}

static inline uint32_t kseq_crc32 (kseq1_t * kseq)
{
  return crc32 (0, (char*)kseq, sizeof(kseq1_t));
}

int chop_contig_seqs2kmers (mp_t(ctg) * seqs, int n_thread, int kmer_len);

int put_contig_kmers2hashs (xh_t ** khashs, mp_t(ctg) * seqs, int n_thread);

void kmer_set_dump (xh_t * khash, int kmer_len);

xh_t ** kmer_hash_init (int n_thread);
void kmer_hash_clear (xh_t ** khashs, int n_thread);
void kmer_hash_free (xh_t ** khashs, int n_thread);

void kmer_stat (xh_t ** ctg_khashs, int n_thread, int kmer_len, const char * mesg);

void kmer_stat2 (xh_set_t(kmer) ** khashs, int n_thread, int kmer_len, const char * mesg);

#endif
