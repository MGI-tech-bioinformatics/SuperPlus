/*************************************************
 * File Name: kmer.c
 * Description: 
 * Author: Chen Xi
 * Mail: chenxi1@genomics.cn
 * Created Time: 2018/12/09
 * Edit History: 
 *************************************************/

#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "mp.h"
#include "str.h"
#include "hash.h"
#include "crc32.h"
#include "utils.h"
#include "contig.h"
#include "kmer.h"

typedef struct {
  int kmer_len;
  int n_thread;
	int32_t * cur_idx;
  mp_t(ctg) * seqs;
} chop_para_t;

typedef struct {
  int hs_id;
  mp_t(ctg) * seqs;
  xh_t * khash;
} hash_para_t;

static void *
chop_kmer_core (void * data)
{
	char * s;
  int kmer_len;
  int n_thread;
  int32_t i, j;
  int32_t kend;
	int32_t l_seq;
  int32_t n_seqs;
  int32_t pth_cur_idx;
  int32_t * cur_idx;
	uint64_t bi; // base index
	kseq1_t kseq;
	kseq1_t rv_kseq;
  kseq1_t kseq_mask;
	kmer_t * sk;
  ctg_t * seq;
  chop_para_t * para;
  mp_t(ctg) * seqs;

  para = (chop_para_t *) data;
  cur_idx = para->cur_idx;
  kmer_len = para->kmer_len;
  n_thread = para->n_thread;
  seqs = para->seqs;
  n_seqs = mp_cnt (seqs);
  create_kseq1_mask (&kseq_mask, kmer_len);

  for (;;) {
    pth_cur_idx = __sync_fetch_and_add (cur_idx, 1);
    if (pth_cur_idx >= n_seqs)
      break;
    seq = mp_at (ctg, seqs, pth_cur_idx);

		s = seq->seq->s;
    kend = seq->seq->l - kmer_len + 1;
    if (kend < 1)
      continue;

    kseq1_new_b (&kseq, s, kmer_len);
    kseq1_fast_reverse_comp (&kseq, &rv_kseq, kmer_len);

    sk = seq->kmers + seq->n_kmer++;
    sk->tid = pth_cur_idx;
    sk->pos = 0;
    sk->kmer_len = kmer_len;
    //if (NULL != sk->ont_pos)
    //  mp_clear (kpos, sk->ont_pos, NULL);

    if (kseq1_cmp(&kseq,&rv_kseq) < 0) {
      memcpy (&sk->kseq, &kseq, sizeof(kseq1_t));
      sk->hs_id = kseq_crc32(&kseq) % n_thread;
      sk->flag = 0;
    } else {
      memcpy (&sk->kseq, &rv_kseq, sizeof(kseq1_t));
      sk->hs_id = kseq_crc32(&rv_kseq) % n_thread;
      sk->flag = KMER_REV;
    }

    for (i=1,j=kmer_len; i<kend; ++i,++j) {
      bi = base2int (s[j]);
      kseq1_next (&kseq, bi, &kseq_mask);
      kseq1_prev (&rv_kseq, int_comp(bi), kmer_len);

      sk = seq->kmers + seq->n_kmer++;
      sk->tid = pth_cur_idx;
      sk->pos = i;
      sk->kmer_len = kmer_len;
      //if (NULL != sk->ont_pos)
      //  mp_clear (kpos, sk->ont_pos, NULL);

      if (kseq1_cmp(&kseq,&rv_kseq) < 0) {
        memcpy (&sk->kseq, &kseq, sizeof(kseq1_t));
        sk->hs_id = kseq_crc32(&kseq) % n_thread;
        sk->flag = 0;
      } else {
        memcpy (&sk->kseq, &rv_kseq, sizeof(kseq1_t));
        sk->hs_id = kseq_crc32(&rv_kseq) % n_thread;
        sk->flag = KMER_REV;
      }
    }
  }

  return (void *) 0;
}

static void *
hash_kmer_core (void * data)
{
  int hs_id;
  int32_t i, j;
  int32_t n_seq;
  ctg_t * seq;
  kmer_t * onk;
  hash_para_t * para;
  xh_t * khash;
  mp_t(ctg) * seqs;

  para = (hash_para_t *) data;
  hs_id = para->hs_id;
  khash = para->khash;
  seqs = para->seqs;

  n_seq = mp_cnt (seqs);

  for (i=0; i<n_seq; ++i) {
    seq = mp_at (ctg, seqs, i);
    for (j=0; j<seq->n_kmer; ++j) {
      onk = seq->kmers + j;
      if (onk->hs_id == hs_id)
  			_xh_set_add (khash, (void*)(seq->kmers+j));
    }
  }

  return (void *) 0;
}

int
chop_contig_seqs2kmers (mp_t(ctg) * seqs, int n_thread, int kmer_len)
{
  int i;
  time_t time_beg;
  int32_t cur_idx;
	pthread_t * pids;
  chop_para_t * paras;

  cur_idx = 0;
  time (&time_beg);
  pids = (pthread_t *) ckalloc (n_thread, sizeof(pthread_t));
  paras = (chop_para_t *) ckalloc (n_thread, sizeof(chop_para_t));
  for (i=0; i<n_thread; ++i) {
    paras[i].cur_idx = &cur_idx;
    paras[i].kmer_len = kmer_len;
    paras[i].n_thread = n_thread;
    paras[i].seqs = seqs;
    ckpthread_create (pids+i, NULL, chop_kmer_core, (void*)(paras+i));
  }

  for (i=0; i<n_thread; ++i)
    ckpthread_join (pids[i]);

  fprintf (stdout, "  chop kmers cost: %lds\n", time(NULL)-time_beg);

  free (pids);
  free (paras);

  return 0;
}

int
put_contig_kmers2hashs (xh_t ** khashs, mp_t(ctg) * seqs, int n_thread)
{
  int i;
  time_t time_beg;
  pthread_t * pids;
  hash_para_t * paras;

  time (&time_beg);
  pids = (pthread_t *) ckalloc (n_thread, sizeof(pthread_t));
  paras = (hash_para_t *) ckalloc (n_thread, sizeof(hash_para_t));
  for (i=0; i<n_thread; ++i) {
    paras[i].hs_id = i;
    paras[i].seqs = seqs;
    paras[i].khash = khashs[i];
    ckpthread_create (pids+i, NULL, hash_kmer_core, (void*)(paras+i));
  }

  for (i=0; i<n_thread; ++i)
    ckpthread_join (pids[i]);

  fprintf (stdout, "  hash kmers cost: %lds\n", time(NULL)-time_beg);

  free (pids);
  free (paras);

  return 0;
}

void
kmer_set_dump (xh_t * khash, int kmer_len)
{
  char * kmer_seq;
  int32_t i;

  kmer_seq = ALLOC_LINE;
  for (i=0; i<khash->cnt; ++i) {
    kmer_t * k = (kmer_t*) khash->pool[i].key;
    uint32_t multi = khash->pool[i].multi;
    kseq12seq (&k->kseq, kmer_seq, kmer_len);
    printf (">%lx\n%s\t%c\n", k->kseq, kmer_seq, k->flag&KMER_REV?'R':'F');
  }
  free (kmer_seq);
}

xh_t **
kmer_hash_init (int n_thread)
{
  int i;
  xh_t ** khashs;

  khashs = (xh_t **) ckalloc (n_thread, sizeof(xh_t*));
  for (i=0; i<n_thread; ++i)
    khashs[i] = _xh_init (65536, 0.75, kmer_hash_func, kmer_is_equal);

  return khashs;
}

void
kmer_hash_clear (xh_t ** khashs, int n_thread)
{
  int i;

  for (i=0; i<n_thread; ++i)
    _xh_clear (khashs[i]);
}

void
kmer_hash_free (xh_t ** khashs, int n_thread)
{
  int i;

  for (i=0; i<n_thread; ++i)
    _xh_free (khashs[i]);

  free (khashs);
}

void
kmer_stat (xh_t ** khashs, int n_thread, int kmer_len, const char * mesg)
{
  int64_t i, j;
  int64_t n_total_kmers;
  int64_t n_unique_kmers;
  xh_t * h;
  xh_item_t * item;

  n_total_kmers = n_unique_kmers = 0;
  for (i=0; i<n_thread; ++i) {
    h = khashs[i];
    n_total_kmers += h->cnt;
    for (j=0; j<h->cnt; ++j) {
      item = h->pool + j;
      if (item->multi == 1)
        ++n_unique_kmers;
    }
  }

  printf ("\n  >>> %s Kmer%d Stat <<<\n", mesg, kmer_len);
  printf ("  total kmer count: %ld\n", n_total_kmers);
  printf ("  unique kmer count: %ld\n", n_unique_kmers);
}

void
kmer_stat2 (xh_set_t(kmer) ** khashs, int n_thread, int kmer_len, const char * mesg)
{
  int64_t i, j;
  int64_t n_total_kmers;
  int64_t n_unique_kmers;
  xh_t * h;
  xh_item_t * item;

  n_total_kmers = n_unique_kmers = 0;
  for (i=0; i<n_thread; ++i) {
    h = khashs[i]->hash;
    n_total_kmers += h->cnt;
    for (j=0; j<h->cnt; ++j) {
      item = h->pool + j;
      if (item->multi == 1)
        ++n_unique_kmers;
    }
  }

  printf ("\n  >>> %s Kmer%d Stat <<<\n", mesg, kmer_len);
  printf ("  total kmer count: %ld\n", n_total_kmers);
  printf ("  unique kmer count: %ld\n", n_unique_kmers);
}
