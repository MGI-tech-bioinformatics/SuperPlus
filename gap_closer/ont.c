/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-12-10 17:23:39
  *Edit History: 
***********************************************************/

#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "mp.h"
#include "ont.h"
#include "rseq.h"
#include "utils.h"

#define THREAD_WAIT 0
#define THREAD_STOP 1
#define THREAD_SEARCH_ONT_KMERS 2
#define THREAD_REHASH_ONT_KMERS 3
#define THREAD_FIND_UNANKOR_POS 4

typedef struct {
  int * signal;
  int pid;

  int n_thread;
  int kmer_len;

  int64_t * cur_idx;
  int64_t * n_pth_done;

  mp_t(rs) * ont_seqs;
  mp_t(ctg) * ctg_seqs;
  mp_t(okseq) * okseqs;
  xh_t ** ctg_khashs;
  xh_set_t(kmer) ** anchored_ksets;
} work_para_t;

typedef struct {
  int n_thread;
  int * is_done;
  int64_t * n_pth_done;
  int64_t n_total;
  char * func_mesg;
  char * item_mesg;
} watch_para_t;

/*****************************************************************************/
/*------------------------ Shared Static Functions --------------------------*/
/*****************************************************************************/

static void *
watch_core (void * data)
{
  char * func_mesg;
  char * item_mesg;
  int i;
  int n_thread;
  int * is_done;
  int64_t n_total;
  int64_t n_total_done;
  int64_t * n_pth_done;
  double percent;
  watch_para_t * para;

  para = (watch_para_t *) data;
  n_thread = para->n_thread;
  is_done = para->is_done;
  n_total = para->n_total;
  n_pth_done = para->n_pth_done;
  func_mesg = para->func_mesg;
  item_mesg = para->item_mesg;

  for (;;) {
    for (i=0; i<n_thread; ++i)
      if (!is_done[i])
        break;
    if (i == n_thread)
      break;

    for (i=0,n_total_done=0; i<n_thread; ++i)
      n_total_done += n_pth_done[i];

    percent = 100.0 * ((double)n_total_done) / ((double)n_total);

    printf ("%s, %ld %s addressed, %.2lf%% percent have been done.\n",
        func_mesg, n_total_done, item_mesg, percent);

    sleep (10);
  }

  return (void *) 0;
}

static void
send_work_signal (int * signals, int signal, int n_thread)
{
  int i;

  for (i=0; i<n_thread; ++i)
    signals[i] = signal;

  for (;;) {
    for (i=0; i<n_thread; ++i)
      if (signals[i] != THREAD_WAIT)
        break;
    if (i == n_thread)
      break;

    usleep (100);
  }
}

static void
dump_okseq_info (okseq_t * okseq)
{
  int32_t i;
  int32_t n_okmers;
  ont_kmer_t * okmer;

  n_okmers = mp_cnt (okseq->okmers);
  for (i=0; i<n_okmers; ++i) {
    okmer = mp_at (okmer, okseq->okmers, i);
    printf ("ont_pos: %d\tctg_id: %d\tctg_pos: %d\n",
        okmer->ont_pos, okmer->kmer->tid, okmer->kmer->pos);
  }
  printf ("\n");
}

/*****************************************************************************/
/*---------------------------- Search ONT Kmers -----------------------------*/
/*****************************************************************************/

static void
chop1segm (char * bases, int32_t beg, int32_t end, okseq_t * okseq,
    int kmer_len, int n_thread, xh_t ** ctg_khashs, kseq1_t * kseq_mask)
{
  char bi;
  int ont_flag;
  int hs_id;
  int32_t i, j;
  int32_t l_segm;
  kseq1_t kseq;
  kseq1_t rv_kseq;
  kmer_t kmer;
  ont_kmer_t * ok;
  xh_item_t * item;

  if ((l_segm=end-beg) < kmer_len)
    return;

  kseq1_new_b (&kseq, bases, kmer_len);
  kseq1_fast_reverse_comp (&kseq, &rv_kseq, kmer_len);

  if (kseq1_cmp(&kseq,&rv_kseq) < 0) {
    memcpy (&kmer.kseq, &kseq, sizeof(kseq1_t));
    ont_flag = 0;
  } else {
    memcpy (&kmer.kseq, &rv_kseq, sizeof(kseq1_t));
    ont_flag = ONT_KMER_REV;
  }

  hs_id = kseq_crc32(&kmer.kseq) % n_thread;
  item = _xh_set_search3 (ctg_khashs[hs_id], (void*)(&kmer));
  if (item!=NULL && item->multi==1) {
    //ok = mp_alloc (okmer, okseq->okmers);
    ok = mp_at (okmer, okseq->okmers, beg);
    ok->kmer = (kmer_t *) item->key;
    ok->ont_pos = beg;
    ok->hs_id = hs_id;
    ok->flag = ont_flag;
  }

  for (i=beg+kmer_len,j=1; i<end; ++i,++j) {
    bi = base2int (bases[i]);
    kseq1_next (&kseq, bi, kseq_mask);
    kseq1_prev (&rv_kseq, int_comp(bi), kmer_len);

    if (kseq1_cmp(&kseq,&rv_kseq) < 0) {
      memcpy (&kmer.kseq, &kseq, sizeof(kseq1_t));
      ont_flag = 0;
    } else {
      memcpy (&kmer.kseq, &rv_kseq, sizeof(kseq1_t));
      ont_flag = ONT_KMER_REV;
    }

    hs_id = kseq_crc32(&kmer.kseq) % n_thread;
    item = _xh_set_search3 (ctg_khashs[hs_id], (void*)(&kmer));
    if (item!=NULL && item->multi==1) {
      //ok = mp_alloc (okmer, okseq->okmers);
      ok = mp_at (okmer, okseq->okmers, j);
      ok->kmer = (kmer_t *) item->key;
      ok->ont_pos = j;
      ok->hs_id = hs_id;
      ok->flag = ont_flag;
    }
  }
}

static int
chop1ont2kmer (rseq_t * ont_seq, okseq_t * okseq,
    int kmer_len, int n_thread, xh_t ** ctg_khashs, kseq1_t * kseq_mask)
{
  char * b;
  int32_t i;
  int32_t n_segs;
  ont_seg_t * seg;

  b = ont_seq->b;
  n_segs = mp_cnt (okseq->segs);
  for (i=0; i<n_segs; ++i) {
    seg = mp_at (oseg, okseq->segs, i);
    chop1segm (b, seg->beg, seg->end, okseq, kmer_len, n_thread, ctg_khashs, kseq_mask);
  }

	return 0;
}

/*****************************************************************************/
/*--------------------------- Re-hash ONT Kmers -----------------------------*/
/*****************************************************************************/

static int
rehash1okseq (okseq_t * okseq, xh_set_t(kmer) * anchored_kmers, int pid)
{
  int32_t i;
  int32_t n_okmers;
  ont_kmer_t * okmer;
  kmer_t * kmer;
  //kpos_t * kpos;

  n_okmers = mp_cnt (okseq->okmers);
  for (i=0; i<n_okmers; ++i) {
    okmer = mp_at (okmer, okseq->okmers, i);
    if (NULL == okmer->kmer)
      continue;
    if (okmer->hs_id != pid)
      continue;
    kmer = xh_set_add2 (kmer, anchored_kmers, okmer->kmer);
    //if (kmer->ont_pos == NULL)
    //  kmer->ont_pos = mp_init (kpos, NULL, NULL);
    //kpos = mp_alloc (kpos, kmer->ont_pos);
    //kpos->tid = i;
    //kpos->pos = okmer->ont_pos;
    //kpos->flag = okmer->flag;
    okmer->kmer = kmer;
  }  
}

/*****************************************************************************/
/*------------------------- Find Unankor Segments ---------------------------*/
/*****************************************************************************/

#define ONT_UNANKOR 0
#define ONT_ANKORED 1

static int
find_unankor_segs (okseq_t * okseq)
{
  int prev_stat;
  int curr_stat;
  int64_t i;
  int64_t n_okmers;
  ont_seg_t * seg;
  ont_kmer_t * okmer;

  mp_clear (oseg, okseq->segs, NULL);
  n_okmers = mp_cnt (okseq->okmers);

  prev_stat = ONT_ANKORED;
  for (i=0; i<n_okmers; ++i) {
    okmer = mp_at (okmer, okseq->okmers, i);
    if (NULL == okmer->kmer)
      curr_stat = ONT_UNANKOR;
    else
      curr_stat = ONT_ANKORED;

    switch (prev_stat) {
    case ONT_UNANKOR:
      switch (curr_stat) {
      case ONT_ANKORED:
        seg->end = i;
        break;
      }
      break;
    case ONT_ANKORED:
      switch (curr_stat) {
      case ONT_UNANKOR:
        seg = mp_alloc (oseg, okseq->segs);
        seg->beg = i;
        break;
      }
      break;
    }

    prev_stat = curr_stat;
  }

  if (prev_stat == ONT_UNANKOR)
    seg->end = i;

  return 0;
}

/*****************************************************************************/
/*---------------------------- Thread Routines ------------------------------*/
/*****************************************************************************/

static void *
thread_routines (void * data)
{
  int pid;
  int n_thread;
  int kmer_len;
  int * signal;
  int32_t i;
	int32_t n_ont_reads;
  int32_t n_contigs;
  int32_t n_okseqs;
  int32_t * kstat;
  int64_t pth_cur_idx;
  int64_t * cur_idx;
  int64_t * n_pth_done;
  rseq_t * r;
	kseq1_t kseq_mask;
  okseq_t * okseq;
  work_para_t * para;
  mp_t(rs) * ont_seqs;
  mp_t(ctg) * ctg_seqs;
  mp_t(okseq) * okseqs;
  xh_t ** ctg_khashs;
  xh_set_t(kmer) ** anchored_ksets;

  para = (work_para_t *) data;
  pid = para->pid;
  n_thread = para->n_thread;
  kmer_len = para->kmer_len;
  signal = para->signal;
  cur_idx = para->cur_idx;
  n_pth_done = para->n_pth_done;
  ont_seqs = para->ont_seqs;
  okseqs = para->okseqs;
  ctg_seqs = para->ctg_seqs;
  ctg_khashs = para->ctg_khashs;
  anchored_ksets = para->anchored_ksets;

	create_kseq1_mask (&kseq_mask, kmer_len);

  n_contigs = mp_cnt (ctg_seqs);
  kstat = (int32_t *) ckmalloc (n_contigs * sizeof(int32_t));

  for (;;) {
    if (*signal == THREAD_SEARCH_ONT_KMERS) {
	    n_ont_reads = mp_cnt (ont_seqs);
			for (;;) {
				pth_cur_idx = __sync_fetch_and_add (cur_idx, 1);
				if (pth_cur_idx >= n_ont_reads)
					break;

        okseq = mp_at (okseq, okseqs, pth_cur_idx);
				chop1ont2kmer (okseq->seq, okseq, kmer_len, n_thread, ctg_khashs, &kseq_mask);
			}

      *signal = THREAD_WAIT;
    } else if (*signal == THREAD_REHASH_ONT_KMERS) {
      n_okseqs = mp_cnt (okseqs);
      for (i=0; i<n_okseqs; ++i) {
        okseq = mp_at (okseq, okseqs, i);
        rehash1okseq (okseq, anchored_ksets[pid], pid);
      }

      *signal = THREAD_WAIT;
    } else if (*signal == THREAD_FIND_UNANKOR_POS) {
      n_okseqs = mp_cnt (okseqs);
      for (;;) {
        pth_cur_idx = __sync_fetch_and_add (cur_idx, 1);
        if (pth_cur_idx >= n_okseqs)
          break;

        okseq = mp_at (okseq, okseqs, pth_cur_idx);
        find_unankor_segs (okseq);
      }

      *signal = THREAD_WAIT;
    } else if (*signal == THREAD_STOP) {
      *signal = THREAD_WAIT;
      break;
    }
  }

  free (kstat);

  return (void *) 0;
}

/*****************************************************************************/
/*--------------------------- Extern Functions  -----------------------------*/
/*****************************************************************************/

int
search_kmers_on_ont_reads (mp_t(rs) * ont_seqs, mp_t(ctg) * ctg_seqs,
		xh_t ** ctg_khashs, mp_t(okseq) * okseqs, xh_set_t(kmer) ** anchored_ksets,
    const char * prefix, int n_thread, int kmer_len)
{
  int i;
  int * signals;
  time_t time_beg;
  time_t mod_tbeg;
  int64_t cur_idx;
  int64_t * n_pth_done;
  pthread_t * pids;
  work_para_t * paras;

  time (&time_beg);

  // create working threads
  pids = (pthread_t *) ckmalloc (n_thread * sizeof(pthread_t));
  paras = (work_para_t *) ckmalloc (n_thread * sizeof(work_para_t));
  signals = (int *) ckmalloc (n_thread * sizeof(int));
  n_pth_done = (int64_t *) ckmalloc (n_thread * sizeof(int64_t));
  for (i=0; i<n_thread; ++i) {
    signals[i] = THREAD_WAIT;
    paras[i].signal = signals + i;
    paras[i].pid = i;

    paras[i].n_thread = n_thread;
    paras[i].kmer_len = kmer_len;

    paras[i].cur_idx = &cur_idx;
    paras[i].n_pth_done = n_pth_done + i;

    paras[i].ont_seqs = ont_seqs;
    paras[i].ctg_seqs = ctg_seqs;
    paras[i].okseqs = okseqs;
    paras[i].ctg_khashs = ctg_khashs;
    paras[i].anchored_ksets = anchored_ksets;

    ckpthread_create (pids+i, NULL, thread_routines, (void*)(paras+i));
  }

  // search ont kmers
  time (&mod_tbeg);
  cur_idx = 0;
  memset (n_pth_done, 0, n_thread*sizeof(int64_t));
  send_work_signal (signals, THREAD_SEARCH_ONT_KMERS, n_thread);
  printf ("\n  chop and search ont kmers cost: %lds\n", time(NULL)-mod_tbeg);

  // re-hash ont kmers
  time (&mod_tbeg);
  memset (n_pth_done, 0, n_thread*sizeof(int64_t));
  send_work_signal (signals, THREAD_REHASH_ONT_KMERS, n_thread);
  printf ("\n  re-hash ont kmers cost: %lds\n", time(NULL)-mod_tbeg);

  // find un-ankored positions on onts
  time (&mod_tbeg);
  cur_idx = 0;
  memset (n_pth_done, 0, n_thread*sizeof(int64_t));
  send_work_signal (signals, THREAD_FIND_UNANKOR_POS, n_thread);
  printf ("\n  find un-ankored positions on onts costs: %lds\n", time(NULL)-mod_tbeg);

  send_work_signal (signals, THREAD_STOP, n_thread);

  for (i=0; i<n_thread; ++i)
    ckpthread_join (pids[i]);

	free (pids);
	free (paras);
	free (signals);
	free (n_pth_done);

  printf ("\n  search ONT kmers total cost: %lds\n", time(NULL)-time_beg);

  return 0;
}

int
ont_kseqs_init (mp_t(rs) * ont_seqs, mp_t(okseq) * okseqs)
{
  int64_t i;
  int64_t n_ont_seqs;
  rseq_t * r;
  okseq_t * okseq;
  ont_seg_t * seg;
  
  n_ont_seqs = mp_cnt (ont_seqs);
  mp_resize (okseq, okseqs, n_ont_seqs);
  okseqs->n = mp_cnt (ont_seqs);

  for (i=0; i<n_ont_seqs; ++i) {
    r = mp_at (rs, ont_seqs, i);
    okseq = mp_at (okseq, okseqs, i);

    okseq->seq = r;
    okseq->ont_id = i;

    mp_resize (okmer, okseq->okmers, r->l);
    okseq->okmers->n = r->l;

    mp_clear (oseg, okseq->segs, NULL);
    seg = mp_alloc (oseg, okseq->segs);
    seg->beg = 0;
    seg->end = r->l;
  }

  return 0;
}

int
ont_kseqs_dump (mp_t(okseq) * okseqs)
{
  int64_t i, j;
  int64_t n_okseqs;
  int64_t n_okmers;
  okseq_t * okseq;
  ont_kmer_t * okmer;

  n_okseqs = mp_cnt (okseqs);
  for (i=0; i<n_okseqs; ++i) {
    okseq = mp_at (okseq, okseqs, i);
    printf ("> ont %ld\n", i);
    n_okmers = mp_cnt (okseq->okmers);
    for (j=0; j<n_okmers; ++j) {
      okmer = mp_at (okmer, okseq->okmers, j);
      if (NULL == okmer->kmer)
        continue;
      printf ("%d:%d:%d\n", okmer->ont_pos, okmer->kmer->kmer_len, okmer->kmer->tid);
    }
  }

  return 0;
}
