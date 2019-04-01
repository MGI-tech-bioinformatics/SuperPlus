/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-01-18 14:32:43
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "mp.h"
#include "def.h"
#include "lfr.h"
#include "str.h"
#include "hash.h"
#include "kmer.h"
#include "utils.h"

#define THREAD_WAIT 0
#define THREAD_STOP 1
#define THREAD_CALC_BC 2
#define THREAD_HASH_BC 3
#define THREAD_CLUSTER 4

HASH_SET_DEF (bc, uint64_t);

typedef struct {
  int pid;
  int n_thread;
  int kmer_len;
  int * signal;
  int64_t * cur_idx;

  mp_t(ctg) * ctgs;
  xh_t ** ctg_khashs;
  lfr_set_t * lfr_seqs;

  xh_set_t(bc) ** bc_hashs;
} work_para_t;

/*****************************************************************************/
/*------------------------- Shared Static Functions -------------------------*/
/*****************************************************************************/

static uint64_t
barcode_hash_func (const void * data)
{
  return *((uint64_t*)data);
}

static int
barcode_is_equal (const void * a, const void * b)
{
  return !memcmp (a,b,sizeof(uint64_t));
}

static void
lfr_barcode_init2 (lfr_bc_t * bc, void * data)
{
  bc->reads = mp_init (lfr, NULL, NULL);
}

static void
lfr_barcode_free2 (lfr_bc_t * bc)
{
  mp_free (lfr, bc->reads, NULL);
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
lfr_stat (lfr_set_t * set)
{
  int64_t i;
  lfr_bc_t * lfr_bc;

  printf ("\n");
  printf ("Total LFR reads: %ld\n", mp_cnt(set->lfrs));
  printf ("Total PE reads: %ld\n", mp_cnt(set->pes));
  printf ("Barcode count: %ld\n", mp_cnt(set->bcs));
  for (i=0; i<mp_cnt(set->bcs); ++i) {
    lfr_bc = mp_at (lbc, set->bcs, i);
    printf ("%ldth barcode, %ld pe reads\n", i, mp_cnt(lfr_bc->reads));
  }
}

/*****************************************************************************/
/*--------------------------- Barcode Calculation ---------------------------*/
/*****************************************************************************/

static uint64_t
read_name2barcode (char * read_name, int32_t len)
{
  char * ch;
  int i;
  uint64_t bc;
  uint64_t num;

  for (i=0,ch=read_name; i<len; ++i,++ch)
    if (*ch == '#')
      break;
  if (i >= len)
    err_mesg ("invalid read name '%s'!", read_name);

  bc = 0;
  num = atoi (ch+1);
  bc |= num << 40;

  for (++i,++ch; i<len; ++i,++ch)
    if (*ch == '_')
      break;
  if (i >= len)
    err_mesg ("invalid read name '%s'!", read_name);

  num = atoi (ch+1);
  bc |= num << 20;

  for (++i,++ch; i<len; ++i,++ch)
    if (*ch == '_')
      break;
  if (i >= len)
    err_mesg ("invalid read name '%s'!", read_name);

  num = atoi (ch+1);
  bc |= num;

  return bc;
}

/*****************************************************************************/
/*---------------------------- Cluster PEs by BCs ---------------------------*/
/*****************************************************************************/

static int
add1pe (lfr_t * lfr, xh_set_t(bc) * bc_hash, lfr_bc_t * bcs, int64_t n_bc)
{
  xh_item_t * item;
  lfr_t * new_lfr;
  lfr_bc_t * lfr_bc;

  if ((item=xh_set_search3(bc,bc_hash,&lfr->bc)) == NULL)
    err_mesg ("Read %s's barcode is not found in hash!", lfr->r->n->s);

  if (item->id >= n_bc)
    err_mesg ("Barcode index exceeds barcode count!");

  lfr_bc = bcs + item->id;
  lfr_bc->bc = lfr->bc;
  lfr_bc->bc_id = lfr->bc_id;

  new_lfr = mp_alloc (lfr, lfr_bc->reads);
  memcpy (new_lfr, lfr, sizeof(lfr_t));

  return 0;
}

/*****************************************************************************/
/*----------------------------- Thread Routines -----------------------------*/
/*****************************************************************************/

static void *
thread_routines (void * data)
{
  int pid;  
  int n_thread;
  int kmer_len;
  int * signal;
  int64_t i;
  int64_t n_lfr_pes;
  int64_t pth_cur_idx;
  int64_t * cur_idx;
  int64_t * pth_idx_arr;
  mp_t(ctg) * ctgs;
  xh_t ** ctg_khashs;
  lfr_t * lfr;
  pe_seq_t * pe;
  lfr_set_t * lfr_seqs;
  work_para_t * para;
  mp_t(pe) * pes;
  mp_t(lbc) * bcs;
  mp_t(lfr) * lfrs;
  xh_set_t(bc) ** bc_hashs;

  para = (work_para_t *) data;
  pid = para->pid;
  n_thread = para->n_thread;
  kmer_len = para->kmer_len;
  signal = para->signal;
  cur_idx = para->cur_idx;
  ctgs = para->ctgs;
  ctg_khashs = para->ctg_khashs;
  lfr_seqs = para->lfr_seqs;
  bc_hashs = para->bc_hashs;

  pes = lfr_seqs->pes;
  bcs = lfr_seqs->bcs;
  lfrs = lfr_seqs->lfrs;
  pth_idx_arr = lfr_seqs->pth_idx_arr;
  n_lfr_pes = mp_cnt (lfr_seqs->pes);

  for (;;) {
    if (*signal == THREAD_CALC_BC) {
      for (i=0; i<n_lfr_pes; ++i) {
        if (i % n_thread != pid)
          continue;
        pe = mp_at (pe, pes, i);
        lfr = mp_at (lfr, lfrs, i);
        lfr->r = pe;
        lfr->bc = read_name2barcode (pe->n->s, pe->n->l);
        lfr->hs_id = lfr->bc % n_thread;
      }

      *signal = THREAD_WAIT;
    } else if (*signal == THREAD_HASH_BC) {
      for (i=0; i<n_lfr_pes; ++i) {
        lfr = mp_at (lfr, lfrs, i);
        if (lfr->hs_id != pid)
          continue;
        xh_set_add (bc, bc_hashs[pid], &lfr->bc);
      }

      *signal = THREAD_WAIT;
    } else if (*signal == THREAD_CLUSTER) {
      for (i=0; i<n_lfr_pes; ++i) {
        lfr = mp_at (lfr, lfrs, i);
        if (lfr->hs_id != pid)
          continue;
        add1pe (lfr, bc_hashs[pid], mp_at(lbc,bcs,pth_idx_arr[pid]), xh_set_cnt(bc_hashs[pid]));
      }

      *signal = THREAD_WAIT;
    } else if (*signal == THREAD_STOP) {
      *signal = THREAD_WAIT;
      break;
    }
  }

  return (void *) 0;
}

/*****************************************************************************/
/*---------------------------- Extern Functions -----------------------------*/
/*****************************************************************************/

int
ankor_lfr_reads (mp_t(ctg) * ctgs, xh_t ** ctg_khashs,
    lfr_set_t * lfr_seqs, int n_thread, int kmer_len,
    const char * fq1_file, const char * fq2_file)
{
  int i;
  int * signals;
  time_t time_beg;
  time_t mod_tbeg;
  int64_t cur_idx;
  int64_t * n_pth_done;
  pthread_t * pids;
  xh_set_t(bc) ** bc_hashs;
  work_para_t * paras;

  time (&time_beg);

  // chop contigs to kmers
  contig_seqs_info_clear (ctgs);
  chop_contig_seqs2kmers (ctgs, n_thread, kmer_len);

  // store contig kmers to hash
  kmer_hash_clear (ctg_khashs, n_thread);
  put_contig_kmers2hashs (ctg_khashs, ctgs, n_thread);

  // load fastq_files
  lfr_seqs->pes = pefq_load (fq1_file, fq2_file);

  lfr_seqs->n_thread = n_thread;
  lfr_seqs->pth_idx_arr = (int64_t *) ckalloc (n_thread+1, sizeof(int64_t));
  lfr_seqs->bcs = mp_init (lbc, lfr_barcode_init2, NULL);
  lfr_seqs->lfrs = mp_init (lfr, NULL, NULL);
  mp_resize (lfr, lfr_seqs->lfrs, mp_cnt(lfr_seqs->pes));
  lfr_seqs->lfrs->n = mp_cnt (lfr_seqs->pes);

  // create working threads
  pids = (pthread_t *) ckalloc (n_thread, sizeof(pthread_t));
  paras = (work_para_t *) ckalloc (n_thread, sizeof(work_para_t));
  signals = (int *) ckalloc (n_thread, sizeof(int));
  n_pth_done = (int64_t *) ckalloc (n_thread, sizeof(int64_t));

  bc_hashs = (xh_set_t(bc) **) ckalloc (n_thread, sizeof(xh_set_t(bc) *));
  for (i=0; i<n_thread; ++i) {
    bc_hashs[i] = xh_set_init (bc, 65536, 0.75, NULL, NULL, NULL,
        barcode_hash_func, barcode_is_equal);
  }

  for (i=0; i<n_thread; ++i) {
    signals[i] = THREAD_WAIT;
    paras[i].signal = signals + i;
    paras[i].pid = i;

    paras[i].n_thread = n_thread;
    paras[i].kmer_len = kmer_len;

    paras[i].cur_idx = &cur_idx;
    paras[i].ctgs = ctgs;
    paras[i].ctg_khashs = ctg_khashs;
    paras[i].lfr_seqs = lfr_seqs;
    paras[i].bc_hashs = bc_hashs;

    ckpthread_create (pids+i, NULL, thread_routines, (void*)(paras+i));
  }

  // calculate barcode id
  time (&mod_tbeg);
  memset (n_pth_done, 0, n_thread*sizeof(int64_t));
  send_work_signal (signals, THREAD_CALC_BC, n_thread);
  printf ("\n  calculate barcode id cost: %lds\n", time(NULL)-mod_tbeg);

  // hash barcode id
  time (&mod_tbeg);
  memset (n_pth_done, 0, n_thread*sizeof(int64_t));
  send_work_signal (signals, THREAD_HASH_BC, n_thread);
  printf ("\n  hash barcode id cost: %lds\n", time(NULL)-mod_tbeg);

  // cluster pair end reads by barcodes
  time (&mod_tbeg);
  lfr_seqs->pth_idx_arr[0];
  for (i=0; i<n_thread; ++i)
    lfr_seqs->pth_idx_arr[i+1] = lfr_seqs->pth_idx_arr[i] + xh_set_cnt(bc_hashs[i]);
  mp_resize (lbc, lfr_seqs->bcs, lfr_seqs->pth_idx_arr[n_thread]);
  lfr_seqs->bcs->n = lfr_seqs->pth_idx_arr[n_thread];
  send_work_signal (signals, THREAD_CLUSTER, n_thread);
  printf ("\n  cluster reads cost: %lds\n", time(NULL)-mod_tbeg);

  // lfr reads stat
  //lfr_stat (lfr_seqs);

  return 0;
}

void
lfr_free (lfr_set_t * set)
{

}
