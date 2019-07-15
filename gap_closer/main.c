/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-11-30 17:41:31
  *Edit History: 
***********************************************************/

#include <time.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "mp.h"
#include "def.h"
#include "lfr.h"
#include "str.h"
#include "hash.h"
#include "kseq1.h"
#include "utils.h"
#include "contig.h"
#include "rseq.h"
#include "kmer.h"
#include "ont.h"
#include "gap_closer.h"

static int n_kmers = 1;
static int kmer_lens[] = {25, 21, 19, 17};

typedef struct kmer_map_s {
  int kmer_len;
  xh_set_t(kmer) ** anchored_ksets;
} kmer_map_t;

typedef struct {
  int n_kmers;
  int * kmer_lens;
  kmer_map_t * kmap_infos;
  mp_t(ctg) * ctgs;
  mp_t(rs) * ont_seqs;
  xh_t ** ctg_khashs;
  mp_t(okseq) * okseqs;
  ctg_graph_t * ctg_graph;
  mp_t(sf) * scafs;
  lfr_set_t * lfr_seqs;
} gc_info_t;

static gc_info_t *
gap_closer_info_init (int n_thread)
{
  int i, j;
  gc_info_t * gc_info;
  kmer_map_t * kmap_info;

  gc_info = (gc_info_t *) ckmalloc (sizeof(gc_info_t));
  gc_info->n_kmers = n_kmers;

  gc_info->kmer_lens = (int *) ckmalloc (n_kmers * sizeof(int));
  for (i=0; i<n_kmers; ++i)
    gc_info->kmer_lens[i] = kmer_lens[i];

  gc_info->kmap_infos = (kmer_map_t *) ckmalloc (n_kmers * sizeof(kmer_map_t));
  for (i=0; i<n_kmers; ++i) {
    kmap_info = gc_info->kmap_infos + i;
    kmap_info->kmer_len = kmer_lens[i];
    kmap_info->anchored_ksets = (xh_set_t(kmer) **) ckmalloc (n_thread * sizeof(xh_set_t(kmer) *));
    for (j=0; j<n_thread; ++j) {
      kmap_info->anchored_ksets[j] = xh_set_init (kmer, 65536, 0.75, NULL, NULL,
          kmer_copy_func, kmer_hash_func, kmer_is_equal);
    }
  }

  gc_info->ctg_khashs = kmer_hash_init (n_thread);
  gc_info->ctg_graph = ctg_graph_init ();
  gc_info->okseqs = mp_init (okseq, okseq_init2, NULL);

  gc_info->lfr_seqs = (lfr_set_t *) ckalloc (1, sizeof(lfr_set_t));
  gc_info->scafs = mp_init (sf, scaf_init2, NULL);

  return gc_info;
}

static void
gap_closer_info_free (gc_info_t * gc_info, int n_thread)
{
  int i, j;
  kmer_map_t * kmap_info;

  for (i=0; i<n_kmers; ++i) {
    kmap_info = gc_info->kmap_infos + i;
    for (j=0; j<n_thread; ++j)
      xh_set_free (kmer, kmap_info->anchored_ksets[j], kmer_free2);
    free (kmap_info->anchored_ksets);
  }
  free (gc_info->kmap_infos);

  free (gc_info->kmer_lens);

  contig_seqs_free (gc_info->ctgs);
  kmer_hash_free (gc_info->ctg_khashs, n_thread);

  ctg_graph_free (gc_info->ctg_graph);

  mp_free (sf, gc_info->scafs, scaf_free2);

  mp_free (okseq, gc_info->okseqs, okseq_free2);
  mp_free (rs, gc_info->ont_seqs, rseq_free2);
  lfr_free (gc_info->lfr_seqs);

  free (gc_info);
}

/*****************************************************************************/
/*------------------------- scaff kmer gap closer ---------------------------*/
/*****************************************************************************/

int
main (int argc, char * argv[])
{
  if (argc < 5) {
    fprintf (stderr, "Incomplete input!\n");
    fprintf (stderr, "usage: gc <scaff> <ont> <number of thread> <output prefix>\n");
    return 1;
  }

  char * scaff_file;
  char * ont_file;
  char * lfr_fq_file[2];
  char * prefix;
  int i;
	int n_thread;
  time_t time_beg;
  gc_info_t * gc_info;
  kmer_map_t * kmap_info;

  time (&time_beg);

  scaff_file = argv[1];
  ont_file = argv[2];
  n_thread = atoi (argv[3]);
  //lfr_fq_file[0] = argv[4];
  //lfr_fq_file[1] = argv[5];
  prefix = argv[4];

  hash_func_init ();

  // Init gc info
  gc_info = gap_closer_info_init (n_thread);

  // load contig seq
	gc_info->ctgs = contig_seqs_load (scaff_file, gc_info->ctg_graph, gc_info->scafs);

  // load ont reads
  gc_info->ont_seqs = sefq_load (ont_file);
  ont_kseqs_init (gc_info->ont_seqs, gc_info->okseqs);

  for (i=0; i<gc_info->n_kmers; ++i) {
	  printf ("\n============================\n");
	  printf ("====== Kmer Length %d ======\n", gc_info->kmer_lens[i]);
	  printf ("============================\n\n");
	  kmap_info = gc_info->kmap_infos + i;
	
	  // chop contigs to kmers
	  printf ("  >>> Chop Contig Seqs <<<\n");
	  contig_seqs_info_clear (gc_info->ctgs);
	  chop_contig_seqs2kmers (gc_info->ctgs, n_thread, gc_info->kmer_lens[i]);
	  printf ("\n");
	
	  // store contig kmers to hash
	  printf ("  >>> Hash Scaf Kmers <<<\n");
    kmer_hash_clear (gc_info->ctg_khashs, n_thread);
	  put_contig_kmers2hashs (gc_info->ctg_khashs, gc_info->ctgs, n_thread);
	  printf ("\n");
	
	  // search kmers in ont reads
	  printf ("  >>> Search ONT Kmers <<<\n");
	  search_kmers_on_ont_reads (gc_info->ont_seqs, gc_info->ctgs,
	      gc_info->ctg_khashs, gc_info->okseqs, kmap_info->anchored_ksets,
	      prefix, n_thread, gc_info->kmer_lens[i]);
	  printf ("\n");
	
	  // kmer stat
	  printf ("  >>> Kmer Stat <<<\n");
	  kmer_stat (gc_info->ctg_khashs, n_thread, gc_info->kmer_lens[i], "Scaffold");
	  kmer_stat2 (kmap_info->anchored_ksets, n_thread, gc_info->kmer_lens[i], "ONT");
	  printf ("\n");
  }

  // ankor lfr reads to contigs
  //ankor_lfr_reads (gc_info->ctgs, gc_info->ctg_khashs, gc_info->lfr_seqs,
  //    n_thread, 25, lfr_fq_file[0], lfr_fq_file[1]);

  // gap closing
	printf ("\n============================\n");
  printf ("======= Gap Closing ========\n");
	printf ("============================\n\n");
  gap_closing (gc_info->ctgs, gc_info->ont_seqs, gc_info->lfr_seqs,
      gc_info->okseqs, gc_info->scafs, gc_info->ctg_graph,
      kmap_info->anchored_ksets, n_thread);
  printf ("\n");

  hash_func_free ();

  // free
  gap_closer_info_free (gc_info, n_thread);

  printf ("\nProgram Cost: %lds\n", time(NULL)-time_beg);

  return 0;
}
