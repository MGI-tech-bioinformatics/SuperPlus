/*************************************************
 * File Name: contig.h
 * Description: 
 * Author: Chen Xi
 * Mail: chenxi1@genomics.cn
 * Created Time: 2018/12/09
 * Edit History: 
 *************************************************/

#ifndef SEQ_H
#define SEQ_H

#include <stdint.h>

#include "mp.h"
#include "def.h"
#include "ctg_graph.h"

mp_t(ctg) * contig_seqs_load (const char * scaffoldseq_file, ctg_graph_t * ctg_graph, mp_t(sf) * scafs);

void contig_seqs_info_clear (mp_t(ctg) * seqs);

void contig_seqs_free (mp_t(ctg) * seqs);

#endif
