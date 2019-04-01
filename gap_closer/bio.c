/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-05-22 11:38:27
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "str.h"
#include "bmp.h"
#include "utils.h"
#include "bio.h"
#include "hash.h"

const char base_rc_tbl[128] = {
  'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
  'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
  'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
  'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
  'N','T','N','G', 'N','N','N','C', 'N','N','N','N', 'N','N','N','N',
  'N','N','N','N', 'A','N','N','N', 'N','N','N','N', 'N','N','N','N',
  'N','T','N','G', 'N','N','N','C', 'N','N','N','N', 'N','N','N','N',
  'N','N','N','N', 'A','N','N','N', 'N','N','N','N', 'N','N','N','N',
  };

const char base2int_tbl[128] = {
	4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
	4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
	4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
	4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
	4,0,4,1, 4,4,4,3, 4,4,4,4, 4,4,4,4,
	4,4,4,4, 2,4,4,4, 4,4,4,4, 4,4,4,4,
	4,0,4,1, 4,4,4,3, 4,4,4,4, 4,4,4,4,
	4,4,4,4, 2,4,4,4, 4,4,4,4, 4,4,4,4,
	};

HASH_MAP_DEF (gchr, str_t, int32_t);

/**********************************************************
 ***************** Seq Address Functions ******************
 **********************************************************/

str_t *
reverse_complement1seq_i (str_t * seq, str_t * rv_seq)
{
	int32_t i, j;

	if (seq->l <= 0)
		err_mesg ("[%s] seq is empty!");

	if (rv_seq == NULL)
		rv_seq = str_init ();
	str_resize (rv_seq, seq->l);

	for (i=seq->l-1,j=0; i>=0; --i,++j)
		rv_seq->s[j] = int_comp (seq->s[i]);

	rv_seq->l = seq->l;

	return 0;
}

str_t *
reverse_complement1seq_b (str_t * seq, str_t * rv_seq)
{
	int32_t i, j;

	if (seq->l <= 0)
		err_mesg ("[%s] seq is empty!");

	if (rv_seq == NULL)
		rv_seq = str_init ();
	str_resize (rv_seq, seq->l);

	for (i=seq->l-1,j=0; i>=0; --i,++j)
		rv_seq->s[j] = int2base (int_comp(base2int(seq->s[i])));

	rv_seq->l = seq->l;

	return 0;
}

/**********************************************************
 ******************** Static Functions ********************
 **********************************************************/

static int
validate_ref_idx (const char * ref_idx)
{
	return 0;
}

/**********************************************************
 **************** Genome Address Functions ****************
 **********************************************************/

genome_info_t *
load_ref_genome_info (const char * ref)
{
	char * chr_name;
	char * line;
	char * fai;
	int32_t i;
	FILE * fp;
	genome_info_t * ginfo;

	chr_name = ALLOC_LINE;
	line     = ALLOC_PATH;
	fai      = ALLOC_PATH;

	if (check_ref_idx(ref, NULL) != 0)
		err_mesg ("[%s] '%s.fai' do not exist and fail to create it from '%s' in file %s at line %d!",
				__func__, ref, ref, __LINE__, __LINE__);

	ginfo = (genome_info_t *) ckalloc (1, sizeof(genome_info_t));

	sprintf (fai, "%s.fai", ref);
	fp = ckopen (fai, "r");
	ginfo->n_targets = 0;
	while (fgets(line, LINE_MAX, fp))
		++ginfo->n_targets;

	ginfo->chr_idx_array = (uint32_t *) ckalloc (ginfo->n_targets+1, sizeof(uint32_t));
	ginfo->target_len = (uint32_t *) ckalloc (ginfo->n_targets, sizeof(uint32_t));
	ginfo->target_name = (str_t *) ckalloc (ginfo->n_targets, sizeof(str_t));

	rewind(fp), i=0;
	ginfo->max_chr_len = 0;
	ginfo->chr_idx_array[0] = 0;
	while (fgets(line, LINE_MAX, fp)) {
		sscanf (line, "%s %u", chr_name, ginfo->target_len+i);

		str_init2 (ginfo->target_name+i, NULL);
		str_assign (ginfo->target_name+i, chr_name);

		if (ginfo->target_len[i] > ginfo->max_chr_len)
			ginfo->max_chr_len = ginfo->target_len[i];

		ginfo->chr_idx_array[i+1] = ginfo->chr_idx_array[i] + ginfo->target_len[i];

		++i;
	}

	ginfo->genome_len = ginfo->chr_idx_array[i];
	ginfo->gseq = NULL;
	ginfo->chr_hash = NULL;

	fclose (fp);

	free (chr_name);
	free (line);
	free (fai);

	return ginfo;
}

void
free_genome_info (genome_info_t * ginfo)
{
	int32_t i;

	for (i=0; i<ginfo->n_targets; ++i)
		str_free2 (ginfo->target_name+i);
	free (ginfo->target_name);

	free (ginfo->chr_idx_array);
	free (ginfo->target_len);

	if (ginfo->gseq != NULL)
		free (ginfo->gseq);

	if (ginfo->chr_hash != NULL)
		xh_map_free (gchr, ginfo->chr_hash, str_free2, NULL);

	free (ginfo);
}

int
check_ref_idx (const char * ref, const char * samtools)
{
	char * cmd;
	char * ref_idx;

	cmd     = ALLOC_LINE;
	ref_idx = ALLOC_PATH;

	sprintf (ref_idx, "%s.fai", ref);
	if (validate_ref_idx (ref_idx) != 0) {
		warn_mesg ("[%s] No reference index file exist in file at line %d!",
				__func__, __FILE__, __LINE__);

		if (samtools == NULL)
			err_mesg ("[%s] samtools==NULL, fail to  create reference index file in file %s at line %d!",
					__func__, __FILE__, __LINE__);

		warn_mesg ("[%s] Reference index file will be created for you ...", __func__);
		sprintf (cmd, "%s faidx %s", samtools, ref);
		system (cmd);

		if (validate_ref_idx (ref_idx) != 0) {
			free (cmd);
			free (ref_idx);
			err_mesg ("[%s] Fail to create reference index file in file %s at line %d!",
					__func__, __FILE__, __LINE__);
		} else
			warn_mesg ("[%s] Create reference index in file %s line %d!",
					__func__, __FILE__, __LINE__);
	}

	free (cmd);
	free (ref_idx);

	return 0;
}

char * load_genome_seq
(const char * ref, const genome_info_t * ginfo)
{
	char * ch;
	char * ch2;
	char * gseq;
	char * line;
	FILE * fp;

	line = ALLOC_LINE;
	fp = ckopen (ref, "r");
	gseq = (char *) ckalloc (ginfo->genome_len+1, sizeof(char));
	ch = gseq;
	while (fgets(line, LINE_MAX, fp)) {
		if (*line == '>')
			continue;
		for (ch2=line; *ch2!='\0'; ++ch2,++ch)
			*ch = *ch2;
		if (*(ch-1) == '\n')
			--ch;
	}

	fclose (fp);
	free (line);

	return gseq;
}

char *
get_ref_seq (genome_info_t * ginfo, int32_t tid, uint32_t beg)
{
	if (ginfo->gseq == NULL)
		err_mesg ("genome sequence are not available!");

	return ginfo->gseq + ginfo->chr_idx_array[tid] + beg;
}

int
get_ref_seq2 (genome_info_t * ginfo, int32_t tid, uint32_t beg, uint32_t len, str_t * s)
{
	char * ch;

	ch = get_ref_seq (ginfo, tid, beg);

	if (s == NULL)
		s = str_init ();

	str_resize (s, len+1);
	s->l = len;
	memcpy (s->s, ch, len);
	s->s[s->l] = '\0';

	return 0;
}

int
get_ref_seq3 (genome_info_t * ginfo, uint32_t gbeg, uint32_t len, str_t * s)
{
	char * ch;

	ch = ginfo->gseq + gbeg;

	if (s == NULL)
		s = str_init ();

	str_resize (s, len+1);
	s->l = len;
	memcpy (s->s, ch, len);
	s->s[s->l] = '\0';

	return 0;
}

xh_map_t(gchr) *
hash_ref_chr (const char * ref, const char * samtools)
{
	char * line;
	char * chr_name;
	char * ref_idx;
	int32_t i;
	FILE * fp;
	str_t * s;
	xh_map_t(gchr) * set;

	s = str_init ();
	set = xh_map_init (gchr, 256, 0.75, str_init2, NULL, NULL, NULL, str_copy, NULL, gchr_key_hash_func, gchr_key_equal_func);

	ref_idx  = ALLOC_PATH;
	line     = ALLOC_LINE;
	chr_name = ALLOC_LINE;

	sprintf (ref_idx, "%s.fai", ref);
	if (validate_ref_idx(ref_idx) != 0)
		err_mesg ("[%s] '%s.fai' do not exist or is invalid in file %s at line %d!",
				__func__, ref, __FILE__, __LINE__);

	i = 0;
	fp = ckopen (ref_idx, "r");
	while (fgets(line, LINE_MAX, fp)) {
		sscanf (line, "%s", chr_name);
		str_assign (s, chr_name);
		xh_map_add (gchr, set, s, &i);
		++i;
	}

	free (ref_idx);
	free (line);
	free (chr_name);

	fclose (fp);

	str_free (s);

	return set;
}

xh_map_t(gchr) *
hash_info_chr (const genome_info_t * ginfo)
{
	int32_t i;
	xh_map_t(gchr) * set;

	set = xh_map_init (gchr, 256, 0.75, str_init2, NULL, NULL, NULL, str_copy, NULL, gchr_key_hash_func, gchr_key_equal_func);

	for (i=0; i<ginfo->n_targets; i++)
		xh_map_add (gchr, set, ginfo->target_name+i, &i);

	return set;
}

int32_t
get_chr_idx (xh_map_t(gchr) * chr_hash, str_t * chr_name)
{
	int32_t * val_ptr;

	if (chr_hash == NULL)
		return -1;

	val_ptr = xh_map_search (gchr, chr_hash, chr_name);

	if (val_ptr == NULL)
		return -1;

	return *val_ptr;
}
