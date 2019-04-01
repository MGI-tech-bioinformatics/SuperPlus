/*************************************************
 * File Name: bmp.h
 * Description: Block Memory Pool
 * Author: Chen Xi
 * Mail: chenxi1@genomics.cn
 * Created Time: 2017/10/06
 * Edit History: 
 *************************************************/

#ifndef XDK_BMP_H
#define XDK_BMP_H

#include <stdint.h>

#include "bmp.h"
#include "utils.h"

/*
 * todo:
 * init_type_f==NULL && free_type_f!=NULL
 */

#define BMP_BLK_SIZE 128

#define BMP_DEF2(SCOPE, name, type_t) \
	struct bmp_blk_##name##_s { \
		type_t * datas; \
		int32_t n, m; \
		struct bmp_blk_##name##_s * prev; \
		struct bmp_blk_##name##_s * next; \
	}; \
	typedef struct bmp_blk_##name##_s bmp_blk_##name##_t; \
	struct bmp_##name##_s { \
		bmp_blk_##name##_t * head; \
		bmp_blk_##name##_t * tail; \
		bmp_blk_##name##_t * iter; \
		void (*init_type_f) (type_t*,void*); \
		void * data4init_type; \
		int32_t iter_next_idx; \
		int32_t n_item; \
		int32_t blk_size; \
	}; \
	typedef struct bmp_##name##_s bmp_##name##_t; \
	SCOPE bmp_blk_##name##_t * bmp_blk_init_##name ( \
			void (*init_type_f) (type_t*,void*), \
			int blk_size, \
			void * data4init_type) { \
		int32_t i; \
		bmp_blk_##name##_t * blk; \
		blk = (bmp_blk_##name##_t *) ckalloc (1, sizeof(bmp_blk_##name##_t)); \
		blk->n = 0; \
		blk->m = blk_size>128 ? blk_size : BMP_BLK_SIZE; \
		blk->datas = (type_t *) ckalloc (blk->m, sizeof(type_t)); \
		if (init_type_f != NULL) \
			for (i=0; i<blk->m; i++) \
				init_type_f (blk->datas+i, data4init_type); \
		blk->prev = blk->next = NULL; \
		return blk; \
	} \
	SCOPE void bmp_blk_free_##name ( \
			bmp_blk_##name##_t * blk, \
			void (*free_type_f) (type_t*)) { \
		int32_t i; \
		if (free_type_f != NULL) \
			for (i=0; i<blk->m; i++) \
				free_type_f (blk->datas + i); \
		free (blk->datas); \
		free (blk); \
	} \
	SCOPE void bmp_blk_clear_##name ( \
			bmp_blk_##name##_t * blk, \
			void (*clear_type_f) (type_t *)) { \
		int32_t i; \
		if (clear_type_f != NULL) \
			for (i=0; i<blk->n; i++) \
				clear_type_f (blk->datas + i); \
		blk->n = 0; \
	} \
	SCOPE bmp_##name##_t * bmp_init_##name ( \
			int blk_size, \
			void (*init_type_f) (type_t*,void*), \
			void * data4init_type) { \
		bmp_##name##_t * bmp; \
		bmp = (bmp_##name##_t *) ckalloc (1, sizeof(bmp_##name##_t)); \
		bmp->n_item = 0; \
		bmp->blk_size = blk_size; \
		bmp->init_type_f = init_type_f; \
		bmp->data4init_type = data4init_type; \
		bmp->head = bmp->tail = bmp_blk_init_##name (init_type_f, blk_size, data4init_type); \
		return bmp; \
	} \
	SCOPE void bmp_free_##name ( \
			bmp_##name##_t * bmp, \
			void (*free_type_f) (type_t*)) { \
		bmp_blk_##name##_t * blk; \
		bmp_blk_##name##_t * next_blk; \
    if (NULL!=free_type_f && NULL==bmp->init_type_f) \
      warn_mesg ("init_type_f == NULL, while free_type_f != NULL!"); \
		blk = bmp->head; \
		while (blk != NULL) { \
      next_blk = blk->next; \
			bmp_blk_free_##name (blk, free_type_f); \
			blk = next_blk; \
		} \
		free (bmp); \
	} \
	SCOPE void bmp_clear_##name ( \
			bmp_##name##_t * bmp, \
			void (*clear_type_f) (type_t *)) { \
		bmp_blk_##name##_t * blk; \
		blk = bmp->head; \
		while (blk != NULL) { \
			bmp_blk_clear_##name (blk, clear_type_f); \
			blk = blk->next; \
		} \
		bmp->n_item = 0; \
		bmp->tail = bmp->head; \
	} \
	SCOPE type_t * bmp_alloc_##name ( \
			bmp_##name##_t * bmp) { \
		bmp_blk_##name##_t * blk; \
		blk = bmp->tail; \
		++bmp->n_item; \
		if (blk->n < blk->m) { \
			return blk->datas + blk->n++; \
		} \
		if (blk->next != NULL) { \
			blk = blk->next; \
			bmp->tail = blk; \
			return blk->datas + blk->n++; \
		} \
		blk = bmp_blk_init_##name (bmp->init_type_f, bmp->blk_size, bmp->data4init_type); \
		bmp->tail->next = blk; \
    blk->prev = bmp->tail; \
		bmp->tail = blk; \
		return blk->datas + blk->n++; \
	} \
  SCOPE type_t * bmp_pop_##name ( \
      bmp_##name##_t * bmp) { \
    bmp_blk_##name##_t * blk; \
    blk = bmp->tail; \
    --blk->n; \
    if (blk->n == 0) \
      bmp->tail = bmp->tail->prev; \
    return blk->datas; \
  } \
	SCOPE int bmp_iter_init_##name ( \
			bmp_##name##_t * bmp) { \
		if (bmp->head == NULL) \
			return -1; \
		bmp->iter = bmp->head; \
		bmp->iter_next_idx = 0; \
		return 0; \
	} \
	SCOPE type_t * bmp_iter_next_##name ( \
			bmp_##name##_t * bmp) { \
		if (bmp->iter_next_idx < bmp->iter->n) \
			return bmp->iter->datas + bmp->iter_next_idx++; \
    if (bmp->iter == bmp->tail) \
      return NULL; \
		bmp->iter = bmp->iter->next; \
		bmp->iter_next_idx = 0; \
		if (bmp->iter == NULL) \
			return NULL; \
		return bmp->iter->datas + bmp->iter_next_idx++; \
	} \
	SCOPE type_t * bmp_at_##name ( \
			bmp_##name##_t * bmp, \
			int32_t idx) { \
		bmp_blk_##name##_t * blk; \
		if (idx >= bmp->n_item) \
			return NULL; \
		blk = bmp->head; \
		while (idx >= bmp->blk_size) { \
			idx -= bmp->blk_size; \
			if (blk == NULL) \
				err_mesg ("[%s] fail to locate item."); \
			blk = blk->next; \
		} \
		return blk->datas + idx; \
	} \
	static int __bmp_##name##_def_end = 1

#define BMP_DEF(name, type) \
	BMP_DEF2(static inline, name, type)

#define bmp_t(name) bmp_##name##_t
#define bmp_cnt(bmp) ((bmp)->n_item)

#define bmp_init(name,blk_size,func,data) bmp_init_##name((blk_size),(func),(data))
#define bmp_free(name,bmp,func) bmp_free_##name(bmp,func)
#define bmp_clear(name,bmp,func) bmp_clear_##name(bmp,func)
#define bmp_alloc(name,bmp) bmp_alloc_##name(bmp)
#define bmp_pop(name,bmp) bmp_pop_##name(bmp)
#define bmp_iter_init(name,bmp) bmp_iter_init_##name(bmp)
#define bmp_iter_next(name,bmp) bmp_iter_next_##name(bmp)
#define bmp_at(name,bmp,idx) bmp_at_##name(bmp,idx)

#endif
