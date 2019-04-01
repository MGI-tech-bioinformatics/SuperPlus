/*************************************************
 * File Name: ../../ds/array.h
 * Description: 
 * Author: Chen Xi
 * Mail: chenxi1@genomics.cn
 * Created Time: 2018/11/30
 * Edit History: 
 *************************************************/

#ifndef XDK_ARRAY_H
#define XDK_ARRAY_H

#include <stdint.h>

#include "utils.h"

#define ARR_INIT_SIZE 16

#define ARR_DEF(name, type_t) \
	struct arr_##name##_s { \
		int32_t n, m; \
		type_t * arr; \
	}; \
	typedef struct arr_##name##_s arr_##name##_t; \
	\
	static inline arr_##name##_t * arr_init_##name (void) { \
		arr_##name##_t * arr; \
		arr = (arr_##name##_t *) ckalloc (1, sizeof(arr_##name##_t)); \
		arr->n = 0; \
		arr->m = ARR_INIT_SIZE; \
		arr->arr = (type_t *) ckalloc (ARR_INIT_SIZE, sizeof(type_t)); \
		return arr; \
	} \
	\
	static inline void arr_free_##name ( \
			arr_##name##_t * arr) { \
		free (arr->arr); \
		free (arr); \
	} \
	\
	static inline void arr_clear_##name ( \
			arr_##name##_t * arr) { \
		arr->n = 0; \
		memset (arr->arr, 0, arr->m*sizeof(type_t)); \
	} \
	\
	static inline void arr_resize_##name ( \
			arr_##name##_t * arr, \
			int32_t new_size) { \
		int32_t old_max; \
		if (new_size <= arr->m) \
			return; \
		old_max = arr->m; \
		while (new_size > arr->m) { \
			if (arr->m < 0x10000) \
				arr->m <<= 1; \
			else \
				arr->m += 0x10000; \
		} \
		arr->arr = (type_t *) ckrealloc (arr->arr, arr->m*sizeof(type_t)); \
		memset (arr->arr+old_max, 0, (arr->m-old_max)*sizeof(type_t)); \
	} \
	\
	static inline type_t * arr_at_##name ( \
			arr_##name##_t * arr, \
			int32_t idx) { \
		arr_resize_##name (arr, idx+1); \
		return arr->arr + idx; \
	} \
  \
  static inline type_t * arr_alloc_##name ( \
      arr_##name##_t * arr) { \
    arr_resize_##name (arr, arr->n+1); \
    return arr->arr + arr->n++; \
  } \
  \
  static inline void arr_copy_##name ( \
      arr_##name##_t * dst, \
      arr_##name##_t * src) { \
    arr_resize_##name (dst, src->n); \
    dst->n = src->n; \
    memcpy (dst->arr, src->arr, src->n*sizeof(type_t)); \
  } \
  \
	static int __array_##name##_def_end__ = 1

#define arr_t(name) arr_##name##_t
#define arr_cnt(arr) (arr->n)

#define arr_init(name) arr_init_##name()
#define arr_free(name,arr) arr_free_##name(arr)
#define arr_clear(name,arr) arr_clear_##name(arr)
#define arr_resize(name,arr,new_size) arr_resize_##name(arr,new_size)
#define arr_at(name,arr,idx) arr_at_##name(arr,idx)
#define arr_alloc(name,arr) arr_alloc_##name(arr)
#define arr_copy(name,dst,src) arr_copy_##name((dst),(src))

#endif
