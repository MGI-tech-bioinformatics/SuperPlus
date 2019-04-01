/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-05-16 15:07:45
  *Edit History: 
***********************************************************/

#ifndef XDK_MP_H
#define XDK_MP_H

#include <stdint.h>
#include <stdlib.h>

#include "utils.h"

#define MP_INIT_SIZE 2

#define MP_DEF2(SCOPE, name, type_t) \
	struct mp_##name##_s { \
		int64_t n, m; \
    int64_t mn; \
		type_t * pool; \
		void (*init_type_f) (type_t*,void*); \
		void * data4init_type; \
	}; \
	typedef struct mp_##name##_s mp_##name##_t; \
	\
	SCOPE mp_##name##_t * mp_init_##name ( \
			void (*init_type_f)(type_t*,void*), \
			void * data4init_type) { \
		int64_t i; \
		mp_##name##_t * mp; \
		mp = (mp_##name##_t *) ckalloc (1, sizeof(mp_##name##_t)); \
		mp->n = 0; \
    mp->mn = 0; \
		mp->m = MP_INIT_SIZE; \
		mp->init_type_f = init_type_f; \
		mp->data4init_type = data4init_type; \
		mp->pool = (type_t *) ckalloc (mp->m, sizeof(type_t)); \
		if (init_type_f != NULL) \
			for (i=0; i<mp->m; i++) \
				init_type_f (mp->pool+i, data4init_type); \
		return mp; \
	} \
	\
	SCOPE void mp_free_##name (mp_##name##_t * mp, void (*free_type_f)(type_t*)) { \
		int64_t i; \
    int64_t max; \
    if (mp->init_type_f == NULL) \
      max = mp->n>mp->mn ? mp->n : mp->mn; \
    else \
      max = mp->m; \
		if (free_type_f != NULL) { \
			for (i=0; i<max; i++) \
				free_type_f (mp->pool + i); \
		} \
		free (mp->pool); \
		free (mp); \
	} \
	SCOPE void mp_clear_##name (mp_##name##_t * mp, void (*clear_type_f)(type_t*)) { \
		int64_t i; \
		if (clear_type_f != NULL) { \
			for (i=0; i<mp->n; i++) \
				clear_type_f (mp->pool + i); \
		} \
    if (mp->n > mp->mn) \
      mp->mn = mp->n; \
		mp->n = 0; \
	} \
	SCOPE void mp_resize_##name (mp_##name##_t * mp, int64_t new_size) { \
		int64_t i, old_max; \
		if (mp->m >= new_size) \
			return; \
		old_max = mp->m; \
		if (mp->m <= 0) \
			mp->m = MP_INIT_SIZE; \
		while (mp->m < new_size) { \
			if (mp->m < 0x10000) \
				mp->m <<= 1; \
			else \
				mp->m += 0x10000; \
		} \
		mp->pool = (type_t *) ckrealloc (mp->pool, mp->m*sizeof(type_t)); \
		if (mp->init_type_f == NULL) \
			memset (mp->pool+old_max, 0, (mp->m-old_max)*sizeof(type_t)); \
		else { \
			for (i=old_max; i<mp->m; i++) \
				mp->init_type_f (mp->pool+i, mp->data4init_type); \
		} \
	} \
	SCOPE type_t * mp_alloc_##name (mp_##name##_t * mp) { \
		mp_resize_##name (mp, mp->n+1); \
		return mp->pool + mp->n++; \
	} \
	SCOPE int64_t mp_add_##name (mp_##name##_t * mp, type_t * src, \
			void (*copy_type_f)(type_t*,type_t*)) { \
		type_t * dst; \
		dst = mp_alloc (name, mp); \
		if (src == NULL) \
			err_mesg ("src == NULL"); \
		if (copy_type_f == NULL) \
			memcpy (dst, src, sizeof(type_t)); \
		else \
			copy_type_f (dst, src); \
		return mp->n - 1; \
	} \
	SCOPE type_t * mp_at_##name (mp_##name##_t * mp, int64_t idx) { \
		if (idx >= mp->n) \
			return NULL; \
		return mp->pool + idx; \
	} \
	SCOPE void mp_dump_##name (mp_##name##_t * mp, const char * file, \
			void(*type2str_f)(type_t*,char*,int64_t)) { \
		char * buf; \
		int64_t i; \
		FILE * fp; \
		buf = ALLOC_LINE; \
		fp = ckopen (file, "w"); \
		for (i=0; i<mp->n; i++) { \
			type2str_f (mp->pool+i, buf, LINE_MAX); \
			fprintf (fp, "%s", buf); \
		} \
		fclose (fp); \
		free (buf); \
	} \
	SCOPE void mp_copy_##name (mp_##name##_t * dst, mp_##name##_t * src, \
			void (*copy_type_f)(type_t*,type_t*)) { \
		int64_t i; \
		mp_resize_##name (dst, src->n); \
		dst->n = src->n; \
		if (copy_type_f == NULL) { \
			memcpy (dst->pool, src->pool, src->n*sizeof(type_t)); \
		} else { \
			for (i=0; i<src->n; ++i) \
				copy_type_f (dst->pool+i, src->pool+i); \
		} \
	} \
  SCOPE void mp_append_##name (mp_##name##_t * dst, mp_##name##_t * src, \
      void (*copy_type_f)(type_t*,type_t*)) { \
    int64_t i; \
    int64_t n_old; \
    type_t * dst_ptr; \
    type_t * src_ptr; \
    n_old = dst->n; \
    dst->n += src->n; \
    mp_resize_##name (dst, dst->n); \
    dst_ptr = dst->pool + n_old; \
    src_ptr = src->pool; \
    if (copy_type_f == NULL) { \
      memcpy (dst_ptr, src_ptr, src->n*sizeof(type_t)); \
    } else { \
      for (i=0; i<src->n; ++i,++dst_ptr,++src_ptr) \
        copy_type_f (dst_ptr, src_ptr); \
    } \
  } \
	static int __mp_##name##_def_end__ = 1

#define MP_DEF(name, type_t) \
	MP_DEF2(static inline, name, type_t)

#define mp_t(name) mp_##name##_t
#define mp_cnt(mp) ((mp)->n)

#define mp_init(name,func,data) mp_init_##name((func), (data))
#define mp_free(name,mp,f) mp_free_##name((mp), (f))
#define mp_clear(name,mp,f) mp_clear_##name((mp), (f))
#define mp_resize(name,mp,size) mp_resize_##name((mp), (size))
#define mp_alloc(name,mp) mp_alloc_##name(mp)
#define mp_add(name,mp,ptr,copy_f) mp_add_##name((mp), (ptr), (copy_f))
#define mp_at(name,mp,idx) mp_at_##name((mp), (idx))
#define mp_dump(name,mp,file,f) mp_dump_##name((mp), (file), (f))
#define mp_copy(name,dst,src,f) mp_copy_##name((dst), (src), (f))
#define mp_append(name,dst,src,f) mp_append_##name((dst), (src), (f))

#endif
