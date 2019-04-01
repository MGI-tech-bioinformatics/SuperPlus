/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-05-13 19:07:11
  *Edit History: 
***********************************************************/

#ifndef XDK_UTILS_H
#define XDK_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <stdint.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/types.h>

#define XUTILS_STATUS "utested"

/**********************************************************
 ******************* Macro Definitions ********************
 **********************************************************/

#define IS_N 1
#define NOT_N 0

#ifndef LINE_MAX
#define LINE_MAX 4096
#endif

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#ifndef LONG_LINE_MAX
#define LONG_LINE_MAX 16777216
#endif

#define ckopen(f,m) (err_ckopen((f),(m), \
			__func__, __FILE__, __LINE__))
#define ckpopen(f,m) (err_ckpopen((f),(m), \
			__func__, __FILE__, __LINE__))
#define xopen(f,m) (err_xopen((f),(m), \
			__func__, __FILE__, __LINE__))
#define ckfwrite(ptr,size,n,fp) (err_ckfwrite((ptr),(size),(n),(fp), \
			__func__, __FILE__, __LINE__))
#define ckfflush(fp) (err_ckfflush((fp), \
			__func__, __FILE__, __LINE__))

#define ckopendir(path) (err_ckopendir((path), \
			__func__, __FILE__, __LINE__))
#define ckcreate_dir(path) (err_ckcreate_dir((path), \
			__func__, __FILE__, __LINE__))

#define rm1file(path) (err_rm1file((path), \
			__func__, __FILE__, __LINE__))
#define rm1dir(path) (err_rm1dir((path), \
			__func__, __FILE__, __LINE__))

#define ckalloc(n,size) (err_ckalloc((n),(size), \
			__func__, __FILE__, __LINE__))
#define ckmalloc(size) (err_ckmalloc((size), \
			__func__, __FILE__, __LINE__))
#define ckrealloc(old,new_size) (err_ckrealloc((old),(new_size), \
			__func__, __FILE__, __LINE__))
//#define ckaligned_malloc(align,size) (err_ckaligned_malloc((align),(size), \
			__func__, __FILE__, __LINE__))

#define str_cmp_tail(s1,s2,n) (err_str_cmp_tail((s1),(s2),(n), \
			__func__, __FILE__, __LINE__))
#define str_cut_tail(s,t) (err_str_cut_tail((s),(t), \
			__func__, __FILE__, __LINE__))
#define chomp(s) (err_chomp((s), \
			__func__, __FILE__, __LINE__))
#define chop(s,c) (err_chop((s),(c), \
			__func__, __FILE__, __LINE__))
#define split_string(s,m,c) (err_split_string((s),(m),(c), \
			__func__, __FILE__, __LINE__))
#define is_empty_line(s) (err_is_empty_line((s), \
			__func__, __FILE__, __LINE__))

#define get_abs_path(p) (err_get_abs_path((p), \
			__func__, __FILE__, __LINE__))
#define get_program_path(p) (err_get_program_path((p), \
			__func__, __FILE__, __LINE__))

#define ALLOC_LINE ((char*)ckalloc(LINE_MAX,sizeof(char)))
#define ALLOC_PATH ((char*)ckalloc(PATH_MAX,sizeof(char)))
#define ALLOC_LONG_LINE ((char*)ckalloc(LONG_LINE_MAX,sizeof(char)))

#define xswap(a,b,t) ((t)=(a),(a)=(b),(b)=(t))
#define xmax(a,b) ((a)>(b)?(a):(b))
#define xmin(a,b) ((a)<(b)?(a):(b))

#define XDK_TRUE 1
#define XDK_FALSE 0

/**********************************************************
 *************** IO, FileSystem Functions *****************
 **********************************************************/

#ifdef __cplusplus
extern "C" {
#endif

	FILE * err_ckopen (const char * file, const char * mode,
			const char * c_func, const char * c_file, int line_c);

	FILE * err_ckpopen (const char * cmd, const char * mode,
			const char * c_func, const char * c_file, int line_c);

	int err_xopen (const char * file, int flags,
			const char * c_func, const char * c_file, int line_c);

	size_t err_ckfwrite (const void * ptr, size_t size, size_t n_elem, FILE * fp,
			const char * c_func, const char * c_file, int line_c);

	int err_ckfflush (FILE * fp,
			const char * c_func, const char * c_file, int line_c);

	DIR * err_ckopendir (const char * path,
			const char * c_func, const char * c_file, int line_c);

	int err_ckcreate_dir (const char * path,
			const char * c_func, const char * c_file, int line_c);

	int err_rm1file (const char * path,
			const char * c_func, const char * c_file, int line_c);

	int err_rm1dir (const char * path,
			const char * c_func, const char * c_file, int line_c);

	void * err_ckalloc (size_t n_elem, size_t size_elem,
			const char * c_func, const char * c_file, int line_c);

	void * err_ckmalloc (size_t size,
			const char * c_func, const char * c_file, int line_c);

	void * err_ckrealloc (void * old, size_t new_size,
			const char * c_func, const char * c_file, int line_c);

//	void * err_ckaligned_malloc (size_t alignment, size_t size,
//		const char * c_func, const char * c_file, int line_c);

#ifdef __cplusplus
}
#endif

/**********************************************************
 ********************* Log Functions **********************
 **********************************************************/

#ifdef __cplusplus
extern "C" {
#endif

	void err_mesg (const char * format, ...);

	void warn_mesg (const char * format, ...);

#ifdef __cplusplus
}
#endif

/**********************************************************
 ************* String Manipulation Functions **************
 **********************************************************/

#ifdef __cplusplus
extern "C" {
#endif

	int err_str_cmp_tail (const char * s1, const char * s2, int len,
			const char * c_func, const char * c_file, int line_c);

	int err_str_cut_tail (char * str, const char * tail,
			const char * c_func, const char * c_file, int line_c);

	int err_chomp (char * str,
			const char * c_func, const char * c_file, int line_c);

	int err_chop (char * str, char c,
			const char * c_func, const char * c_file, int line_c);

	char ** err_split_string (char * string, char * markers, int * item_c,
			const char * c_func, const char * c_file, int line_c);

	int err_is_empty_line (char * line,
			const char * c_func, const char * c_file, int line_c);

#ifdef __cplusplus
}
#endif

/**********************************************************
 ******************** Misc Functions **********************
 **********************************************************/

#ifdef __cplusplus
extern "C" {
#endif

	char * err_get_abs_path (const char * path,
			const char * c_func, const char * c_file, int line_c);

	char * err_get_program_path (const char * program,
			const char * c_func, const char * c_file, int line_c);

	int cksystem (const char * cmd);

	uint64_t next_prime (uint64_t num);

	int int2deci_nbits (int num);

	void time_diff (struct timeval * beg, struct timeval * end, struct timeval * diff);

	void time_div (struct timeval * all, struct timeval * unit, int n_part);

  int ckpthread_create (pthread_t * pid, const pthread_attr_t * attr,
      void* (*start_routine)(void*), void * arg);
  int ckpthread_join (pthread_t pid);

#ifdef __cplusplus
}
#endif

#endif
