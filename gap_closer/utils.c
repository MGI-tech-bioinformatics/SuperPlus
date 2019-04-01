/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2017-05-13 19:52:51
  *Edit History: 
***********************************************************/

#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <stdio.h>
#include <fcntl.h>
#include <dirent.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "utils.h"

/**********************************************************
 *************** IO, FileSystem Functions *****************
 **********************************************************/

FILE *
err_ckopen (const char * file, const char * mode,
		const char * c_func, const char * c_file, int line_c)
{
	FILE * fp;

	if (strcmp(file, "-") == 0)
		return (strstr(mode,"r") ? stdin : stdout);

	if ((fp = fopen(file,mode)) == NULL)
		err_mesg ("[%s] fail to open file '%s' in file %s at line %d!",
				c_func, file, c_file, line_c);

	return fp;
}

FILE *
err_ckpopen (const char * cmd, const char * mode,
		const char * c_func, const char * c_file, int line_c)
{
	FILE * fp;

	if ((fp = popen(cmd,mode)) == NULL)
		err_mesg ("[%s] fail to run cmd '%s' in file %s at line %d!",
				c_func, cmd, c_file, line_c);

	return fp;
}

int
err_xopen (const char * file, int flags,
		const char * c_func, const char * c_file, int line_c)
{
	int fd;

	if ((fd = open(file,flags)) < 0)
		err_mesg ("[%s] fail to open file '%s' in file %s at line %d!",
				c_func, file, c_file, line_c);

	return fd;
}

size_t
err_ckfwrite (const void * ptr, size_t size, size_t n_elem, FILE * fp,
		const char * c_func, const char * c_file, int line_c)
{
	size_t ret;

	if ((ret = fwrite(ptr,size,n_elem,fp)) != n_elem)
		err_mesg ("[%s] fail to call fwrite in file %s at line %d: %s!",
				c_func, c_file, line_c, strerror(errno));

	return ret;
}

int
err_ckfflush (FILE * fp,
		const char * c_func, const char * c_file, int line_c)
{
	int ret;

	if ((ret = fflush(fp)) != 0)
		err_mesg ("[%s] fail to call fflush in file %s at line %d: %s!",
				c_func, c_file, line_c, strerror(errno));

#ifdef FSYNC_ON_FLUSH
	struct stat st;
	if (fstat(fileno(fp),&st) != 0)
		err_mesg ("[%s] fail to call fstat in file %s at line %d: %s!",
				c_func, c_file, line_c, strerror(errno));

	if (S_ISREG(st.st_mode)
			&& fsync(fileno(fp)) != 0)
		err_mesg ("[%s] fail to call fsync in file %s at line %d: %s!",
				c_func, c_file, line_c, strerror(errno));
#endif

	return ret;
}

DIR *
err_ckopendir (const char * path,
		const char * c_func, const char * c_file, int line_c)
{
	DIR * dp;

	if ((dp = opendir(path)) == NULL)
		err_mesg ("[%s] fail to open directory '%s' in file %s at line %d!",
				c_func, path, c_file, line_c);

	return dp;
}

int
err_ckcreate_dir (const char * path,
		const char * c_func, const char * c_file, int line_c)
{
	if (access(path,0)!=0 && mkdir(path,0755)!=0)
		err_mesg ("[%s] fail to create directory '%s' in file %s at line %d!",
				c_func, path, c_file, line_c);

	return 0;
}

int
err_rm1file (const char * path,
		const char * c_func, const char * c_file, int line_c)
{
	if (unlink(path) != 0)
		err_mesg ("[%s] fail to remove file '%s' in file %s at line %d!",
				c_func, path, c_file, line_c);

	return 0;
}

int
err_rm1dir (const char * path,
		const char * c_func, const char * c_file, int line_c)
{
	char cmd[4096];

	sprintf (cmd, "rm -rf %s", path);
	cksystem (cmd);

	return 0;
}

void *
err_ckalloc (size_t n_elem, size_t size_elem,
		const char * c_func, const char * c_file, int line_c)
{
	void * ptr;

	//printf ("n: %ld\tsize: %ld\n", n_elem, size_elem);
	if ((ptr = calloc(n_elem,size_elem)) == NULL)
		err_mesg ("[%s] fail to alloc space in file %s at line %d!",
				c_func, c_file, line_c);

	return ptr;
}

void *
err_ckmalloc (size_t size,
		const char * c_func, const char * c_file, int line_c)
{
	void * ptr;

	if ((ptr = malloc(size)) == NULL)
		err_mesg ("[%s] fail to alloc space in file %s at line %d!",
				c_func, c_file, line_c);

	return ptr;
}

void *
err_ckrealloc (void * old, size_t new_size,
		const char * c_func, const char * c_file, int line_c)
{
	void * ptr;

	if ((ptr = realloc(old,new_size)) == NULL)
		err_mesg ("[%s] fail to realloc space in file %s at line %d!",
				c_func, c_file, line_c);

	return ptr;
}

/*
void *
err_ckaligned_malloc (size_t alignment, size_t size,
		const char * c_func, const char * c_file, int line_c)
{
	void * ptr;

	if ((ptr = aligned_alloc(alignment,size)) == NULL)
		err_mesg ("[%s] fail to alloc aligned memory in file %s at line %d!",
				c_func, c_file, line_c);

	return ptr;
}
*/

/**********************************************************
 ********************* Log Functions **********************
 **********************************************************/

void
err_mesg (const char * format, ...)
{
	static char buf[LINE_MAX];
	va_list arg;

	va_start (arg, format);
	vsnprintf (buf, LINE_MAX, format, arg);
	va_end (arg);

	fprintf (stderr, "Error: %s\n", buf);

	abort ();
}

void
warn_mesg (const char * format, ...)
{
	static char buf[LINE_MAX];
	va_list arg;

	va_start (arg, format);
	vsnprintf (buf, LINE_MAX, format, arg);
	va_end (arg);

	fprintf (stderr, "Warning: %s\n", buf);
}

/**********************************************************
 ************* String Manipulation Functions **************
 **********************************************************/

int
err_str_cmp_tail (const char * s1, const char * s2, int len,
		const char * c_func, const char * c_file, int line_c)
{
	int l1, l2;
	int ret;

	l1 = strlen (s1);
	l2 = strlen (s2);

	if (l1<len || l2<len)
		return -1;

	ret = strncmp (s1+l1-len, s2+l2-len, len);

	// make sure return value -1 is for invalid compare
	if (ret == -1)
		ret = -2;

	return ret;
}

int
err_str_cut_tail (char * str, const char * tail,
		const char * c_func, const char * c_file, int line_c)
{
	int l_str;
	int l_tail;

	l_str = strlen (str);
	l_tail = strlen (tail);

	if (l_str < l_tail
			|| strncmp(str+l_str-l_tail,tail,l_tail) != 0)
		return -1;

	*(str+l_str-l_tail) = '\0';

	return 0;
}

int
err_chomp (char * str,
		const char * c_func, const char * c_file, int line_c)
{
	return err_chop (str, '\n', c_func, c_file, line_c);
}

int
err_chop (char * str, char c,
		const char * c_func, const char * c_file, int line_c)
{
	int l;

	if ((l = strlen(str)) <= 0)
		return -1;

	if (str[l-1] == c) {
		str[l-1] = '\0';
		return 0;
	}

	return 1;
}

char **
err_split_string (char * string, char * markers, int * item_c,
		const char * c_func, const char * c_file, int line_c)
{
	char * ch;
	char * str_dup;
	char ** items;
	int count;

	str_dup = strdup (string);
	if ((ch = strtok(str_dup, markers)) == NULL)
		return NULL;

	count = 1;
	while ((ch = strtok(NULL, markers)) != NULL)
		count++;

	items = (char **) ckalloc (count, sizeof(char *));

	count = 0;
	items[count++] = string;
	ch = strtok (string, markers);
	while ((ch = strtok(NULL, markers)) != NULL)
		items[count++] = ch;

	*item_c = count;

	free (str_dup);

	return items;
}

int
err_is_empty_line (char * line,
		const char * c_func, const char * c_file, int line_c)
{
	char * ch;

	if (*line == '#')
		return 1;

	for (ch=line; *ch!='\0'; ch++)
		if (!isspace(*ch))
			break;
	if (*ch == '\0')
		return 1;

	return 0;
}

/**********************************************************
 ******************** Misc Functions **********************
 **********************************************************/

char *
err_get_abs_path (const char * path,
			const char * c_func, const char * c_file, int line_c)
{
	char * c1, * c2, * c3;
	char * cwd, * abs_path;
	int len;

	abs_path = ALLOC_PATH;

	if (*path == '/')
		strcpy (abs_path, path);
	else {
		cwd = ALLOC_PATH;
		if ((cwd = getcwd(cwd,PATH_MAX)) == NULL)
			err_mesg ("[%s] fail to get current working directory in file %s at line %d!",
					c_func, c_file, line_c);
		sprintf (abs_path, "%s/%s", cwd, path);
		free (cwd);
	}

	// merge '//' to '/'
	for (c1 = abs_path; *(c1+1)!='\0'; c1++) {
		if (*c1=='/' && *(c1+1)=='/') {
			for (c2=c1+1; *c2!='\0'; c2++)
				*(c2-1) = *c2;
			*(c2-1) = '\0';
		}
	}

	// merge '/xxx/../' to '/'
	c1 = abs_path;
	while ((c1 = strstr(abs_path, "/../")) != NULL) {
		for (c3=c1-1; c3>=abs_path; c3--)
			if (*c3 == '/')
				break;

		if (c3 < abs_path)
			err_mesg ("[%s] invalid path '%s' in file %s at line %d!",
					c_func, path, c_file, line_c);

		len = c1 + 3 - c3;
		for (c2=c1+3; *c2!='\0'; c2++)
			*(c2-len) = *c2;
		*(c2-len) = '\0';
	}

	if (str_cmp_tail(abs_path,"/..",3) == 0) {
		for (c1=abs_path+strlen(abs_path)-4; c1>=abs_path; c1--)
			if (*c1 == '/')
				break;

		if (c1 < abs_path)
			err_mesg ("[%s] invalid path '%s' in file %s at line %d!",
					c_func, path, c_file, line_c);

		*c1 = '\0';
	}

	chop (abs_path, '.');
	chop (abs_path, '/');

	return abs_path;
}

char *
err_get_program_path (const char * program,
			const char * c_func, const char * c_file, int line_c)
{
	char * path;
	char * dir;
	char * env_dirs;
	DIR * dp;
	struct dirent * dirp;

	env_dirs = getenv ("PATH");

	if ((dir = strtok(env_dirs, ":")) == NULL)
		return NULL;

	path = ALLOC_PATH;

	dp = ckopendir (dir);
	while (NULL == (dirp=readdir(dp))) {
		if (strcmp(program, dirp->d_name) != 0)
			continue;
		sprintf (path, "%s/%s", dir, dirp->d_name);
		return path;
	}
	closedir (dp);

	while (NULL == (dir=strtok(NULL, ":"))) {
		dp = ckopendir (dir);
		while (NULL == (dirp=readdir(dp))) {
			if (strcmp(program, dirp->d_name))
				continue;
			sprintf (path, "%s/%s", dir, dirp->d_name);
			return path;
		}
	}

	return NULL;
}

int
cksystem (const char * cmd)
{
	int ret;

	if ((ret = system(cmd)) == -1)
		return -1;

	if (WIFEXITED(ret) && 0==WEXITSTATUS(ret))
		return 0;
	else
		return 1;
}

uint64_t
next_prime (uint64_t num)
{
	int flag;
	uint64_t i, max;

	if ((num & 0x01) == 0)
		num++;

	if (num == 1)
		num = 3;

	while (1) {
		if (num < 4)
			return num;

		max = (uint64_t) sqrt ((double)num);
		flag = 0;
		for (i=3; i<=max; i+=2) {
			if (num % i == 0) {
				flag = 1;
				break;
			}
		}

		if (!flag)
			return num;

		num += 2;
	}
}

int
int2deci_nbits (int num)
{
	int is_negative;

	if (num == 0)
		return 1;

	if (num < 0) {
		num = -num;
		is_negative = 1;
	} else
		is_negative = 0;

	return (int)log10((double)num+1e-2) + 1 + (is_negative?1:0);
}

void
time_diff (struct timeval * beg, struct timeval * end, struct timeval * diff)
{
	if (end->tv_usec < beg->tv_usec) {
		end->tv_usec += 1000000;
		end->tv_sec -= 1;
	}

	diff->tv_sec = end->tv_sec - beg->tv_sec;
	diff->tv_usec = end->tv_usec - beg->tv_usec;
}

void
time_div (struct timeval * all, struct timeval * unit, int n_part)
{
	long usec;

	if (n_part <= 0)
		err_mesg ("[%s] 'n_part' must > 0!", __func__);

	usec = all->tv_sec*1000000 + all->tv_usec;
	usec /= n_part;

	unit->tv_sec = usec / 1000000;
	unit->tv_usec = usec % 1000000;
}

int
ckpthread_create (pthread_t * pid, const pthread_attr_t * attr,
    void* (*start_routine)(void*), void * arg)
{
  if (pthread_create(pid, attr, start_routine, arg) != 0)
    err_mesg ("fail to create thread!");

  return 0;
}

int
ckpthread_join (pthread_t pid)
{
  void * status;

  if (pthread_join(pid, &status) != 0)
    err_mesg ("fail to join thread");

  return 0;
}
