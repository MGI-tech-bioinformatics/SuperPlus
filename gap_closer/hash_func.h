/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-11-02 18:29:34
  *Edit History: 
***********************************************************/

#ifndef XDK_HASH_FUNC_H
#define XDK_HASH_FUNC_H

void hash_func_init (void); // not thread safte
void hash_func_free (void); // not thread safte

// Blizzard hash function
uint64_t blizzard_hash_func (const char * key, int key_len, int dwHashType); // thread safe

#endif
