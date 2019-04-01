/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-11-05 11:21:59
  *Edit History: 
***********************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "hash_func.h"

static uint64_t * blizzard_crypt_table;

void
hash_func_init (void)
{
	uint64_t seed = 0x00100001;
	uint64_t index1 = 0;
	uint64_t index2 = 0;
	uint64_t i;

	blizzard_crypt_table = (uint64_t *) ckalloc (0x500, sizeof(uint64_t));

	for (index1=0; index1<0x100; index1++) {
		for (index2=index1,i=0; i<5; i++,index2+=0x100) {
			uint64_t temp1, temp2;

			seed = (seed * 125 + 3) % 0x2AAAAB;
			temp1 = (seed & 0xFFFF) << 0x10;
			seed = (seed * 125 + 3) % 0x2AAAAB;
			temp2 = (seed & 0xFFFF);

			blizzard_crypt_table[index2] = (temp1 | temp2);
		}
	}
}

void
hash_func_free (void)
{
  free (blizzard_crypt_table);
}

/*----------------------------------------------------------------------------*/
/*------------------------------ Blizzard Hash -------------------------------*/
/*----------------------------------------------------------------------------*/

uint64_t
blizzard_hash_func (const char * key, int key_len, int dwHashType)
{
	uint64_t ch;
	uint64_t seed1 = 0x7FED7FED;
	uint64_t seed2 = 0xEEEEEEEE;

	while (--key_len >= 0) {
		ch = toupper (*key++);
		seed1 = blizzard_crypt_table[(dwHashType<<8) + ch] ^ (seed1 + seed2);
		seed2 = ch + seed1 + seed2 + (seed2 << 5) + 3;
	}

	return seed1;
}

/*----------------------------------------------------------------------------*/
/*------------------------------  Hash -------------------------------*/
/*----------------------------------------------------------------------------*/

