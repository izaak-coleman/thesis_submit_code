#ifndef SUFFIX__HASH_GENHASH_PRIVATE_H
#define SUFFIX__HASH_GENHASH_PRIVATE_H

#ifndef SUFFIX__HASH_HASH_H
#include "hash.h"
#endif

/* The options for determining the index of the string */
#define IDXT_BINARY 0
#define IDXT_DIRECT 1

struct gen_hash
{
	char *concat;					/* The concatenated string */
	int concat_len;					/* The length of the concatenated string */

	struct kmer_hash_table *ht;		/* The extented suffix array */

	int n;							/* Number of strings */

	int index_type; 				/* eighter IDXT_BINARY or IDXT_DIRECT */
	int *pos_array;					/* the position of the strings (contains n elements) */
	int *idx_array;					/* the indicies of the strings (contains concat_len elements), maybe NULL */

	int (*hitcallback)(int,int,void*);	/* temporary storage for the user supplied hitcallback */
	void *userdata;						/* temporary storage for the userdata */
	int hits;							/* temporary storgae for the hit counter */
};

#endif /* SUFFIX__HASH_GENHASH_PRIVATE_H */
