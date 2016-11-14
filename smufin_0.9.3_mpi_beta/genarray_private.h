#ifndef SUFFIX__GENARRAY_PRIVATE_H
#define SUFFIX__GENARRAY_PRIVATE_H

#ifndef SUFFIX__ARRAY_H
#include "array.h"
#endif

#include "array8.h"

#define IDXT_BINARY 0
#define IDXT_DIRECT 1

struct GenArray
{
	int *concat;				/* The concatenated string */
	unsigned char *concat8;
	int concat_len;				/* The length of the concatenated string */

	struct exSuffixArray *ex;	/* The extented suffix array */
	struct exSuffixArray8 *ex8;	/* The extented suffix array */

	int n;						/* Number of strings */

	int index_type; 				/* eighter IDXT_BINARY or IDXT_DIRECT */
	int *pos_array;					/* the position of the strings (contains n elements) */
	int *idx_array;					/* the indicies of the strings (contains concat_len elements), maybe NULL */

	int (*hitcallback)(int,int,void*);	/* temporary storage for the hitcallback */
	void *userdata;						/* temporary storage for the userdata */
	int hits;							/* temporary storgae for the hit counter */
};

#endif /*SUFFIX__GENARRAY_PRIVATE_H*/
