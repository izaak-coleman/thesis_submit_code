#ifndef __ARRAY8_PRIVATE_H
#define __ARRAY8_PRIVATE_H

#ifdef USE_COMPRESSED_ARRAYS

struct keyval
{
	int key;
	int val;
};

struct keyval_vector
{
	struct keyval *table;
	int table_allocated;
	int table_size;
};

#endif

struct exSuffixArray8
{
	int n;						/* length of the string */
	const unsigned char *s;		/* The string */
	int *sa;					/* The suffix array */
	int *lcp;					/* The lcp array */
	struct childTabEntry *ct;	/* The child table (n+1 elements) */

#ifdef USE_COMPRESSED_ARRAYS
	struct keyval_vector ct_up_vector;
	struct keyval_vector ct_down_vector;
	struct keyval_vector ct_nextIndex_vector;
#endif

};

#endif /*__ARRAY_H*/
