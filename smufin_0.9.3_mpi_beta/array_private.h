#ifndef __ARRAY_PRIVATE_H
#define __ARRAY_PRIVATE_H

struct exSuffixArray
{
	int n;							/* length of the string */
	const int *s;					/* The string */
	const int *sa;					/* The suffix array */
	const int *lcp;					/* The lcp array */
	const struct childTabEntry *ct;	/* The child table (n+1 elements) */
};

#endif /*__ARRAY_H*/
