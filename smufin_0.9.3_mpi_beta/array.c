/*
 * array.c
 */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "gmem.h"
#include "stack.h"
#include "trace.h"

#include "array.h"
#include "array_private.h"

/********************************************************************/

#define MIN(a,b) ((a)<(b)?(a):(b))

/********************************************************************/

/* Helper functions */

/***************************************************************
 Lexicographical order for pairs
****************************************************************/
static inline int leq_pair(int a1, int a2, int b1, int b2)
{
  return (a1 < b1 || (a1 == b1 && a2 <= b2));
}

/***************************************************************
 Lexicographical order for trippels
****************************************************************/
static inline int leq_triple(int a1, int a2, int a3, int b1, int b2, int b3)
{
  return (a1 < b1 || (a1 == b1 && leq_pair(a2,a3, b2,b3)));
}

/***************************************************************
 Stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
****************************************************************/
static void radixPass(const int *a, int *b, const int *r, int n, int K) 
{
  int i,sum;
  int *c = gsuffix_malloc((K+1)*sizeof(int));  // counter array

  // count occurrences
  for (i = 0;  i <= K;  i++) c[i] = 0;     // reset counters
  for (i = 0;  i < n;  i++) c[r[a[i]]]++;  // count occurrences
  for (i = 0, sum = 0;  i <= K;  i++)      // exclusive prefix sums
  {
     int t = c[i];  c[i] = sum;  sum += t;
  }
  for (i = 0;  i < n;  i++) b[c[r[a[i]]]++] = a[i];      // sort

  gsuffix_free(c);
}

/********************************************************************/

/* Core algorithms */

/***************************************************************
 Constructs the suffix array SA of s[0..n-1] in {1..K}^n.
 Requires that s[n]=s[n+1]=s[n+2]=0, n>=2
 
 Algorithm based on "Linear Work Suffix Array Constuction"
****************************************************************/
static void suffixArray(const int *s, int *SA, int n, int K)
{
  int i,j,t,k,p;

  int n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n12=n0+n2; /* was orginally called n02 using n1 + n2 doesn't work */

  int *s12  = gsuffix_malloc((n12 + 3)*sizeof(int)); s12[n12] = s12[n12+1] = s12[n12+2]=0; 
  int *SA12 = gsuffix_malloc((n12 + 3)*sizeof(int)); SA12[n12] = SA12[n12+1] = SA12[n12+2]=0;
  int *s0   = gsuffix_malloc(n0*sizeof(int));
  int *SA0  = gsuffix_malloc(n0*sizeof(int));

  //******* Step 0: Construct sample *******
  // generate positions of mod 1 and mod  2 suffixes
  // the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
  for (i=0, j=0;  i < n+(n0-n1);  i++)
  {
  	if (i%3 != 0)
  		s12[j++] = i;
  }

  //******* Step 1: Sort sample suffixes *******
  // lsb radix sort the mod 1 and mod 2 triples
  radixPass(s12 , SA12, s+2, n12, K);
  radixPass(SA12, s12 , s+1, n12, K);  
  radixPass(s12 , SA12, s  , n12, K);

  // find lexicographic names of triples
  int name = 0, c0 = -1, c1 = -1, c2 = -1;
  for (i = 0;  i < n12;  i++) {
    if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) { 
      name++;  c0 = s[SA12[i]];  c1 = s[SA12[i]+1];  c2 = s[SA12[i]+2];
    }
    if (SA12[i] % 3 == 1) { s12[SA12[i]/3]      = name; } // left half
    else                  { s12[SA12[i]/3 + n0] = name; } // right half
  }

  // recurse if names are not yet unique
  if (name < n12) {
    suffixArray(s12, SA12, n12, name);
    // store unique names in s12 using the suffix array 
    for (i = 0;  i < n12;  i++) s12[SA12[i]] = i + 1;
  } else // generate the suffix array of s12 directly
    for (i = 0;  i < n12;  i++) SA12[s12[i] - 1] = i; 

  //******* Step 2: Sort nonsample suffixes *******
  // stably sort the mod 0 suffixes from SA12 by their first character
  for (i=0, j=0;  i < n12;  i++) if (SA12[i] < n0) s0[j++] = 3*SA12[i];
  radixPass(s0, SA0, s, n0, K);

  //******* Step 3: Merge *******
  // merge sorted SA0 suffixes and sorted SA12 suffixes
  for (p=0,  t=n0-n1,  k=0;  k < n;  k++)
  {
#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
    int i = GetI(); // pos of current offset 12 suffix
    int j = SA0[p]; // pos of current offset 0  suffix
    if (SA12[t] < n0 ? 
        leq_pair(s[i],       s12[SA12[t] + n0], s[j],       s12[j/3]) :
        leq_triple(s[i],s[i+1],s12[SA12[t]-n0+1], s[j],s[j+1],s12[j/3+n0]))
    { // suffix from SA12 is smaller
      SA[k] = i;  t++;
      if (t == n12)
      {
      	// done --- only SA0 suffixes left
        for (k++;  p < n0;  p++, k++) SA[k] = SA0[p];
      }
    } else
    { 
      SA[k] = j;  p++; 
      if (p == n0) 
      { 
      	// done --- only SA12 suffixes left
        for (k++;  t < n12;  t++, k++)
        {
        	SA[k] = GetI(); 
        }
      }
    }  
  }
  gsuffix_free(s12);
  gsuffix_free(SA12);
  gsuffix_free(s0);
  gsuffix_free(SA0);
}

/***************************************************************
 Constructs the lcp array belonging to the suffix array
 SA and the text s.
 
 Algorithm based on "Linear-Time Longest-Common-Prefix
 Computation in Suffix Arrays and Its Applications" by
 Kasai et al.
****************************************************************/
static void lcpArray(const int SA[], const int s[], int lcp[], int n)
{
	int h,i,k;
	int *rank;

	/* Allocate auxiliary space for the suffix's ranks */
	rank = (int*)gsuffix_malloc(sizeof(int)*n);

	/* Build inverse suffix array (the ranks of the suffixes) */
	for (i=0;i<n;i++)
		rank[SA[i]] = i;

	/* Is by definition 0 and not touched in the following loop */
	lcp[0] = 0 ;

	/* now build the lcp array */
	for (i=0,h=0;i<n;i++)
	{
		int r = rank[i];

		if (rank[i] > 0)
		{
			k = SA[r-1];

			while (s[i+h] == s[k+h])
				h++;

			lcp[r] = h;
			if (h>0) h--;
		}
	}

	gsuffix_free(rank);
}

struct childTabEntry
{
	int up;
	int down;
	int nextIndex;
};

/***************************************************************
 Constructs the up and down values which can be used for
 simulating top down traversals.

 Implemented algorithm is based on algorithm described in
 "Replacing suffix trees with enhanced suffix arrays" by
 Abouelhoda, Kurtz, Ohlebusch.
 
 What is done is basically a simulation of a bottom up
 traversal within a suffix tree by only using the lcp array
****************************************************************/
static void childArray(const int lcp [], struct childTabEntry ct[], int n)
{
	int i,top,lastIndex = -1;
	struct stack *st;

	st = stack_new();
	stack_push(st,0);

	for (i=0;i<n;i++)
	{
		top = stack_top(st);

		while (lcp[i] < lcp[top])
		{
			lastIndex = stack_pop(st);
			top = stack_top(st);

			if (lcp[i] <= lcp[top] && (lcp[top] != lcp[lastIndex]))
				ct[top].down = lastIndex;
		}

		/* now lcptab[i] >= lcptab[top] */
		if (lastIndex != -1)
		{
			ct[i].up = lastIndex;
			lastIndex = -1;
		}
		stack_push(st,i);
	}

	stack_delete(st);

	/* Calculate the nextLIndex values */
	st = stack_new();
	stack_push(st,0);

	for (i=0;i<n;i++)
	{
		while (lcp[i] < lcp[stack_top(st)])
			stack_pop(st);
		
		if (lcp[i] == lcp[stack_top(st)])
			ct[stack_top(st)].nextIndex = i;

		stack_push(st,i);
	}

	stack_delete(st);
}


/********************************************************************/

typedef struct {int l; int r;} interval_t;

/***************************************************************
 Create the extended suffix array. K means the maximum value
 a character in S has. If K == -1, it is determined by the
 function. Note that n have to be a multiple of 3!
****************************************************************/
struct exSuffixArray *exSuffixArray_create(const int s[], int n, int K)
{
	struct exSuffixArray *es;
	struct childTabEntry *ct;
	int *sa, *lcp;

	/** BEGIN_DEBUG **/ 
	ENTER("s=%p,n=%d,K=%d",s,n,K);
	SHOWINTSTRING(LOG4C_PRIORITY_TRACE,s,0,n);
	/** END_DEBUG **/ 

	/* Determine K if requested */
	if (K == -1)
	{
		int i;

		for (i=0;i<n;i++)
			if (K < s[i]) K = s[i];

		/** BEGIN_DEBUG **/ 
		LOG(LOG4C_PRIORITY_TRACE,"K will be %d",K);
		/** END_DEBUG **/ 
	}

	/* allocate necessary memory */
	sa = gsuffix_malloc((n+3)*sizeof(int));
	lcp = gsuffix_malloc((n+3)*sizeof(int));
	ct = gsuffix_malloc(sizeof(ct[0])*(n+3));

	memset(sa,0,sizeof((n+3)*sa[0]));
	memset(lcp,0,sizeof((n+3)*lcp[0]));
	memset(ct,0xff,sizeof(ct[0])*(n+3));

	/** BEGIN_DEBUG **/ 
	LOG(LOG4C_PRIORITY_TRACE,"Creating suffix array",K);
	/** END_DEBUG **/ 

	suffixArray(s,sa,n,K);
	
	/** BEGIN_DEBUG **/ 
	LOG(LOG4C_PRIORITY_TRACE,"Creating lcp array",K);
	/** END_DEBUG **/ 

	lcpArray(sa,s,lcp,n);
	
	/** BEGIN_DEBUG **/ 
	LOG(LOG4C_PRIORITY_TRACE,"Creating child array",K);
	/** END_DEBUG **/ 

	childArray(lcp,ct,n);
	ct[n].up = ct[0].nextIndex; /* to make getInterval() work also for [0,n] interval */

	es = gsuffix_malloc(sizeof(*es));
	es->s = s;
	es->sa = sa;
	es->lcp = lcp;
	es->ct = ct;
	es->n = n;

	/** BEGIN_DEBUG **/ 
	LEAVE("RC=%p",es);
	/** END_DEBUG **/ 

	return es;
}

/***************************************************************
 Complete free the resources of the given extendend suffix
 array.
****************************************************************/
void exSuffixArray_delete(struct exSuffixArray *ex)
{
	ENTER("ex=%p",ex);

	if (!ex) return;

	gsuffix_free((void*)ex->lcp);
	gsuffix_free((void*)ex->sa);
	gsuffix_free((void*)ex->ct);
	gsuffix_free(ex);
	
	LEAVE("");
}

/***************************************************************
 Determine the lcp of the suffixes in the given interval
****************************************************************/
static int exSuffixArray_getlcp(const struct exSuffixArray *ex, const interval_t lcp_interval)
{
	int i, j;

	const struct childTabEntry *ct = ex->ct;
	const int *lcp = ex->lcp;

	i = lcp_interval.l;
	j = lcp_interval.r;

	if (i < ct[j+1].up && ct[j+1].up <= j) return lcp[ct[j+1].up];
	return lcp[ct[i].down];
}

/***************************************************************
 Determine the embedded lcp interval for the given character c
****************************************************************/
static interval_t exSuffixArray_getInterval(struct exSuffixArray *ex, const interval_t lcp_interval, int c, int l)
{
	int i,j,i1,i2;

	interval_t child_interval;

	ENTER("ex=%p,lcp_interval.l=%d,lcp_interval.r=%d,c=%d,l=%d",ex,lcp_interval.l,lcp_interval.r,c,l);

	i = lcp_interval.l;
	j = lcp_interval.r;

	if (i < ex->ct[j+1].up && ex->ct[j+1].up <= j) i1 = ex->ct[j+1].up;
	else i1 = ex->ct[i].down;

//	printf("i1 = %d\n",i1);

	/* Child interval is [i..i-1]
	 * 
	 * s[SA[i]+l] contains the branching characters,
	 * Note that s[SA[i]+l] = s[SA[j]+l] for every  i<j<=i-1
	 */
	if (ex->s[ex->sa[i]+l] == c)
	{
		child_interval.l = i;
		child_interval.r = i1 - 1;

		assert(child_interval.l <= child_interval.r);
		assert(child_interval.l >= i && child_interval.r <= j);

		LEAVE("l=%d,r=%d",child_interval.l,child_interval.r);
		return child_interval;
	}

//	printf("%c\n",ex->s[ex->sa[i]+l]);
//	printf("i1=%d, i2=%d\n",i1,i2);

	while (ex->ct[i1].nextIndex != -1)
	{
		i2 = ex->ct[i1].nextIndex;

		/* Child interval is [i1..i2-1] */
		if (ex->s[ex->sa[i1]+l] == c)
		{
			child_interval.l = i1;
			child_interval.r = i2 - 1;
			assert(child_interval.l <= child_interval.r);
			assert(child_interval.l >= i && child_interval.r <= j);

			LEAVE("l=%d,r=%d",child_interval.l,child_interval.r);
			return child_interval;
		}

		i1 = i2;
	}

	if (ex->s[ex->sa[i1]+l] == c)
	{
		child_interval.l = i1;
		child_interval.r = j;
		assert(child_interval.l <= child_interval.r);
		assert(child_interval.l >= i && child_interval.r <= j);

		LEAVE("l=%d,r=%d",child_interval.l,child_interval.r);
		return child_interval;
	}

	/* Nothing found */
	child_interval.l = -1;
	child_interval.r = -1;

	LEAVE("l=%d,r=%d",child_interval.l,child_interval.r);
	return child_interval;
}

/***************************************************************
 Searches substring p with length m within the string
 represented by the extended suffix array.

 Returns the number of hits (also if the callback cancels the
 search)
****************************************************************/
int exSuffixArray_subString_callback(struct exSuffixArray *ex, const int p[], int m, int (*hitcallback)(int pos, void *userdata), void *userdata)
{
	interval_t ij;
	int c;

	/** BEGIN_DEBUG **/ 
	ENTER("ex=%p,p=%p,m=%d,cb=%p,userdata=%p",ex,p,m,hitcallback,userdata);
	SHOWINTSTRING(LOG4C_PRIORITY_TRACE,p,0,m);
	/** END_DEBUG **/ 

	c = 0;

	ij.l = 0;
	ij.r = ex->n - 1;
	ij = exSuffixArray_getInterval(ex,ij,p[c],exSuffixArray_getlcp(ex,ij));

	while (ij.l != -1 && c < m)
	{
		if (ij.l != ij.r)
		{
			int l = exSuffixArray_getlcp(ex,ij);
			int min = MIN(l,m);
			int j;

			for (j=c;j<min;j++)
			{
				if (ex->s[ex->sa[ij.l]+j] != p[j])
				{
					LEAVE("not found");
					return 0;
				}
			}

			c = min;

			if (c == m) break;
			ij = exSuffixArray_getInterval(ex,ij,p[c],exSuffixArray_getlcp(ex,ij));
			if (ij.l == -1)
			{
				LEAVE("not found");
				return 0;
			}
		} else
		{
			int j;

			/* Found a singleton interval */

			for (j=c;j<m;j++)
			{
				if (ex->s[ex->sa[ij.l]+j] != p[j])
				{
					LEAVE("not found");
					return 0;
				}
			}

			break;
		}
	}
	
	if (ij.l == -1)
	{
		LEAVE("not found");
		return 0;
	}

	/* All occurrences of the substring can be found in the current
	 * lcp interval */
	if (hitcallback)
	{
		int j;

		/** BEGIN_DEBUG **/ 
		LOG(LOG4C_PRIORITY_TRACE,"found lcp interval (%d,%d)",ij.l,ij.r);
		/** END_DEBUG **/ 
		
		for (j=ij.l;j<=ij.r;j++)
		{
			/** BEGIN_DEBUG **/ 
			LOG(LOG4C_PRIORITY_TRACE,"found string at pos %d",ex->sa[j]);
			/** END_DEBUG **/ 

			/* Returns 0, in case search should be aborted */
			if (!hitcallback(ex->sa[j],userdata))
				break;
		}	
	}
	LEAVE("found %d occurrences",ij.r - ij.l + 1);
	return ij.r - ij.l + 1;
}

/***************************************************************
 The basic callback for exSuffixArray_subString
****************************************************************/
static int exSuffixArray_subString_basic_callback(int pos, void *userdata)
{
	int *ptr = (int*)userdata;
	*ptr = pos;
	
	/* We are only interested in the first occurrence */
	return 0;
}

/***************************************************************
 Searches substring p with length m within the string
 represented by the extended suffix array
****************************************************************/
int exSuffixArray_subString(struct exSuffixArray *ex, const int p[], int m)
{
	int pos = -1;

	exSuffixArray_subString_callback(ex,p,m,exSuffixArray_subString_basic_callback, &pos);

	return pos;
}

/*******************************************************************/

/* Some questions which need further research:
 *   Can the child table be constructed without the previous
 *   construction of the lcp array?
 *
 *   Can Ukkonen (or any other suffix tree build algorithm) be changed
 *   to build up the child table array instead of whole suffix tree?
 */
