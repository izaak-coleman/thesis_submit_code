#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "array8.h"
#include "gmem.h"
#include "stack.h"
#include "trace.h"

#include "ds_ssort.h"

#include "array8_private.h"

#if 0
int exSuffixArray8_enumerate(struct exSuffixArray8 *ex, int k);
#endif

/********************************************************************/

#define MIN(a,b) ((a)<(b)?(a):(b))

/********************************************************************/

#ifdef USE_COMPRESSED_ARRAYS

static void keyval_init(struct keyval_vector *vec)
{
	vec->table = gsuffix_malloc(2000*sizeof(struct keyval));
	vec->table_size = 0;
	vec->table_allocated = 2000;
}

static void keyval_add(struct keyval_vector *vec, int key, int val)
{
	if (vec->table_size == vec->table_allocated)
	{
		vec->table_allocated *= 2;
		if (!(vec->table = myrealloc(vec->table,vec->table_allocated*sizeof(struct keyval))))
		{
			fprintf(stderr,"Not enough memory for increasing vector\n");
			exit(-1);
		}
	}
	vec->table[vec->table_size].key = key;
	vec->table[vec->table_size].val = val;
	vec->table_size++;
}

static int keyval_compare_func(const void *arg1, const void *arg2)
{
	struct keyval *a1 = (struct keyval*)arg1;
	struct keyval *a2 = (struct keyval*)arg2;
	return a1->key - a2->key;
}

static void keyval_finish(struct keyval_vector *vec)
{
	qsort(vec->table,vec->table_size,sizeof(vec->table[0]),keyval_compare_func);
}

static int keyval_find(const struct keyval_vector *vec, int key)
{
	int low, high, mid;

	low = 0;
	high = vec->table_size - 1;

	while(low <= high)
	{
		mid = (low+high)/2;
		if(key < vec->table[mid].key) high = mid - 1;
		else if(key > vec->table[mid].key) low = mid + 1;
		else return vec->table[mid].val;
	}

/*
	int i;
	
	for (i=0;i<vec->table_size;i++)
	{
		if (vec->table[i].key == key)
			return vec->table[i].val;
	}*/
	return -1;
}

#endif

/********************************************************************/


struct childTabEntry
{
#ifdef USE_COMPRESSED_ARRAYS
	unsigned char up_diff;
	unsigned char down_diff;
	unsigned char nextIndex_diff;
#else
	int up;
	int down;
	int nextIndex;
#endif
};

/********************************************************************/
static inline int exSuffixArray8_get_child_up(const struct exSuffixArray8 *ex, int index)
{
	int up;

#ifdef USE_COMPRESSED_ARRAYS
	up = ex->ct[index].up_diff;
	if (up != 255) up = index - up;
	else up = keyval_find(&ex->ct_up_vector,index);
#else
	up = ex->ct[index].up;
#endif

	return up;
}

static inline int exSuffixArray8_get_child_down(const struct exSuffixArray8 *ex, int index)
{
	int down;

#ifdef USE_COMPRESSED_ARRAYS
	down = ex->ct[index].down_diff;
	if (down != 255) down = index + down;
	else down = keyval_find(&ex->ct_down_vector,index);
#else
	down = ex->ct[index].down;
#endif

	return down;
}

static void childArray(struct exSuffixArray8 *es)
{
	int i,top,lastIndex = -1;
	struct stack *st;

	int *lcp = es->lcp;
	struct childTabEntry *ct = es->ct;
	int n = es->n;

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
			{
#ifdef USE_COMPRESSED_ARRAYS
				int diff = lastIndex - top;

				if (diff >= 0 && diff < 255)
					ct[top].down_diff = diff;
				else
				{
					ct[top].down_diff = 255;
					keyval_add(&es->ct_down_vector,top,lastIndex);
				}
#else
				ct[top].down = lastIndex;
#endif
			}
		}

		if (lastIndex != -1)
		{
#ifdef USE_COMPRESSED_ARRAYS
			int diff = i - lastIndex;

			if (diff >= 0 && diff < 255)
				ct[i].up_diff = diff;
			else
			{
				ct[i].up_diff = 255;
				keyval_add(&es->ct_up_vector,i,lastIndex);
			}
#else
			ct[i].up = lastIndex;
#endif
			lastIndex = -1;
		}
		stack_push(st,i);
	}

	stack_delete(st);

	st = stack_new();
	stack_push(st,0);

	for (i=0;i<n;i++)
	{
		while (lcp[i] < lcp[stack_top(st)])
			stack_pop(st);
		
		if (lcp[i] == lcp[stack_top(st)])
		{
			int top = stack_top(st);
#ifdef USE_COMPRESSED_ARRAYS
			int diff = i - top;
			if (diff >= 0 && diff < 254)
				ct[top].nextIndex_diff = diff;
			else
			{
				ct[top].nextIndex_diff = 254;
				keyval_add(&es->ct_nextIndex_vector,top,i);
			}
#else
			ct[top].nextIndex = i;
#endif
		}

		stack_push(st,i);
	}

	stack_delete(st);
}


static void lcpArray8(const int SA[], const unsigned char s[], int lcp[], int n)
{
	int h,i,k;
	int *rank;

	rank = (int*)gsuffix_malloc(sizeof(int)*n);

	for (i=0;i<n;i++)
		rank[SA[i]] = i;

	lcp[0] = 0 ;

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

#if 0
static void lcpArray8_space_saving(const int SA[], const char s[], int lcp[], int n)
{
  int i,h,j,k,nextk=-1;
  int *rn;

  rn = lcp;
  k = _bw_sa2ranknext(t, n, sa, occ, rn); // k is the rank of s
  for(h=i=0; i<n; i++,k=nextk) {
    assert(k>0 && i==sa[k]); 
    nextk=rn[k];          // read nextk before it is overwritten
    if(k>1) {
      // compute lcp between suffixes of rank k and k-1 (recall i==sa[k])
      j = sa[k-1];
      while(i+h<n && j+h<n && t[i+h]==t[j+h])
	h++;
      lcp[k]=h;           // h is the lcp, write it to lcp[k];  
    }
    if(h>0) h--;
  }
  assert(nextk==0);        // we have reached the last suffix s[n-1]
  return lcp;
}
#endif


typedef struct {int l; int r;} interval_t;

struct exSuffixArray8 *exSuffixArray8_create(unsigned char s[], int n)
{
	struct exSuffixArray8 *es;
	struct childTabEntry *ct;
	int *sa, *lcp;

	/** BEGIN_DEBUG **/ 
	ENTER("s=%p,n=%d,K=%d",s,n,K);
	/** END_DEBUG **/ 

	sa = gsuffix_malloc(sizeof(sa[0])*(n+3));
	lcp = gsuffix_malloc(sizeof(lcp[0])*(n+3));
	ct = gsuffix_malloc(sizeof(ct[0])*(n+3));

	memset(sa,0,sizeof(sa[0])*(n+3));
	memset(lcp,0,sizeof(lcp[0])*(n+3));
	memset(ct,0xff,sizeof(ct[0])*(n+3));

	/** BEGIN_DEBUG **/ 
	LOG(LOG4C_PRIORITY_TRACE,"Creating suffix array",K);
	/** END_DEBUG **/ 

	ds_ssort(s,sa,n);

	/** BEGIN_DEBUG **/ 
	LOG(LOG4C_PRIORITY_TRACE,"Creating lcp array",K);
	/** END_DEBUG **/ 

	lcpArray8(sa,s,lcp,n);
	
	/** BEGIN_DEBUG **/ 
	LOG(LOG4C_PRIORITY_TRACE,"Creating child array",K);
	/** END_DEBUG **/ 

	es = gsuffix_malloc(sizeof(*es));
	es->s = s;
	es->sa = sa;
	es->lcp = lcp;
	es->ct = ct;
	es->n = n;

#ifdef USE_COMPRESSED_ARRAYS
	keyval_init(&es->ct_up_vector);
	keyval_init(&es->ct_down_vector);
	keyval_init(&es->ct_nextIndex_vector);
#endif

	childArray(es);

#ifdef USE_COMPRESSED_ARRAYS
	keyval_finish(&es->ct_up_vector);
	keyval_finish(&es->ct_down_vector);
	keyval_finish(&es->ct_nextIndex_vector);
#endif

#ifdef USE_COMPRESSED_ARRAYS
	{
		int diff, nextIndex;

		nextIndex = es->ct[0].nextIndex_diff;
		if (nextIndex != 254) nextIndex += 0;
		else nextIndex = keyval_find(&es->ct_nextIndex_vector,0);

		diff = n - nextIndex;
		if (diff >= 0 && diff < 255)
			ct[n].up_diff = diff;
		else
		{
			keyval_add(&es->ct_up_vector,n,nextIndex);
			ct[n].up_diff = 255;
		}
	}
#else
	ct[n].up = ct[0].nextIndex;
#endif


	/** BEGIN_DEBUG **/ 
	LEAVE("RC=%p",es);
	/** END_DEBUG **/ 

	return es;
}

void exSuffixArray8_delete(struct exSuffixArray8 *ex)
{
	ENTER("ex=%p",ex);

	if (!ex) return;

	gsuffix_free((void*)ex->lcp);
	gsuffix_free((void*)ex->sa);
	gsuffix_free((void*)ex->ct);
	gsuffix_free(ex);
	
	LEAVE("");
}

static int exSuffixArray8_getlcp(const struct exSuffixArray8 *ex, const interval_t lcp_interval)
{
	int i, j;

	const int *lcp = ex->lcp;

	int up;
	int down;

	i = lcp_interval.l;
	j = lcp_interval.r;

	up = exSuffixArray8_get_child_up(ex,j+1);

	if (i < up && up <= j)
		return lcp[up];

	down = exSuffixArray8_get_child_down(ex,i);

	return lcp[down];
}


static interval_t exSuffixArray8_getInterval(struct exSuffixArray8 *ex, const interval_t lcp_interval, int c, int l)
{
	int i,j,i1,i2;
	int up, nextIndex;

	interval_t child_interval;

	ENTER("ex=%p,lcp_interval.l=%d,lcp_interval.r=%d,c=%d,l=%d",ex,lcp_interval.l,lcp_interval.r,c,l);

	i = lcp_interval.l;
	j = lcp_interval.r;

	up = exSuffixArray8_get_child_up(ex,j+1);

	if (i < up && up <= j) i1 = up;
	else i1 = exSuffixArray8_get_child_down(ex,i);

	if (ex->s[ex->sa[i]+l] == c)
	{
		child_interval.l = i;
		child_interval.r = i1 - 1;

		assert(child_interval.l <= child_interval.r);
		assert(child_interval.l >= i && child_interval.r <= j);

		LEAVE("l=%d,r=%d",child_interval.l,child_interval.r);
		return child_interval;
	}

#ifdef USE_COMPRESSED_ARRAYS
	while ((nextIndex = ex->ct[i1].nextIndex_diff) != 0xff)
#else
	while ((nextIndex = ex->ct[i1].nextIndex) != -1)
#endif
	{
#ifdef USE_COMPRESSED_ARRAYS
		nextIndex = ex->ct[i1].nextIndex_diff;
		if (nextIndex != 254) nextIndex += i1;
		else nextIndex = keyval_find(&ex->ct_nextIndex_vector,i1);
#else
		nextIndex = ex->ct[i1].nextIndex;
#endif
		i2 = nextIndex;

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

	child_interval.l = -1;
	child_interval.r = -1;

	LEAVE("l=%d,r=%d",child_interval.l,child_interval.r);
	return child_interval;
}

int exSuffixArray8_subString_callback(struct exSuffixArray8 *ex, const char p[], int m, int (*hitcallback)(int pos, void *userdata), void *userdata)
{
	interval_t ij;
	int c;

	/** BEGIN_DEBUG **/ 
	ENTER("ex=%p,p=%p,m=%d,cb=%p,userdata=%p",ex,p,m,hitcallback,userdata);
	/** END_DEBUG **/ 

	c = 0;

	ij.l = 0;
	ij.r = ex->n - 1;
	ij = exSuffixArray8_getInterval(ex,ij,p[c],exSuffixArray8_getlcp(ex,ij));

	while (ij.l != -1 && c < m)
	{
		if (ij.l != ij.r)
		{
			int l = exSuffixArray8_getlcp(ex,ij);
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
			ij = exSuffixArray8_getInterval(ex,ij,p[c],exSuffixArray8_getlcp(ex,ij));
			if (ij.l == -1)
			{
				LEAVE("not found");
				return 0;
			}
		} else
		{
			int j;

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

			if (!hitcallback(ex->sa[j],userdata))
				break;
		}	
	}
	LEAVE("found %d occurrences",ij.r - ij.l + 1);
	return ij.r - ij.l + 1;
}

#if 0

static interval_t get_first_interval(const struct exSuffixArray8 *ex, interval_t ij)
{
	int i, j, i1, up;
	interval_t child_interval;

	i = ij.l;
	j = ij.r;

	up = exSuffixArray8_get_child_up(ex,j+1);

	if (i < up && up <= j) i1 = up;
	else i1 = exSuffixArray8_get_child_down(ex,i);
	
	child_interval.l = i;
	child_interval.r = i1 - 1;

	return child_interval;
}

static interval_t get_next_interval(const struct exSuffixArray8 *ex, interval_t last_child_interval, interval_t parent_interval)
{
	int nextIndex;
	int i1,i2;
	
	i1 = last_child_interval.r + 1;
	
#ifdef USE_COMPRESSED_ARRAYS
	if ((nextIndex = ex->ct[i1].nextIndex_diff) != 0xff)
#else
	if ((nextIndex = ex->ct[i1].nextIndex) != -1)
#endif
	{
#ifdef USE_COMPRESSED_ARRAYS
		nextIndex = ex->ct[i1].nextIndex_diff;
		if (nextIndex != 254) nextIndex += i1;
		else nextIndex = keyval_find(&ex->ct_nextIndex_vector,i1);
#else
		nextIndex = ex->ct[i1].nextIndex;
#endif
		i2 = nextIndex;

		last_child_interval.l = i1;
		last_child_interval.r = i2 - 1;
		assert(last_child_interval.l <= last_child_interval.r);
		assert(last_child_interval.l >= parent_interval.l && last_child_interval.r <= parent_interval.r);

//		printf("%d %d %c\n", child_interval.l, child_interval.r,ex->s[ex->sa[i1]+l]);

		return last_child_interval;
	}

	last_child_interval.l = i1;
	last_child_interval.r = parent_interval.r;

	assert(last_child_interval.l <= last_child_interval.r);
	assert(last_child_interval.l >= parent_interval.l && last_child_interval.r <= parent_interval.r);

//	printf("%d %d %c\n", child_interval.l, child_interval.r, ex->s[ex->sa[i1]+l]);
	return last_child_interval;
}


static int exSuffixArray8_enumerate_recursive(struct exSuffixArray8 *ex, interval_t interval)
{
	interval_t child_interval;

	int l = exSuffixArray8_getlcp(ex,interval);

	ENTER("ex=%p,lcp_interval.l=%d,lcp_interval.r=%d,c=%d,l=%d",ex,lcp_interval.l,lcp_interval.r,c,l);

	printf("ENTER left=%d right=%d lcp=%d\n",interval.l, interval.r, l);
	if (interval.l == interval.r)
	{
		printf("LEAVE\n");
		return 0;
	}

	child_interval = get_first_interval(ex, interval);
	printf("left=%d right=%d %c\n", child_interval.l, child_interval.r, ex->s[ex->sa[child_interval.l]+l]);
	exSuffixArray8_enumerate_recursive(ex,child_interval);

	while (child_interval.r != interval.r)
	{
		child_interval = get_next_interval(ex,child_interval,interval);
		printf("left=%d right=%d %c\n", child_interval.l, child_interval.r, ex->s[ex->sa[child_interval.l]+l]);
		exSuffixArray8_enumerate_recursive(ex,child_interval);
	}
	printf("LEAVE\n");
#if 0
	if (lcp_interval.l == lcp_interval.r) return 0;

	up = exSuffixArray8_get_child_up(ex,j+1);

	if (i < up && up <= j) i1 = up;
	else i1 = exSuffixArray8_get_child_down(ex,i);

	child_interval.l = i;
	child_interval.r = i1 - 1;

	printf("%d %d %c\n", child_interval.l, child_interval.r, ex->s[ex->sa[i]+l]);

	exSuffixArray8_enumerate_recursive(ex,child_interval);

#ifdef USE_COMPRESSED_ARRAYS
	while ((nextIndex = ex->ct[i1].nextIndex_diff) != 0xff)
#else
	while ((nextIndex = ex->ct[i1].nextIndex) != -1)
#endif
	{
#ifdef USE_COMPRESSED_ARRAYS
		nextIndex = ex->ct[i1].nextIndex_diff;
		if (nextIndex != 254) nextIndex += i1;
		else nextIndex = keyval_find(&ex->ct_nextIndex_vector,i1);
#else
		nextIndex = ex->ct[i1].nextIndex;
#endif
		i2 = nextIndex;

		/* Child interval is [i1..i2-1] */
//		if (ex->s[ex->sa[i1]+l] == c)
		{
			child_interval.l = i1;
			child_interval.r = i2 - 1;
			assert(child_interval.l <= child_interval.r);
			assert(child_interval.l >= i && child_interval.r <= j);

			printf("%d %d %c\n", child_interval.l, child_interval.r,ex->s[ex->sa[i1]+l]);
			exSuffixArray8_enumerate_recursive(ex,child_interval);
		}

		i1 = i2;
	}

//	if (ex->s[ex->sa[i1]+l] == c)
	{
		child_interval.l = i1;
		child_interval.r = j;
		assert(child_interval.l <= child_interval.r);
		assert(child_interval.l >= i && child_interval.r <= j);

		printf("%d %d %c\n", child_interval.l, child_interval.r, ex->s[ex->sa[i1]+l]);
		exSuffixArray8_enumerate_recursive(ex,child_interval);
	}
#endif
}

int exSuffixArray8_enumerate(struct exSuffixArray8 *ex, int k)
{
	int c;
	interval_t ij;

printf("enumerate %d %s\n",ex->n,ex->s);

	c = 0;

	ij.l = 0;
	ij.r = ex->n - 1;


	exSuffixArray8_enumerate_recursive(ex,ij);

	return 0;
}

#endif
