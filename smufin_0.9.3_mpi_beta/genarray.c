#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "array.h"
#include "gmem.h"
#include "trace.h"
#include "util.h"

#include "genarray_private.h"
#include "genarray.h"

#include "ds_ssort.h"

static int intstrlen(const int *str)
{
	const int *begin = str;

	while (*str) str++;

	return str - begin;
}

struct GenArray *genArray_create(int **strings, int nstrings, int options)
{
	struct GenArray *ga;
	int i,len;
	int *ptr,*idx_ptr;
	int index_type = IDXT_DIRECT;

	ENTER("strings=%p, nstrings=%d, options=%d",strings,nstrings,options);

	if (options & GENSARRAY_OPT_BINARY_INDEX_SEARCH) index_type = IDXT_BINARY;

	len = 0;
	for (i=0;i<nstrings;i++)
		len += intstrlen(strings[i]);

	LOG(LOG4C_PRIORITY_DEBUG,"%d strings with %d chars",nstrings,len);

	ga = gsuffix_malloc(sizeof(struct GenArray));
	memset(ga,0,sizeof(*ga));

	ga->concat = gsuffix_malloc(sizeof(ga->concat[0])*(len+4+nstrings));
	ga->n = nstrings;
	ga->index_type = index_type;

	ga->pos_array = gsuffix_malloc(sizeof(ga->pos_array[0])*(nstrings+1));
	if (index_type == IDXT_DIRECT) 
		ga->idx_array = gsuffix_malloc(sizeof(ga->idx_array[0])*(len+4+nstrings));

	for (i=0,ptr=ga->concat;i<nstrings;i++)
	{
		int *cur = strings[i];
		int c;

		ga->pos_array[i] = ptr - ga->concat;

		LOG(LOG4C_PRIORITY_DEBUG,"String %d starts at %d",i,ga->pos_array[i]);

		while ((c = *cur))
		{
			*ptr++ = c;
			if (idx_ptr) *idx_ptr++= i;
			cur++;
		}
		*ptr++ = 254; /* string delimiter */
		if (idx_ptr) *idx_ptr++= i;
	}

	ga->pos_array[i] = ptr - ga->concat; /* Hold last position */
	*ptr++ = 255; 	/* End the string with an high value */
	ptr[0] = 0; /* And then a nul "byte" */
	ptr[1] = 0; /* And then a nul "byte" */
	ptr[2] = 0; /* And then a nul "byte" */

	ga->concat_len = ptr - ga->concat;

	ga->ex = exSuffixArray_create(ga->concat,ga->concat_len,-1);

	LEAVE("RC=%p",ga);
	return ga;
}

struct GenArray *genArray_create_from_char(char **strings, int nstrings, int options)
{
	struct GenArray *ga;
	int i,len;
	int *idx_ptr;
	int index_type = IDXT_DIRECT;
	int eightbit = !!(options & GENSARRAY_OPT_8BIT);

	ENTER("strings=%p, nstrings=%d",strings,nstrings);


	if (options & GENSARRAY_OPT_BINARY_INDEX_SEARCH) index_type = IDXT_BINARY;

	len = 0;
	for (i=0;i<nstrings;i++)
		len += strlen(strings[i]);

	LOG(LOG4C_PRIORITY_DEBUG,"%d strings with %d chars",nstrings,len);

	ga = gsuffix_malloc(sizeof(struct GenArray));
	memset(ga,0,sizeof(*ga));

	ga->n = nstrings;
	ga->index_type = index_type;

	ga->pos_array = gsuffix_malloc(sizeof(ga->pos_array[0])*(nstrings+1));
	if (index_type == IDXT_DIRECT) 
		ga->idx_array = gsuffix_malloc(sizeof(ga->idx_array[0])*(len+4+nstrings));

	if (eightbit)
	{
		unsigned char *ptr;
		int overshoot = init_ds_ssort(500,2000);

		if (!overshoot)
			return NULL;

		ga->concat8 = gsuffix_malloc(sizeof(ga->concat8[0])*(len+4+nstrings+overshoot));

		for (i=0,ptr=ga->concat8,idx_ptr=ga->idx_array;i<nstrings;i++)
		{
			char *cur = strings[i];
			int c;

			ga->pos_array[i] = ptr - ga->concat8;

			LOG(LOG4C_PRIORITY_DEBUG,"String %d starts at %d",i,ga->pos_array[i]);

			while ((c = *cur))
			{
				*ptr++ = c;
				if (idx_ptr) *idx_ptr++= i;
				cur++;
			}
			*ptr++ = 254; /* string delimiter */
			if (idx_ptr) *idx_ptr++= i;
		}

		ga->pos_array[i] = ptr - ga->concat8; /* Hold last position */
		*ptr++ = 255; 	/* End the string with an high value */
		ptr[0] = 0; /* And then a nul "byte" */
		ptr[1] = 0; /* And then a nul "byte" */
		ptr[2] = 0; /* And then a nul "byte" */

		ga->concat_len = ptr - ga->concat8;

		ga->ex8 = exSuffixArray8_create(ga->concat8,ga->concat_len);

	} else
	{
		int *ptr;

		ga->concat = gsuffix_malloc(sizeof(ga->concat[0])*(len+4+nstrings));

		for (i=0,ptr=ga->concat,idx_ptr=ga->idx_array;i<nstrings;i++)
		{
			char *cur = strings[i];
			int c;

			ga->pos_array[i] = ptr - ga->concat;

			LOG(LOG4C_PRIORITY_DEBUG,"String %d starts at %d",i,ga->pos_array[i]);

			while ((c = *cur))
			{
				*ptr++ = c;
				if (idx_ptr) *idx_ptr++= i;
				cur++;
			}
			*ptr++ = 254; /* string delimiter */
			if (idx_ptr) *idx_ptr++= i;
		}

		ga->pos_array[i] = ptr - ga->concat; /* Hold last position */
		*ptr++ = 255; 	/* End the string with an high value */
		ptr[0] = 0; /* And then a nul "byte" */
		ptr[1] = 0; /* And then a nul "byte" */
		ptr[2] = 0; /* And then a nul "byte" */

		ga->concat_len = ptr - ga->concat;

		ga->ex = exSuffixArray_create(ga->concat,ga->concat_len,-1);
	}
	LEAVE("RC=%p",ga);
	return ga;
}

static int genArray_lookup_callback_binary(int pos, void *userdata)
{
	struct GenArray *ga = (struct GenArray*)userdata;
	int i;
	
	i = find_interval(ga->pos_array,ga->n,pos);
	assert(i <= ga->n);

	ga->hits++;
	ga->hitcallback(i, pos - ga->pos_array[i], ga->userdata);
	return 1;	
}

static int genArray_lookup_callback_direct(int pos, void *userdata)
{
	struct GenArray *ga = (struct GenArray*)userdata;
	int i;

	i = ga->idx_array[pos];

	ga->hits++;
	ga->hitcallback(i, pos - ga->pos_array[i], ga->userdata);
	return 1;	
}

int genArray_lookup(struct GenArray *ga, const int *p, int m, int (*hitcallback)(int index, int pos, void *userdata), void *userdata)
{
	int (*genArray_lookup_callback)(int pos, void *userdata);

	assert(ga->ex != NULL);

	/** BEGIN_DEBUG **/
	ENTER("ga=%p,p=%p,m=%d,cb=%p,userdata=%p",ga,p,m,hitcallback,userdata);
	SHOWINTSTRING(LOG4C_PRIORITY_TRACE,p,0,m);
	/** END_DEBUG **/

	if (ga->index_type == IDXT_DIRECT) genArray_lookup_callback = genArray_lookup_callback_direct;
	else genArray_lookup_callback = genArray_lookup_callback_binary;

	ga->hitcallback = hitcallback;
	ga->userdata = userdata;
	ga->hits = 0;
	exSuffixArray_subString_callback(ga->ex,p,m,genArray_lookup_callback,ga);

	LEAVE("hits=%d",ga->hits);
	return ga->hits;
}

int genArray_lookup_from_char(struct GenArray *ga, const char *p, int m, int (*hitcallback)(int index, int pos, void *userdata), void *userdata)
{
	int (*genArray_lookup_callback)(int pos, void *userdata);

	if (ga->index_type == IDXT_DIRECT) genArray_lookup_callback = genArray_lookup_callback_direct;
	else genArray_lookup_callback = genArray_lookup_callback_binary;

	ga->hitcallback = hitcallback;
	ga->userdata = userdata;
	ga->hits = 0;

	if (ga->ex)
	{
		int ipat[m];
		int i;

		for (i=0;i<m;i++)
			ipat[i] = (unsigned char)p[i];

		exSuffixArray_subString_callback(ga->ex,ipat,m,genArray_lookup_callback,ga);
	} else
	{
		exSuffixArray8_subString_callback(ga->ex8,p,m,genArray_lookup_callback,ga);
	}

	LEAVE("hits=%d",ga->hits);
	return ga->hits;
}

void genArray_delete(struct GenArray *ga)
{
	exSuffixArray_delete(ga->ex);
	exSuffixArray8_delete(ga->ex8);

	gsuffix_free(ga->concat);
	gsuffix_free(ga->concat8);
	gsuffix_free(ga->pos_array);
	gsuffix_free(ga->idx_array);
	gsuffix_free(ga);
}

