/*
 * genarray.c
 */


#include <assert.h>
#include <string.h>

#include "gmem.h"
#include "hash.h"
#include "trace.h"
#include "util.h"

#include "genhash_private.h"
#include "genhash.h"

/*************************************************************/


/*******************************************************
 Construct a generalised kmer hash index containing
 nstrings strings. Note that the lineendings of the
 strings in the array strings have to be a 0 "byte".
*********************************************************/
struct gen_hash *gen_hash_create(char **strings, int nstrings, int query_length, int options)
{
	struct gen_hash *gh;
	int i,len;
	char *ptr;
	int *idx_ptr;

	int index_type = IDXT_DIRECT;
	
	if (options & GENHASH_OPT_BINARY_INDEX_SEARCH) index_type = IDXT_BINARY;

	ENTER("strings=%p, nstrings=%d, query_length=%d",strings,nstrings,query_length);

	/* Determine size of the concatenated string */
	len = 0;
	for (i=0;i<nstrings;i++)
		len += strlen(strings[i]);

	LOG(LOG4C_PRIORITY_DEBUG,"%d strings with %d chars",nstrings,len);

	gh = gsuffix_malloc(sizeof(struct gen_hash));
	memset(gh,0,sizeof(*gh));

	/* Allocate mem for concatenated string including additional space
	 * for string delimiter, line ending */
	gh->concat = gsuffix_malloc(sizeof(gh->concat[0])*(len+nstrings+1));
	gh->n = nstrings;
	gh->index_type = index_type;

	/* Allocate mem for position array */
	gh->pos_array = gsuffix_malloc(sizeof(gh->pos_array[0])*(nstrings+1));
	
	/* If direct index type is requested, allocate the needed resources as well */
	if (index_type == IDXT_DIRECT) 
		gh->idx_array = gsuffix_malloc(sizeof(gh->idx_array[0])*(len+1+nstrings));
	
	/* Fill in the concatenated string */
	for (i=0, ptr=gh->concat, idx_ptr=gh->idx_array; i<nstrings; i++)
	{
		char *cur = strings[i];
		int c;

		gh->pos_array[i] = ptr - gh->concat;

		while ((c = *cur))
		{
			*ptr++ = c;
			if (idx_ptr) *idx_ptr++= i;
			cur++;
		}
		*ptr++ = 254; /* string delimiter */
		if (idx_ptr) *idx_ptr++= i;
	}
	gh->pos_array[i] = ptr - gh->concat; /* Don't forget the last position */

/*// checks if the interval is set up properly
	{
		int q;
		for (q=0;q<len+nstrings;q++)
		{
			int interval = find_interval(gh->pos_array,gh->n,q);
			int index = gh->idx_array[q];

			if (interval != index)
			{
				printf("q=%d interval=%d index=%d\n",q,interval,index);		
			}
		}
		
	}
*/

	/* Store the real length of the string */
	gh->concat_len = ptr - gh->concat;

	/* And now construct the hash from the concatenated string */
	if (!(gh->ht = kmer_hash_table_create((unsigned char*)gh->concat,gh->concat_len,-1,query_length)))
	{
		gen_hash_delete(gh);
		gh = NULL;
	}

	LEAVE("RC=%p",gh);
	return gh;
}

/*******************************************************
 Frees all resources taken by the hash
*********************************************************/
void gen_hash_delete(struct gen_hash *gh)
{
	if (gh->ht) kmer_hash_table_delete(gh->ht);
	gsuffix_free(gh->pos_array);
	gsuffix_free(gh->idx_array);
	gsuffix_free(gh->concat);
	gsuffix_free(gh);
}

/*******************************************************
 This is the callback for genArray_lookup. It determines
 the index of the string via a binary search
*********************************************************/
static int gen_hash_lookup_callback_binary(int pos, void *userdata)
{
	struct gen_hash *ga = (struct gen_hash*)userdata;
	int i;

	i = find_interval(ga->pos_array,ga->n,pos);

	ga->hits++;
	ga->hitcallback(i, pos - ga->pos_array[i], ga->userdata);
	return 1;	
}

/*******************************************************
 This is the callback for genArray_lookup. It determines
 the index of the string directly from the idx array.
*********************************************************/
static int gen_hash_lookup_callback_direct(int pos, void *userdata)
{
	struct gen_hash *ga = (struct gen_hash*)userdata;
	int i;

	i = ga->idx_array[pos];

	ga->hits++;
	ga->hitcallback(i, pos - ga->pos_array[i], ga->userdata);
	return 1;	
}

/*******************************************************
 Looks up a given p with length m. On every hit,
 hit_callback is called with arguments "index" and "pos"
 The index argument determines the string index as upon
 creation of the genarray and the argument "pos"
 determines the position of the hit within that string.

 Returns the number of hits.
*********************************************************/
int gen_hash_lookup(struct gen_hash *ga, const char *p, int (*hitcallback)(int index, int pos, void *userdata), void *userdata)
{
	int (*gen_hash_lookup_callback)(int pos, void *userdata);
	
	/** BEGIN_DEBUG **/
	ENTER("ga=%p,p=%p,cb=%p,userdata=%p",ga,p,hitcallback,userdata);
	/** END_DEBUG **/

	ga->hitcallback = hitcallback;
	ga->userdata = userdata;
	ga->hits = 0;

	/* Choose the proper callback function */
	if (ga->index_type == IDXT_DIRECT) gen_hash_lookup_callback = gen_hash_lookup_callback_direct;
	else gen_hash_lookup_callback = gen_hash_lookup_callback_binary;

	kmer_hash_table_subString_callback(ga->ht,p,gen_hash_lookup_callback,ga);

	LEAVE("hits=%d",ga->hits);
	return ga->hits;
}
