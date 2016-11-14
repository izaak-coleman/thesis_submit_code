#include <string.h>

#include "direct.h"
#include "gmem.h"

struct entry
{
	int str_index;	
	int pos;
	struct entry *next;
};

struct gen_direct
{
	int hits;
	int query_length;

	int entry_table_length; 
	struct entry **entry_table;

	struct item_pool *entry_pool;	
};

#define convert(c)	\
switch (c) \
{ \
	case	'A': c = 0; break; \
	case	'C': c = 1; break; \
	case	'G': c = 2; break; \
	case	'T': c = 3; break; \
	default: \
		goto err; \
}

struct gen_direct *gen_direct_create(char **strings, int nstrings, int query_length, int options)
{
	int i;
	unsigned int mask;

	struct gen_direct *gd;

	if (query_length > 14) return NULL;

	gd = gsuffix_malloc(sizeof(*gd));
	memset(gd,0,sizeof(*gd));

	gd->query_length = query_length;
	gd->entry_table_length = 1;
	gd->entry_pool = item_pool_create(sizeof(struct entry));
	for (i=0;i<query_length;i++)
		gd->entry_table_length *= 4;

	gd->entry_table = gsuffix_malloc(gd->entry_table_length*sizeof(gd->entry_table[0]));
	memset(gd->entry_table,0,gd->entry_table_length*sizeof(gd->entry_table[0]));

	mask = gd->entry_table_length - 1;

	for (i=0;i<nstrings;i++)
	{
		int j;

		char *string = strings[i];
		char *ptr = string;
		unsigned char c;
		unsigned int val = 0;

		for (j=0;j<query_length-1;j++)
		{
			c = *ptr;
			if (!c) break;

			convert(c);
			val = (val << 2) | c;
			ptr++;
		}

		while ((c = *ptr++))
		{
			struct entry *old_entry, *new_entry;

			convert(c);
			val = (val << 2) | c;
			val &= mask;

			old_entry = gd->entry_table[val];
			new_entry = item_pool_alloc(gd->entry_pool);
			new_entry->str_index = i+1;
			new_entry->pos = ptr - query_length - string;
			new_entry->next = old_entry;
			gd->entry_table[val] = new_entry;
		}
	}

	return gd;

err:
	gen_direct_delete(gd);
	return NULL;
}

void gen_direct_delete(struct gen_direct *gd)
{
	item_pool_delete(gd->entry_pool);
	gsuffix_free(gd->entry_table);
	gsuffix_free(gd);
}

int gen_direct_lookup(struct gen_direct *gd, const char *p, int (*hit_callback)(int index, int pos, void *userdata), void *userdata)
{
	struct entry *entry;
	int i;
	unsigned int val = 0;

	for (i=0;i<gd->query_length;i++)
	{
		unsigned char c = *p++;
		convert(c);
		val = (val << 2) | c;
	}

	entry = gd->entry_table[val];

	while (entry)
	{
		gd->hits++;
		hit_callback(entry->str_index - 1, entry->pos, userdata);
		entry = entry->next;
	}

	return gd->hits;

err:
	return 0;
}
