/*
 * hash.c
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "gmem.h"

#include "hash_private.h"
#include "hash.h"

static unsigned long sdbm(const unsigned char *str, int len)
{
	unsigned long hash = 0;
	unsigned int c;

  	while (len && (c = *str++))
  	{
		hash = c + (hash << 6) + (hash << 16) - hash;
		len--;
  	}

	return hash;
}

struct kmer_hash_table *kmer_hash_table_create(const unsigned char s[],
				int n, int K, int query_length)
{
	struct kmer_hash_table *ht;
	const unsigned char *end;
	int table_size;

	if (n < query_length) return NULL;

	ht = gsuffix_malloc(sizeof(*ht));
	memset(ht,0,sizeof(*ht));

	table_size = 8192*2;

	ht->entry_pool = item_pool_create(sizeof(struct kmer_hash_entry));
	ht->query_length = query_length;
	ht->table_size = table_size;
	ht->table = gsuffix_malloc(ht->table_size * sizeof(ht->table[0]));
	memset(ht->table,0,ht->table_size * sizeof(ht->table[0]));
	ht->str = (char*)s;

	for (end = s + n - query_length + 1;s<end;s++)
	{
		unsigned long key = sdbm(s,query_length);

		struct kmer_hash_entry *entry = &ht->table[key % table_size];
		if (!entry->str)
			entry->str = (char*)s;
		else
		{
			struct kmer_hash_entry *new_entry;

			new_entry = item_pool_alloc(ht->entry_pool);
			new_entry->str = (char*)s;
			new_entry->next = entry->next;
			entry->next = new_entry;
		}
	}
	return ht;
}

void kmer_hash_table_delete(struct kmer_hash_table *ht)
{
	item_pool_delete(ht->entry_pool);
	gsuffix_free(ht->table);
	gsuffix_free(ht);
}

int kmer_hash_table_subString_callback(struct kmer_hash_table *ht,
		const char p[],
		int (*hitcallback)(int pos, void *userdata),
		void *userdata)
{
	struct kmer_hash_entry *entry;
	unsigned long val;
	int hits;

	val = sdbm((const unsigned char*)p,ht->query_length);
	entry = &ht->table[val % ht->table_size];
	hits = 0;

	while (entry && entry->str)
	{
		if (!strncmp(p,entry->str,ht->query_length))
		{
			hits++;
			if (!(hitcallback(entry->str - ht->str,userdata)))
				break;
		}
		entry = entry->next;
	}

	return hits;
}
