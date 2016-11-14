#ifndef SUFFIX__HASH_HASH_PRIVATE_H
#define SUFFIX__HASH_HASH_PRIVATE_H

struct kmer_hash_entry
{
	const char *str;
	struct kmer_hash_entry *next;
};

struct kmer_hash_table
{
	int query_length;	

	int table_size;
	struct kmer_hash_entry *table;
	const char *str;

	struct item_pool *entry_pool;
};

#endif
