#ifndef SUFFIX__HASH_HASH_H
#define SUFFIX__HASH_HASH_H

struct kmer_hash_table;

struct kmer_hash_table *kmer_hash_table_create(const unsigned char s[], int n, int K, int query_length);
void kmer_hash_table_delete(struct kmer_hash_table *ht);

int kmer_hash_table_subString_callback(struct kmer_hash_table *ht,
		const char p[],
		int (*hitcallback)(int pos, void *userdata),
		void *userdata);

#endif
