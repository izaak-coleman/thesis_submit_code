#ifndef KMER__HASH_GENHASH_H
#define KMER__HASH_GENHASH_H

struct gen_hash;

struct gen_hash *gen_hash_create(char **strings, int nstrings, int query_length, int options);
int gen_hash_lookup(struct gen_hash *ga, const char *p, int (*hit_callback)(int index, int pos, void *userdata), void *userdata);
void gen_hash_delete(struct gen_hash *gh);

/* The following flags can be used for options parameter of gen_hash_create() */

/* (default) index of strings is determined directly which
 * is fast but requires more memory */
#define GENHASH_OPT_DIRECT_INDEX_SEARCH (1<<0)

/* index of strings is determined via a binary search which
 * is slower but requires less memory */
#define GENHASH_OPT_BINARY_INDEX_SEARCH (1<<1)


#endif /* KMER__HASH_GENHASH_H */
