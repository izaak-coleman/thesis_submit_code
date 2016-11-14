#ifndef KMER__DIRECT_DIRECT_H
#define KMER__DIRECT_DIRECT_H

struct gen_direct;

struct gen_direct *gen_direct_create(char **strings, int nstrings, int query_length, int options);
int gen_direct_lookup(struct gen_direct *gd, const char *p, int (*hit_callback)(int index, int pos, void *userdata), void *userdata);
void gen_direct_delete(struct gen_direct *gd);


#endif /* KMER__DIRECT_DIRECT_H */
