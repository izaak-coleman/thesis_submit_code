#ifndef KMER__GENARRAY_H
#define KMER__GENARRAY_H

struct GenArray;

#define GENSARRAY_OPT_DIRECT_INDEX_SEARCH (1<<0)

#define GENSARRAY_OPT_BINARY_INDEX_SEARCH (1<<1)

#define GENSARRAY_OPT_8BIT (1<<2)


struct GenArray *genArray_create(int **strings, int nstrings, int options);
struct GenArray *genArray_create_from_char(char **strings, int nstrings, int options);
int genArray_lookup(struct GenArray *ga, const int *p, int m, int (*hit_callback)(int index, int pos, void *userdata), void *userdata);
int genArray_lookup_from_char(struct GenArray *ga, const char *p, int m, int (*hitcallback)(int index, int pos, void *userdata), void *userdata);
void genArray_delete(struct GenArray *ga);



#endif /* KMER__GENARRAY_H */
