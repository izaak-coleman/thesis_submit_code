#ifndef __ARRAY_H
#define __ARRAY_H

struct exSuffixArray;

struct exSuffixArray *exSuffixArray_create(const int s[], int n, int K);
void exSuffixArray_delete(struct exSuffixArray *ex);
int exSuffixArray_subString(struct exSuffixArray *ex, const int p[], int m);
int exSuffixArray_subString_callback(struct exSuffixArray *ex, const int p[], int m, int (*hitcallback)(int pos, void *userdata), void *userdata);

#endif /*__ARRAY_H*/
