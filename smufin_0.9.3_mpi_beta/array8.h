#ifndef __ARRAY8BIT_H
#define __ARRAY8BIT_H

struct exSuffixArray8;

struct exSuffixArray8 *exSuffixArray8_create(unsigned char s[], int n);
void exSuffixArray8_delete(struct exSuffixArray8 *ex);
int exSuffixArray8_subString_callback(struct exSuffixArray8 *ex, const char p[], int m, int (*hitcallback)(int pos, void *userdata), void *userdata);

#endif /*__ARRAY8BIT_H*/
