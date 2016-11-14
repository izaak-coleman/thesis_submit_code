#ifndef SEQGEN_H_
#define SEQGEN_H_

struct alphabet;

struct alphabet *gsuffix_create_alphabet(const char *alphabet);
void gsuffix_free_alphabet(struct alphabet *a);
int gsuffix_next_sequence(struct alphabet *a, unsigned char *sequence, int length);

#endif /*SEQGEN_H_*/


