#ifndef __AMINO_H_
#define __AMINO_H_

extern unsigned char __amino[256];
extern unsigned char __ext_amino[256];
extern unsigned char __dna[256];

#define is_amino(c) (__amino[(unsigned char)c])
#define is_extamino(c) (__ext_amino[(unsigned char)c])
#define is_dna(c) (__dna[(unsigned char)(c)])

void init_amino(void);

#endif /*__AMINO_H_*/
