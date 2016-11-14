#ifndef ALPHABET_H_
#define ALPHABET_H_


#define DNA 0
#define RNA 1
#define DNA_SZ  4
#define MAP_SZ  128
#define A_BASE  0
#define C_BASE  1
#define G_BASE 2
#define T_BASE 3
#define U_BASE 3

extern char mapchar(char ch);

extern char unmapchar(char ch);

extern void init_alphabet(int alphabet_type);

extern void free_alphabet(void);

extern void map_string(char *p, int len);

extern void unmap_string(char *p, int len);


extern char *create_mapped_string(char *p, int len);
	



#endif /*ALPHABET_H_*/
