#ifndef __MYSTRING_H
#define __MYSTRING_H

#include <stddef.h>

size_t strlcpy(char *dst, const char *src, size_t size);
size_t strlcat(char *dst, const char *src, size_t size);



/* Basic string functions */
typedef struct
{
	char *str;
	int len;
	int allocated;
} string;

int string_initialize(string *string, unsigned int size);
int string_append(string *string, char *appstr);
int string_append_part(string *string, char *appstr, int bytes);
void string_crop(string *string, int startpos, int endpos);
void string_free(string *string);

#endif
