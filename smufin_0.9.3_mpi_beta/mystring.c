/*
 * mystring.c
 */

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mystring.h"

size_t strlcpy(char *dst, const char *src, size_t siz)
{
	register char *d = dst;
	register const char *s = src;
	register size_t n = siz;

	if(n != 0 && --n != 0)
	{
		do
		{
			if(((*d++) = (*s++)) == '\0')
				break;
		}
		while(--n != 0);
	}

	if(n == 0)
	{
		if(siz != 0)
			(*d) = '\0'; /* NUL-terminate dst */

		while (*s++);
	}

	return (s - src - 1);	/* count does not include NUL */
}

size_t strlcat(char *dst, const char *src, size_t siz)
{
	register char *d = dst;
	register const char *s = src;
	register size_t n = siz;
	size_t result;
	size_t dlen;

	while(n-- != 0 && (*d) != '\0')
		d++;

	dlen = d - dst;
	n = siz - dlen;

	if (n == 0)
	{
		result = dlen + strlen(s);
	}
	else
	{
		while((*s) != '\0')
		{
			if(n != 1)
			{
				(*d++) = (*s);

				n--;
			}

			s++;
		}

		(*d) = '\0';

		result = dlen + (s - src); 
	}

	return result;
}

int string_initialize(string *string, unsigned int size)
{
	if (!size) size = 1;
	if (!(string->str = (char*)malloc(size)))
		return 0;

	string->str[0] = 0;
	string->allocated = size;
	string->len = 0;
	return 1;
}

int string_append_part(string *string, char *appstr, int bytes)
{
	int alloclen;

	if (!appstr || !bytes) return 1;

	alloclen = string->allocated;

	while (bytes + string->len >= alloclen) 
		alloclen *= 2;

	if (alloclen != string->allocated)
	{
		char *newstr;

		newstr = realloc(string->str,alloclen);
		if (!newstr) return 0;
		string->allocated = alloclen;
		string->str = newstr;
	}

	strncpy(&string->str[string->len],appstr,bytes);
	string->len += bytes;
	string->str[string->len] = 0;
	return 1;

}

int string_append(string *string, char *appstr)
{
	if (!appstr) return 1;
	return string_append_part(string,appstr,strlen(appstr));
}

void string_crop(string *string, int startpos, int endpos)
{
	if (startpos == 0)
	{
		string->len = endpos;
		string->str[endpos] = 0;
	} else
	{
		int newlen = endpos - startpos + 1;
		memmove(string->str, &string->str[startpos], newlen);
		string->len = newlen;
		string->str[newlen] = 0;
	}
}

void string_free(string *str)
{
	free(str->str);
}
