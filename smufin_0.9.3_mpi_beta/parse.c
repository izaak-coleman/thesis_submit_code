#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "parse.h"

char *skip_spaces(char *line)
{
	unsigned char c;

	while ((c = *line))
	{
		if (!isspace(c))
			break;
		line++;
	}

	return line;
}

char *skip_simple_string(char *line)
{
	char *ptr;
	unsigned char c;

	ptr = line;

	while ((c = *ptr))
	{
		if (isspace(c))
			break;
		ptr++;
	}

	return ptr;
}

char *parse_simple_string(char *line, char **dest_ptr)
{
	char *ptr;
	char *dest;
	int len;
	unsigned char c;

	ptr = line;

	while ((c = *ptr))
	{
		if (isspace(c))
			break;
		ptr++;
	}

	if (dest_ptr)
	{
		len = ptr - line;
		if (!(dest = (char*)malloc(len + 1)))
			return NULL;
		strncpy(dest, line, len);
		dest[len] = 0;
		*dest_ptr = dest;
	}

	return ptr;
}

char *parse_double(char *line, double *val_ptr)
{
	double val;
	char *end;

	val = strtod(line,&end);
	if (line == end) return NULL;

	*val_ptr = val;
	return end;
}
