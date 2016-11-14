#include <ctype.h>
#include <stdlib.h>
#include <string.h>

struct alphabet_internal
{
	unsigned int map_char[256];
	unsigned char *alphabet_string;
	int alphabet_length;
};

struct alphabet *gsuffix_create_alphabet(const char *alphabet)
{
	struct alphabet_internal *ai;
	int i;
	int length;

	if (!(ai = (struct alphabet_internal*)malloc(sizeof(*ai))))
		return NULL;

	length = strlen(alphabet);
	if (!(ai->alphabet_string = (unsigned char*)malloc(sizeof(unsigned char)*length)))
	{
		free(ai);
		return NULL;
	}

	for (i=0; i < 256; i++)
		ai->map_char[i] = -1;

	for (i=0;i<length;++i)
	{
		ai->alphabet_string[i] = alphabet[i];
		ai->map_char[(alphabet[i])] = i;
	}

    ai->alphabet_length = length;

	return (struct alphabet*)ai;
}

void gsuffix_free_alphabet(struct alphabet *a)
{
	struct alphabet_internal *ai;
	if (!(ai = (struct alphabet_internal*)a)) return;
	free(ai->alphabet_string);
	free(ai);
}

int gsuffix_next_sequence(struct alphabet *a, unsigned char *sequence, int length)
{
	struct alphabet_internal *ai;
	int pos,done;

	if (!(ai = (struct alphabet_internal*)a)) return 1;
	if (length == 0) return 1;

	if (sequence[0] == 0)
	{
		int i;

		for (i=0;i<length;i++)
			sequence[i] = ai->alphabet_string[i];
	}

	pos = 0;
	done = 0;

	while (pos < length)
	{
		if (ai->map_char[(sequence[pos])] < ai->alphabet_length - 1)
		{
			sequence[pos] = ai->alphabet_string[(ai->map_char[(sequence[pos])])+1];
			return 0;
		} else
		{
			sequence[pos] = ai->alphabet_string[0];
			++pos;
		}
	}
	/* we only exit the while loop here, if we cannot add our sequence */
	sequence[0] = 0;
	return 1;
}

