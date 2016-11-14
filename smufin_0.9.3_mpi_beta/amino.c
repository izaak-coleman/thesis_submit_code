/**
 * @file amino.c
 */

#include <ctype.h>

#include "amino.h"

static const char amino[] = {'A','C','D','E','F',
			     'G','H','I','K','L',
			     'M','N','P','Q','R',
			     'S','T','U','V','W','Y'};

static const char ext_amino[] = {'A','B', 'C','D','E','F',
			     'G','H','I','K','L',
			     'M','N','P','Q','R',
			     'S','T','U','V','W',
			     'X','Y'};
			     
static const char dna[] = { 'A','C','G', 'T' };

unsigned char __amino[256];
unsigned char __ext_amino[256];
unsigned char __dna[256];

void init_amino(void)
{
  int i;

  for (i=0;i<sizeof(amino)/sizeof(amino[0]);i++)
    {
      __amino[toupper(amino[i])] = 1;
      __amino[tolower(amino[i])] = 1;
    }
  for (i=0;i<sizeof(ext_amino)/sizeof(ext_amino[0]);i++)
    {
      __ext_amino[toupper(amino[i])] = 1;
      __ext_amino[tolower(amino[i])] = 1;
    }
    for (i=0;i<sizeof(dna)/sizeof(dna[0]);i++)
    {
    	__dna[toupper(dna[i])] = 1;	
    	__dna[tolower(dna[i])] = 1;	
    	
    }
}

int are_sequences_dna(const char **seq,int nseq)
{
	int i;
	char *p;
	for (i=0;i<nseq;i++)
	{
		p = seq[i];
		while (*p)
		{
			if (!is_dna(*p))
				return 0;
			p++;
		}
	}
	return 1;
}
