/*
 * alphabet.c
 */

#include <stdio.h>
#include <stdlib.h>

#include "alphabet.h"
#include "gmem.h"

static int alpha_type_;

static char map_[MAP_SZ];
static char unmap_[MAP_SZ];

static int alphabet_initialized = 0;

void init_alphabet(int alphabet_type)
{
  int i;

  if (alphabet_initialized) return;
 
  if (alphabet_type == DNA) alpha_type_ = DNA;
  else if (alphabet_type == RNA) alpha_type_ = RNA;
  else
  {
    fprintf(stderr,"Error: alphabet type of %d not supported.\n",alphabet_type);
    exit(1);
  }

  for (i=0; i < MAP_SZ; i++)
  {
    map_[i] = -1;
    unmap_[i] = -1;
  }

  switch (alpha_type_)
  {
  	case	DNA:
    	      	map_['a'] = map_['A'] = A_BASE;
    	   	map_['c'] = map_['C'] = C_BASE;
    		map_['g'] = map_['G'] = G_BASE;
    		map_['t'] = map_['T'] = T_BASE;
    		unmap_[A_BASE] = 'A';
    		unmap_[C_BASE] = 'C';
    		unmap_[G_BASE] = 'G';
    		unmap_[T_BASE] = 'T';
    		break;

  	case RNA:
    		map_['a'] = map_['A'] = A_BASE;
    		map_['c'] = map_['C'] = C_BASE;
    		map_['g'] = map_['G'] = G_BASE;
    		map_['u'] = map_['U'] = U_BASE;
    		unmap_[A_BASE] = 'A';
    		unmap_[C_BASE] = 'C';
    		unmap_[G_BASE] = 'G';
    		unmap_[U_BASE] = 'U';
    		break; 

 	default:
    		fprintf(stderr,"Error: invalid alphabet.\n");
    		exit(-1); 
  }
  alphabet_initialized = 1;
}


void free_alphabet(void) 
{
	
}

char mapchar(char ch)
{
  if (!alphabet_initialized)
  {
    fprintf(stderr,"Alphabet not initialized! Aborting\n");
    exit(-1);
  }

  return map_[(unsigned char)ch];
}



char unmapchar(char ch)
{
  if (!alphabet_initialized)
  {
    fprintf(stderr,"Alphabet not initialized! Aborting\n");
    exit(-1);
  }

  return unmap_[(unsigned char)ch];
}


void map_string(char *p, int len)
{
  int i;

  if (!alphabet_initialized)
  {
    fprintf(stderr,"Alphabet not initialized! Aborting\n");
    exit(-1);
  }

  for (i=0;i<len;++i)
  {
  	char s = mapchar(p[i]);
  	if (s == -1)
  	{
  		fprintf(stderr,"Unable to map character '%c'\n",p[i]);
  		exit(-1);
  	}
    p[i] = mapchar(p[i]);
  }
}

void unmap_string(char *p, int len)
{
  int i;

  if (!alphabet_initialized)
  {
    fprintf(stderr,"Alphabet not initialized! Aborting\n");
    exit(-1);
  }

  for (i=0;i<len;++i)
    p[i] = unmapchar(p[i]);
}

char *create_mapped_string(char *p, int len)
{
  int i;
  char *s = (char*)gsuffix_malloc(sizeof(char)*(len));	
  for (i=0;i<len;++i)
  {
    if ((s[i] = mapchar(p[i]))==-1)
    {
    	fprintf(stderr,"Unable to map character '%c'\n",p[i]);
    	exit(-1);
    }
  }
  return s;
}
	
