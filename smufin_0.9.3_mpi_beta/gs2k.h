#ifndef GS_2_K__
#define GS_2_K__

#include <stdio.h>

extern unsigned char __amino[256];

extern unsigned char __ext_amino[256];

extern unsigned char __dna[256];

#define is_amino(c) (__amino[(unsigned char)c])

#define is_extamino(c) (__ext_amino[(unsigned char)c])

#define is_dna(c) (__dna[(unsigned char)(c)])

void init_amino(void);

int are_sequences_dna(const char **seq,int nseq);

struct plain_sequence
{
	int len;    /**< @brief length of the sequence in characters */
	char *name; /**< @brief id of the sequence */
	char *seq;  /**< @brief the sequence */
};

struct plain_fasta
{
	int num; /**< @brief number of sequences */
	int allocated; /**< @brief number of entries that can be placed in the seqs array */
	struct plain_sequence **seqs; /**< @brief the sequences */
};

struct labeled_sequence
{
	int len;		/**< @brief length of the sequence in characters (nt or amino acids) */
	char *name;		/**< @brief sequence's name */
	char *seq;		/**< @brief sequence */
	double label;	/**< @brief the associated label */
};


struct labeled_fasta
{
	int num; /**< @brief number of sequences */
	int allocated; /**< @brief number of entries that can be placed in the seqs array */
	struct labeled_sequence **seqs; /**< @brief the sequences */
};

int fasta_read(FILE *fh, int (*sequence_function)(char *name, char *sequence, int sequence_len, void *userdata), void *userdata);

struct plain_fasta *fasta_read_plain(FILE *fh);
void fasta_free_plain(struct plain_fasta *fasta);

struct labeled_fasta *fasta_read_labeled(FILE *fh);
void fasta_free_labeled(struct labeled_fasta *fasta);

#define GENSARRAY_OPT_DIRECT_INDEX_SEARCH (1<<0)

#define GENSARRAY_OPT_BINARY_INDEX_SEARCH (1<<1)

#define GENSARRAY_OPT_8BIT (1<<2)

#define GENHASH_OPT_DIRECT_INDEX_SEARCH (1<<0)

#define GENHASH_OPT_BINARY_INDEX_SEARCH (1<<1)

enum gsuffix_type
{
    GSUFFIX_SUFFIXARRAY,      //!< A generalized suffix array.
    GSUFFIX_TRUNCSUFFIXTREE,//!< A truncated suffix tree.
//    GSUFFIX_SHTRUNCSUFFIXTREE,//!< GSUFFIX_SHTRUNCSUFFIXTREE
    GSUFFIX_HASH,             //!< A basic hash index.
    GSUFFIX_SUFFIXTREE,     //!< A suffix tree.
//    GSUFFIX_SHSUFFIXTREE,     //!< GSUFFIX_SHSUFFIXTREE
    GSUFFIX_DIRECT            //!< A direct approach (similar to hash)
};


#define GSUFFIX_LOTRUNCSUFFIXTREE GSUFFIX_TRUNCSUFFIXTREE

#define GSUFFIX_LOSUFFIXTREE GSUFFIX_SUFFIXTREE


struct gsuffix { };

struct gsuffix *gsuffix_create(char **strings, int nstrings, enum gsuffix_type type, int query_length, int options);
void gsuffix_delete(struct  gsuffix *gs);

/*********************************************************************
 *  There are three kinds of lookup function.                        *
 *  1) gsuffix_lookup is used to look up an individual query sequence*
 *  in one of the datastructures. gs is the datastructure created as *
 *  above, p is the pattern to be sought, m is the length of p , hit *
 *  callback is a callback function that expects to receive the      *
 *  arguments index for the index of the sequence in which the       *
 *  pattern is found, and pos for the position of the pattern in the *
 *  sequence, and userdata, a user-defined structure in which to     *
 *  store the results. Note that this callback function will be      *
 * called once for each match for p in all of the input sequences.   *
 *  2) gsuffix_lookup_exists looks to see if the sequence 'p' exists *
 *     in the input data represented by the suffix datastructure.    *
 *     After the first match, the function returns. The callback     *
 *     function does not get information about the sequence and      *
 *     position where 'p' was found.                                 *
 *  3) gsuffix_enum_kmers. Finds each ocurrence of every substring   *
 *       of length k in the input sequences. The callback function   *
 *      is called once for each hit  and gets information about the  *
 * sequence and position within the sequence of the hit.             *
 *********************************************************************/
int gsuffix_lookup(struct  gsuffix *gs,
		   const char *p,
		   int m,
		   int (*hit_callback)(int index, int pos, void *userdata),
		   void *userdata);
int gsuffix_lookup_exists(struct gsuffix *gs, const char *p, int m);
int gsuffix_enum_k_mers(struct gsuffix *gs, int k,
			int (*callback)(char *kmer, int id, int pos, void *userdata),
			void *userdata);

#endif
