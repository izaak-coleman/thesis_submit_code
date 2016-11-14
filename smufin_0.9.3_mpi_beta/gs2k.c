#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <ctype.h> /* isspace */

#include "mystring.h"
//#include "mem.h"
#include "parse.h"

#include "direct.h"
#include "genarray.h"
#include "genhash.h"
#include "stree.h"
#include "gmem.h"

#include "gs2k.h"

void* gsuffix_routine_thread(void* param);
int gsuffix_fi_asm_impl(char* mem, size_t size_mem, FILE* stream_out);
int gsuffix_co_asm_impl(char* mem, size_t size_mem, FILE* out);
int gsuffix_so_asm_impl(char* mem, size_t size_mem, FILE* stream_out);
int gsuffix_somatic_bp_asm_impl(char* ref_file, FILE* normal_mem, FILE* cancer_mem, FILE* out_file);

#define LINE_LENGTH 60000

int fasta_read(FILE *fh,
	       int (*sequence_function)(char *name,
					char *sequence,
					int sequence_len,
					void *userdata),
	       void *userdata)
{
	char *line;
	char name[1024];
	int line_count;

	if (!(line = (char*)malloc(LINE_LENGTH)))
		return 0;

	string cur_seq_str;

	line_count = 0;
	string_initialize(&cur_seq_str,10000);

	name[0] = 0;

	while (1)
	{
		int eof = !fgets(line,LINE_LENGTH,fh);
		line_count++;
		if (line[0] == '>' || eof)
		{
			if (cur_seq_str.len)
			{
				if (!(sequence_function(name,cur_seq_str.str,cur_seq_str.len,userdata)))
					goto bailout;

				string_crop(&cur_seq_str,0,0);
			}
			if (eof) break;

			if (line[0] == '>')
				strlcpy(name,line+1,sizeof(name));
		} else
		{
			if (line[0] != ';' && line[0] != '#')
			{
				static char buf[sizeof(line)];

				char c;
				char *ptr = line;
				char *dest_ptr = buf;

				while ((c = *ptr++))
				{
					if (c == '\r' || c == '\n')
					{
						*dest_ptr++ = 0;
						break;
					}
					*dest_ptr++ = c;
				}
				*dest_ptr++ = 0;
				string_append(&cur_seq_str,buf);
			}
		}
	}

	free(line);
	string_free(&cur_seq_str);
	return 1;

bailout:
	free(line);
	string_free(&cur_seq_str);
	return 0;
}


static int labeled_fasta_callback(char *name, char *sequence, int sequence_len, void *userdata)
{
	struct labeled_fasta *fast; /* int num;truct sequence **seqs; */
	struct labeled_sequence *lseq;
	char *ptr, *tmp;
	double current_label;

	fast = (struct labeled_fasta*) userdata;

	if (!(lseq = malloc(sizeof(struct labeled_sequence))))
	{
		fprintf(stderr,"Could not allocate memory at %s (%d)\n",__FILE__,__LINE__);
		return 0;
	}

	if (fast->num == fast->allocated)
	{
		fast->allocated *= 2;
		if (!(fast->seqs = realloc(fast->seqs,sizeof(struct sequence*)*fast->allocated)))
		{
			fprintf(stderr,"Could not allocate memory at %s (%d)\n",__FILE__,__LINE__);
			return 0;
		}
	}
	lseq->len = sequence_len;
	if (!(lseq->seq = (char*)malloc(sizeof(sequence[0])*(sequence_len+1))))
	{
		fprintf(stderr,"Could not allocate memory at %s (%d)\n",__FILE__,__LINE__);
		return 0;
	}
	tmp = sequence;
	ptr = lseq->seq;
	while (*tmp)
	{
	 	*ptr = toupper(*tmp);
	 	ptr++;tmp++;
	}
	*ptr = 0;

	ptr = name;
	ptr = parse_simple_string(ptr,&tmp);	/* read name */

	lseq->name = tmp;
	ptr = skip_spaces(ptr);			/* skip spaces after  ">name" to get the the label. */
	if (!(ptr = parse_double(ptr,&current_label)))  /* read label value */
	{
		fprintf(stderr,"Could not parse the label in line: >%s\n",name);
		exit(1);
	}

	lseq->label = current_label;
	fast->seqs[fast->num] = lseq;
	fast->num++;
	return 1;
}

struct labeled_fasta *fasta_read_labeled(FILE *fh)
{
	struct labeled_fasta *fasta;

	if (!(fasta = (struct labeled_fasta *)malloc(sizeof(struct labeled_fasta))))
		return NULL;
	memset(fasta,0,sizeof(fasta));

	fasta->num = 0;
	fasta->allocated = 128;
	if (!(fasta->seqs = malloc(sizeof(struct labeled_sequence*) * fasta->allocated)))
	{
		free(fasta);
		return NULL;
	}
	fasta_read(fh, &labeled_fasta_callback, fasta);
	return fasta;
}

void fasta_free_labeled(struct labeled_fasta *fasta)
{
	int i;

	for (i=0;i<fasta->num;i++)
	{
		free(fasta->seqs[i]->name);
		free(fasta->seqs[i]->seq);
		free(fasta->seqs[i]);
	}

	free(fasta->seqs);
	free(fasta);
}

static int plain_fasta_callback(char *name, char *sequence, int sequence_len, void *userdata)
{
	struct plain_fasta *fast; /* int num;truct sequence **seqs; */
	struct plain_sequence *seq;

	char *ptr,*tmp;
	int len;

	fast = (struct plain_fasta*) userdata;
	if (!(seq = malloc(sizeof(struct plain_sequence))))
	{
		fprintf(stderr,"Could not allocate memory at %s (%d)\n",__FILE__,__LINE__);
		return 0;
	}

	if (fast->num == fast->allocated)
	{
		fast->allocated *= 2;
		if (!(fast->seqs = realloc(fast->seqs,sizeof(struct plain_sequence*)*fast->allocated)))
		{
			fprintf(stderr,"Could not allocate memory at %s (%d)\n",__FILE__,__LINE__);
			return 0;
		}
	}
	seq->len = sequence_len;
	if (!(seq->seq = malloc(sizeof(sequence[0])*(sequence_len+1))))
	{
		fprintf(stderr,"Could not allocate memory at %s (%d)\n",__FILE__,__LINE__);
		return 0;
	}

	tmp = sequence;
	ptr = seq->seq;
	while (*tmp)
	{
	 	*ptr = toupper((unsigned char)*tmp);
	 	ptr++;tmp++;
	}
	*ptr = 0;
	seq->seq[sequence_len] = 0;

	ptr = name;
	len = 0;

	while (!isspace(*ptr) && len < strlen(name) ) {++ptr; ++len; }
	if (!(seq->name = malloc(len+1)))
	{
		fprintf(stderr,"Could not allocate memory at %s (%d)\n",__FILE__,__LINE__);
		return 0;
	}
	strncpy(seq->name,name,len);
	seq->name[len] = 0;

	fast->seqs[fast->num] = seq;
	fast->num++;
	return 1;
}

struct plain_fasta *fasta_read_plain(FILE *fh)
{
	struct plain_fasta *fasta;

	if (!(fasta = (struct plain_fasta *)malloc(sizeof(struct plain_fasta))))
		return NULL;
	memset(fasta,0,sizeof(struct plain_fasta));

	fasta->num = 0;
	fasta->allocated = 128;
	if (!(fasta->seqs = (struct labeled_sequence*)malloc(sizeof(struct labeled_sequence*) * fasta->allocated)))
	{
		free(fasta);
		return NULL;
	}
	fasta_read(fh, &plain_fasta_callback, fasta);
	return fasta;
}

void fasta_free_plain(struct plain_fasta *fasta)
{
	int i;
	for (i=0;i<fasta->num;++i)
	{
		free(fasta->seqs[i]->name);
		free(fasta->seqs[i]->seq);
		free(fasta->seqs[i]);
	}
	free(fasta->seqs);
	free(fasta);
}

struct gsuffix_internal
{
  enum  gsuffix_type type; /**< type of data structure */
  union
  {
    struct GenArray *genarray;
    struct gen_hash *genhash;
    struct gen_direct *gendirect;
    struct stree *lotruncsuffixtree;
    struct stree *losuffixtree;

  } ds; /**< data structure */
};

struct gsuffix *gsuffix_create(char **strings, int nstrings, enum gsuffix_type type, int query_length, int options)
{
	struct gsuffix_internal *gs;

	if (query_length == 0) {
		if ( type== GSUFFIX_LOTRUNCSUFFIXTREE) {
			fprintf(stderr, "Error: query_length option must be >0 for truncated suffix tree!\n");
			return NULL;
		} else if (type == GSUFFIX_DIRECT) {
			fprintf(stderr, "Error: query_length option must be >0 for direct search!\n");
			return NULL;
		} else if (type == GSUFFIX_HASH) {
			fprintf(stderr, "Error: query_length option must be >0 for hash search!\n");
			return NULL;
		}
	}

	if (!(gs = malloc(sizeof(*gs))))
		return NULL;

	gs->type = type;
	switch (type)
	{
		case	GSUFFIX_SUFFIXARRAY:
			if (!(gs->ds.genarray = genArray_create_from_char(strings, nstrings, options)))
				goto bailout;
			break;

		case	GSUFFIX_HASH:
			if (!(gs->ds.genhash = gen_hash_create(strings, nstrings, query_length, options)))
				goto bailout;
			break;

		case    GSUFFIX_LOTRUNCSUFFIXTREE:
		  if (!(gs->ds.lotruncsuffixtree = lokTstree_create_from_char(strings,nstrings,query_length)))
				goto bailout;
			break;
/*
 		case    GSUFFIX_SHTRUNCSUFFIXTREE:
			if (!(gs->ds.shtruncsuffixtree = shkTstree_create_from_char(strings,nstrings,query_length)))
				goto bailout;
			break;
*/
		case	GSUFFIX_LOSUFFIXTREE:
			if (!(gs->ds.losuffixtree = lostree_create_from_char(strings, nstrings)))
				goto bailout;
			break;
/*
		case	GSUFFIX_SHSUFFIXTREE:
			if (!(gs->ds.shsuffixtree = shstree_create_from_char(strings, nstrings, query_length)))
				goto bailout;
			break;
*/

		case	GSUFFIX_DIRECT:
			if (!(gs->ds.gendirect = gen_direct_create(strings, nstrings, query_length, 0)))
				goto bailout;
			break;

		default:
			goto bailout;
	}

	return (struct gsuffix*)gs;

bailout:
	free(gs);
	return NULL;
}

void gsuffix_delete(struct gsuffix *gs)
{
	struct gsuffix_internal *gsi = (struct gsuffix_internal*)gs;

	if (!gsi) return;

	switch (gsi->type)
	{
		case	GSUFFIX_SUFFIXARRAY:
			genArray_delete(gsi->ds.genarray);
			break;

		case	GSUFFIX_LOTRUNCSUFFIXTREE:
			lokTstree_delete_tree(gsi->ds.lotruncsuffixtree);
			break;
		case	GSUFFIX_LOSUFFIXTREE:
					lostree_delete_tree(gsi->ds.losuffixtree);
					break;

/*
		case	GSUFFIX_SHTRUNCSUFFIXTREE:
			shkTstree_delete_tree(kmers->ds.shtruncsuffixtree);
			break;
*/


/*
		case	GSUFFIX_SHSUFFIXTREE:
			shstree_delete_tree(kmers->ds.shsuffixtree);
			break;
*/
		case	GSUFFIX_HASH:
			gen_hash_delete(gsi->ds.genhash);
			break;

		case	GSUFFIX_DIRECT:
			gen_direct_delete(gsi->ds.gendirect);
			break;

		default:
			break;
	}
	gsuffix_free(gsi);
}

int gsuffix_load_fastq_sequences(char* mem_fastq_1, size_t size_mem_fastq_1, char* mem_fastq_2, size_t size_mem_fastq_2, FILE* fout)
{
	return my_main_mem_paired(mem_fastq_1, size_mem_fastq_1, mem_fastq_2, size_mem_fastq_2, fout);
}

void* gsuffix_routine(void* param)
{
	return gsuffix_routine_thread(param);
}

int gsuffix_me(FILE* stream_out, int n, char** mem_files, size_t* size_mem_files, int n_threads)
{
	return gsuffix_me_asm_impl(stream_out, n, mem_files, size_mem_files, n_threads);
}

int gsuffix_fget_size(FILE* mem)
{
	return gsuffix_fget_size_asm_impl(mem);
}

int gsuffix_fi(char* mem, size_t size_mem, FILE* stream_out)
{
	return gsuffix_fi_asm_impl(mem, size_mem, stream_out);
}

int gsuffix_co(char* mem, size_t size_mem, FILE* out)
{
	return gsuffix_co_asm_impl(mem, size_mem, out);
}

int gsuffix_so(char* mem, size_t size_mem, FILE* stream_out)
{
	return gsuffix_so_asm_impl(mem, size_mem, stream_out);
}

int gsuffix_somatic_bp(char* ref_file, FILE* normal_mem, FILE* cancer_mem, FILE* out_file)
{
	return gsuffix_somatic_bp_asm_impl(ref_file, normal_mem, cancer_mem, out_file);
}


int gsuffix_lookup(struct gsuffix *gs, const char *p, int m, int (*hit_callback)(int index, int pos, void *userdata), void *userdata)
{
	struct gsuffix_internal *gsi = (struct gsuffix_internal*)gs;

	switch (gsi->type)
	{
		case	GSUFFIX_SUFFIXARRAY:
			return genArray_lookup_from_char(gsi->ds.genarray,p,m,hit_callback,userdata);

		case	GSUFFIX_HASH:
			return gen_hash_lookup(gsi->ds.genhash, p, hit_callback, userdata);

		case 	GSUFFIX_LOTRUNCSUFFIXTREE:
			return lotruncatedSuffixTree_lookup(gsi->ds.lotruncsuffixtree, p, m, hit_callback, userdata);

/*		case	GSUFFIX_SHTRUNCSUFFIXTREE:
			return shtruncatedSuffixTree_lookup(kmers->ds.shtruncsuffixtree, p, m, hit_callback, userdata);
*/
		case	GSUFFIX_LOSUFFIXTREE:
			return lotruncatedSuffixTree_lookup(gsi->ds.losuffixtree, p, m, hit_callback, userdata);
/*
		case	GSUFFIX_SHSUFFIXTREE:
			return shtruncatedSuffixTree_lookup(kmers->ds.shsuffixtree, p, m, hit_callback, userdata);
*/
		case	GSUFFIX_DIRECT:
			return gen_direct_lookup(gsi->ds.gendirect, p, hit_callback, userdata);

		default:
			break;
	}
	return 0;
}

int gsuffix_lookup_exists(struct gsuffix *gs, const char *p, int m)
{
	struct gsuffix_internal *gsi = (struct gsuffix_internal*)gs;
	switch (gsi->type)
	{
		case	GSUFFIX_LOTRUNCSUFFIXTREE:
				return lotruncatedSuffixTree_lookup_exists(gsi->ds.lotruncsuffixtree, p, m);

		default:
				break;
	}
	
	return 0;
}

int gsuffix_enum_k_mers(struct gsuffix *gs, int k,
				       int (*callback)(char *kmer, int id, int pos, void *userdata),
				       void *userdata)
{
	struct gsuffix_internal *gsi = (struct gsuffix_internal*)gs;

	switch (gsi->type)
	{
	case GSUFFIX_LOTRUNCSUFFIXTREE:
	case	GSUFFIX_LOSUFFIXTREE:
			return lostree_enum_k_mers(gsi->ds.lotruncsuffixtree, k, callback, userdata);
			break;
/*
	case GSUFFIX_SHTRUNCSUFFIXTREE:
			return shstree_enum_k_mers(kmers->ds.shtruncsuffixtree, k, callback, userdata);
			break;
*/
		default:
			return 0;
	}
}

