#ifndef KTSTREE_H_
#define KTSTREE_H_

struct stree;


struct stree *shkTstree_create_from_char(char **strings, int nstrings, int query_length);
int shtruncatedSuffixTree_lookup_exists(struct stree *tree, const char *pattern, int pattern_length);
int shstree_enum_k_mers(struct stree *tree, int k,
				       int (*callback)(char *kmer, int id, int pos, void *userdata),
				       void *userdata);
int shsuffixTree_lookup(struct stree *tree,
			const char *pattern,
			int pattern_length,
			int (*hitcallback)(int index,
					   int pos,
					   void *userdata),
			void *userdata);


struct stree *lokTstree_create_from_char(char **strings, int nstrings, int query_length);
void lokTstree_delete_tree(struct stree *tree);
int lotruncatedSuffixTree_lookup(struct stree *tree,
				 const char *pattern,
				 int pattern_length,
				 int (*hitcallback)(int index,
						    int pos,
						    void *userdata),
				 void *userdata);
int lotruncatedSuffixTree_lookup_exists(struct stree *tree, const char *pattern, int pattern_length);
int lostree_enum_k_mers(struct stree *tree, int k,
				       int (*callback)(char *kmer, int id, int pos, void *userdata),
				       void *userdata);

struct stree *lostree_create_from_char (char **strings, int nstrings);
void lostree_delete_tree(struct stree *tree);
int losuffixTree_lookup(struct stree *tree,
			const char *pattern,
			int pattern_length,
			int (*hitcallback)(int index,
					   int pos,
					   void *userdata),
			void *userdata);


#endif
