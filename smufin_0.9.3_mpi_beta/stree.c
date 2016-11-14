#include <assert.h>
#include <string.h>  /* memset */
#include <stdlib.h>
#include <stdio.h>

#include "alphabet.h"
#include "gmem.h"
#include "stree.h"


#define ASSERT_IS_LEAF(node)			(assert((node)->isaleaf))
#define ASSERT_IS_INTERNAL_NODE(node) 	(assert(!((node)->isaleaf)))

#define USE_LINKEDLIST


#define MIN(a,b) ((a)<(b)?(a):(b))

#ifdef SHORTINTLEAF
#define FUNCTIONNAME(x) sh ## x
typedef unsigned short int_type;
#else
#define FUNCTIONNAME(x) lo ## x
typedef unsigned int int_type;
#endif


#ifndef USE_LINKEDLIST
#define DNA_SZ 4
#endif

struct node
{
	int isaleaf;
	char *edgestr;
	int edgelen;
	struct node *parent;
#ifdef USE_LINKEDLIST
	struct node *next; /* Pointer to next sibling */
#endif
};

struct int_leaf
{
  int_type strid, pos;
  struct int_leaf *next;
};

struct internal_node
{
	struct node n;
	struct node *suffix_link;

	struct int_leaf *leaves;

#ifdef USE_LINKEDLIST
	struct node *first; /* Pointer to first child */
#else
	struct node *children[DNA_SZ]; /* Array containing pointers to children */
#endif
};

#define INTERNAL_NODE_SET_SUFFIXLINK(node,sl) ((((struct internal_node*)node)->suffix_link) = (sl))
#define INTERNAL_NODE_GET_SUFFIXLINK(node) (((struct internal_node*)node)->suffix_link)

struct leaf_node
{
	struct node n;
	int strid, pos;
};

struct stree
{
  struct item_pool *node_pool;
  struct freeable_item_pool *leaf_node_pool;
  struct item_pool *int_leaf_pool;

  struct node *root;

  int num_nodes;

  char **strings;
  int *lengths, *ids;
  int alpha_size;
  int nextslot, strsize;

  int tree_size;
  int num_compares, edges_traversed, links_traversed;
  int child_cost, nodes_created, creation_cost;

  int (*hitcallback) (int, int, void *);    /* temporary storage for the hitcallback */
  int (*hitcallback3) (int,int,int,void *); /* Alternatve hit callback with 3 int args */
  void *userdata;			    /* temporary storage for the userdata */
  int hits;				    /* temporary storage for the hit counter */
  int k;				    /* the maximum depth of the tree */
};


#define intleaf_get_strid(ileaf) (ileaf->strid)

#define node_isaleaf(n) ((n)->isaleaf)
#define node_get_edgestr(node)  ((node)->edgestr)
#define node_get_edgelen(node)  ((node)->edgelen)
#define node_get_char(node)  (*((node)->edgestr))
#define node_get_parent(node) ((node)->parent)
#define leaf_get_strid(leaf) ((leaf)->strid)

#define stree_get_string(tree,id)  ((tree)->strings[(id)])
#define stree_get_length(tree,id)  ((tree)->lengths[(id)])
#define kTstree_get_max_depth(tree)  ((tree)->k)

static struct node *stree_edge_split(struct stree *tree,
				     struct node *node,
				     int len);
static struct stree *kTstree_new_tree(int number_of_strings,
				    int alphabet_size,int k);
static void node_reconnect(struct stree *stree, struct internal_node *parent,
                           struct node *oldchild, struct node *newchild);

static struct internal_node *new_internal_node(struct stree *stree,
											   char *edgestr,
											   int_type edgelen)
{
	struct internal_node *node;

	node = item_pool_alloc(stree->node_pool);

	memset(node, 0, sizeof(*node));
	node->n.edgestr = edgestr;
	node->n.edgelen = edgelen;
	return node;
}


static struct leaf_node *new_leaf(struct stree *st, int_type strid, char *edgestr,  int_type edgelen, int_type leafpos)
{
 	struct leaf_node *leaf;

	leaf = freeable_item_pool_alloc(st->leaf_node_pool);

	leaf->n.isaleaf = 1;
	leaf->n.edgestr = edgestr;
 	leaf->n.edgelen = edgelen;
 	leaf->n.parent = NULL;
#ifdef USE_LINKEDLIST
 	leaf->n.next = NULL;
#endif
	leaf->strid = strid;
	leaf->pos = leafpos;

	return leaf;
}

static void free_leaf(struct stree *st, struct leaf_node *leaf)
{
	freeable_item_pool_free(st->leaf_node_pool,leaf);
}

static int node_add_intleaf(struct stree *stree, struct node *n, int_type strid, int_type pos)
{
	struct int_leaf *intleaf;
	struct internal_node *node;

	ASSERT_IS_INTERNAL_NODE(n);

	node = (struct internal_node*)n;
	intleaf = item_pool_alloc(stree->int_leaf_pool);
	intleaf->strid = strid;
	intleaf->pos = pos;

	intleaf->next = node->leaves;
	node->leaves = intleaf;

	return 1;
}

struct node *node_find_child(struct node *node, char ch)
{
	if (node_isaleaf(node))
		return NULL;
	else
	{
#ifdef USE_LINKEDLIST
		struct node *child;
		for (child = ((struct internal_node*)node)->first;child != NULL;child = child->next)
		{
			if (node_get_char(child) == ch)
				return child;
		}
		return NULL;
#else
		struct internal_node *inode = (struct internal_node*)node;

		return inode->children[(unsigned char)ch];
#endif
	}
}

static struct node *convert_leaf_into_internal_node(struct stree *stree, struct node *node)
{
	struct internal_node *newnode;
	struct leaf_node *leaf;
	struct int_leaf *ileaf;

	ASSERT_IS_LEAF(node);
	leaf = (struct leaf_node*)node;
	newnode = new_internal_node(stree,leaf->n.edgestr,leaf->n.edgelen);

	ileaf = item_pool_alloc(stree->int_leaf_pool);
	ileaf->next = NULL;
	ileaf->strid = leaf->strid;
	ileaf->pos = leaf->pos;
	newnode->leaves = ileaf;

	node_reconnect(stree,(struct internal_node*)node->parent, node, &newnode->n);
	free_leaf(stree,leaf);

	return &newnode->n;
}

static struct node *node_connect(struct stree *stree, struct node *p, struct node *child)
{
	struct internal_node *parent;

	if (node_isaleaf(p))
	{
		p = convert_leaf_into_internal_node(stree,p);
	}

	child->parent = p;
	parent = (struct internal_node*)p;

#ifdef USE_LINKEDLIST
	if (parent->first)
	{
		child->next = parent->first;
	}
	parent->first = child;
#else
 	{
		unsigned char ch;
		ch = node_get_char(child);
		parent->children[ch] = child;
	}
#endif

  return p;
}

#ifdef USE_LINKEDLIST
static void node_disconnect(struct internal_node *parent, struct node *to_remove)
{
	struct node *child;
	struct node *prev_child = NULL;

	for (child = parent->first;child != NULL;child = child->next)
	{
		if (child == to_remove)
		{
			if (!prev_child)
				parent->first = to_remove->next;
			else
				prev_child->next = to_remove->next;

			to_remove->next = NULL;
			break;
		}
		prev_child = child;
	}
}
#endif

static void node_reconnect(struct stree *stree, struct internal_node *parent,
                           struct node *oldchild, struct node *newchild)
{
#ifdef USE_LINKEDLIST
	node_disconnect(parent,oldchild);
	node_connect(stree,&parent->n,newchild);
#else
	parent->children[(unsigned char) node_get_char(newchild)] = newchild;
#endif

	newchild->parent = &parent->n;
	oldchild->parent = NULL;
}

static struct node *stree_get_root(struct stree *tree)
{
  return tree->root;
}

static int stree_insert_string(struct stree *tree, char *S,
                            int M, int strid)
{
  int slot;

  slot = tree->nextslot;
  tree->strings[slot] = S;

  tree->lengths[slot] = M;
  tree->ids[slot] = strid;

  tree->nextslot++;

  return slot;
}

static int kTstree_ukkonen_add_string(struct stree *tree,
				      char *S, int M, int strid)
{
	int i, j, g, h, gprime, edgelen, id;
	int k, depth;
	char *edgestr;

	struct node *node, *lastnode, *root, *child, *parent;
	struct leaf_node *leaf;

	if ((id = stree_insert_string(tree, S, M, strid)) == -1)
		return 0;

	k = kTstree_get_max_depth(tree); 
	root = stree_get_root(tree);
	node = lastnode = root;
	g = 0;
	edgelen = 0;
	edgestr = NULL;
	depth = 0;

	for (i=0,j=0; i <= M; i++)
	{
		for ( ; j <= i && j < M; j++)
		{
			if (node == root)
			{
				assert(depth==0);
			}

			if (g == 0 || g == edgelen)
			{
				if (i < M)
				{
					if ((child = node_find_child(node, S[i])) != NULL)
					{
						node = child;
						g = 1;
						edgestr = node_get_edgestr(node);
						edgelen = node_get_edgelen(node);
						depth++;

						if (depth == k && M-i > 1)
						{
							if (node_isaleaf (node)) node = convert_leaf_into_internal_node(tree,node);
							if (!node_add_intleaf(tree, node, id, j)) return 0;

							i++;
						} else
						{
							break;
						}
					} else
					{
						leaf = new_leaf(tree, id, stree_get_string(tree,id)+i, MIN(k-depth, M-i), j);

						if  (node_isaleaf (node))
							node = convert_leaf_into_internal_node(tree,node);
						node = node_connect(tree, node, &leaf->n);
						tree->num_nodes++;
					}
				} else
				{
					if (node_isaleaf(node)) node = convert_leaf_into_internal_node(tree, node);
					if (!node_add_intleaf(tree, node, id, j)) return 0;
				}

				ASSERT_IS_INTERNAL_NODE(lastnode);
				if (lastnode != root && !INTERNAL_NODE_GET_SUFFIXLINK(lastnode))
					INTERNAL_NODE_SET_SUFFIXLINK(lastnode,node);

				lastnode = node;
			} else
			{
				if (i < M && S[i] == edgestr[g])
				{
					g++;
					depth++;

					if (depth == k && M-i > 1)
					{
						if (node_isaleaf(node)) node = convert_leaf_into_internal_node(tree,node);
						if (!node_add_intleaf(tree, node, id, j)) return 0;
						i++;
					} else
					{
						break;
					}
				} else 
				{
					if ((node = stree_edge_split (tree, node, g)) == NULL)
						return 0;

					tree->num_nodes++;
					edgestr = node_get_edgestr(node);
					edgelen = node_get_edgelen(node);

					if (i < M)
					{
						leaf = new_leaf(tree, id, stree_get_string(tree,id)+i, MIN(k-depth, M-i), j);

						if  (node_isaleaf (node))
							node = convert_leaf_into_internal_node(tree,node);

						node = node_connect(tree, node, &leaf->n);
						tree->num_nodes++;
					}   else
					{
						if  (node_isaleaf (node))
							node = convert_leaf_into_internal_node(tree,node);

						if (!node_add_intleaf(tree, node, id, j))
							return 0;
					}
				}

				ASSERT_IS_INTERNAL_NODE(lastnode);
				if (lastnode != root && !INTERNAL_NODE_GET_SUFFIXLINK(lastnode))
					INTERNAL_NODE_SET_SUFFIXLINK(lastnode,node);
				lastnode = node;
			}

			if (node == root);
			else if (g == edgelen && INTERNAL_NODE_GET_SUFFIXLINK(node) && depth < k)
			{
				node = INTERNAL_NODE_GET_SUFFIXLINK(node);
				edgestr = node_get_edgestr(node);
				edgelen = node_get_edgelen(node);
				g = edgelen;
				depth--;
				continue; /* Next extension! */
			} else
			{
				depth--;
				parent = node_get_parent(node);
				if (parent != root)
					node = INTERNAL_NODE_GET_SUFFIXLINK(parent);
				else
				{
					node = root;
					g--;
				}
				edgelen = node_get_edgelen(node);

				h = i - g;

				while (g > 0)
				{
					node = node_find_child(node, S[h]);
					gprime = node_get_edgelen(node);
					if (gprime > g)
						break;

					g -= gprime;
					h += gprime;
				}

				edgestr = node_get_edgestr(node);
				edgelen = node_get_edgelen(node);

				if (g == 0)
				{
					if (lastnode != root && !node_isaleaf(node) && !INTERNAL_NODE_GET_SUFFIXLINK(lastnode))
					{
						INTERNAL_NODE_SET_SUFFIXLINK(lastnode,node);
						lastnode = node;
					}

					if (node != root)
						g = edgelen;
				}
			}
		}
	}
	return 1;
}

static int stree_ukkonen_add_string(struct stree *tree, char *S, int M, int strid)
{
	int i, j, g, h, gprime, edgelen, id;
	char *edgestr;
	struct node *node, *lastnode, *root, *child, *parent;
	struct leaf_node *leaf;

	if ((id = stree_insert_string(tree, S, M, strid)) == -1)
		return 0;

	root = stree_get_root(tree);
	node = lastnode = root;

	g = 0;

	edgelen = 0;
	edgestr = NULL;

	for (i=0,j=0; i <= M; i++)
	{
		for ( ; j <= i && j < M; j++)
		{
			if (g == 0 || g == edgelen)
			{
				if (i < M)
				{
					if ((child = node_find_child(node, S[i])) != NULL)
					{
						node = child;
						g = 1;
						edgestr = node_get_edgestr(node);
						edgelen = node_get_edgelen(node);
						break;
					}

					leaf = new_leaf(tree, id, stree_get_string(tree,id)+i, stree_get_length(tree,id)-i, j);
					node = node_connect(tree, node, (struct node*) leaf);
					tree->num_nodes++;
				} else
				{
           			if (node_isaleaf (node)) node = convert_leaf_into_internal_node(tree, node);
					if (!node_add_intleaf(tree, node, id, j)) return 0;
				}

				ASSERT_IS_INTERNAL_NODE(lastnode);
				if (lastnode != root && !INTERNAL_NODE_GET_SUFFIXLINK(lastnode))
					INTERNAL_NODE_SET_SUFFIXLINK(lastnode,node);
				lastnode = node;
			} else 
			{
				if (i < M && S[i] == edgestr[g])
				{
					g++;
					break;
				}

				if ((node = stree_edge_split (tree, node, g)) == NULL)
					return 0;

				edgestr = node_get_edgestr(node);
				edgelen = node_get_edgelen(node);

				if (i < M)
				{
					leaf = new_leaf(tree, id, stree_get_string(tree,id)+i, stree_get_length(tree,id)-i, j);
					node = node_connect(tree, node, (struct node*) leaf);
					tree->num_nodes++;
				} else
				{
					if (node_isaleaf (node)) node = convert_leaf_into_internal_node(tree, node);
					if (!node_add_intleaf(tree, node, id, j)) return 0;
				}
				ASSERT_IS_INTERNAL_NODE(lastnode);
				if (lastnode != root && !INTERNAL_NODE_GET_SUFFIXLINK(lastnode))
					INTERNAL_NODE_SET_SUFFIXLINK(lastnode,node);
				lastnode = node;
			}

			if (node == root);
			else if (g == edgelen && INTERNAL_NODE_GET_SUFFIXLINK(node) != NULL)
			{
				node = INTERNAL_NODE_GET_SUFFIXLINK(node);
				edgestr = node_get_edgestr(node);
				edgelen = node_get_edgelen(node);
				g = edgelen;
			} else
			{
        			parent = node_get_parent(node);
				if (parent != root)
					node = INTERNAL_NODE_GET_SUFFIXLINK(parent);
				else
				{
					node = root;
					g--;
				}
				edgelen = node_get_edgelen(node);

				h = i - g;

				while (g > 0)
				{
					node = node_find_child(node, S[h]);
					gprime = node_get_edgelen(node);
					if (gprime > g)
						break;

					g -= gprime;
					h += gprime;
				}

				edgestr = node_get_edgestr(node);
				edgelen = node_get_edgelen(node);

				if (g == 0)
				{
					if (lastnode != root &&
					    !node_isaleaf(node) &&
					    INTERNAL_NODE_GET_SUFFIXLINK(lastnode) == NULL)
					{
						INTERNAL_NODE_SET_SUFFIXLINK(lastnode,node);
						lastnode = node;
					}

					if (node != root)
						g = edgelen;
				}
			}
		}
	}
	return 1;
}

static void kTstree_traverse_subtree(struct stree *tree, struct node *node,
                            int (*preorder_fn)(), int (*postorder_fn)());


static struct stree *kTstree_new_tree(int number_of_strings,
									  int alphabet_size, int k)
{
	struct stree *tree;
	int i;

	tree = gsuffix_malloc(sizeof(*tree));
	memset(tree, 0, sizeof(*tree));

	tree->node_pool = item_pool_create(sizeof(struct internal_node));
	tree->int_leaf_pool = item_pool_create(sizeof(struct int_leaf));
	tree->leaf_node_pool = freeable_item_pool_create(sizeof(struct leaf_node));

	tree->k = k;
	tree->alpha_size = alphabet_size;

	tree->strsize = number_of_strings;
	tree->strings = gsuffix_malloc(tree->strsize * sizeof(char *));
	tree->lengths = gsuffix_malloc(tree->strsize * sizeof(int));
	tree->ids = gsuffix_malloc(tree->strsize * sizeof(int));

	for (i = 0; i < tree->strsize; i++)
	{
		tree->strings[i] = NULL;
		tree->lengths[i] = tree->ids[i] = 0;
	}

	tree->nextslot = 0;			
	tree->root = &new_internal_node(tree, NULL, 0)->n;
	tree->num_nodes = 1;		
	return tree;
}


static struct stree *stree_new_tree(int number_of_strings,int alphabet_size)
{
  return kTstree_new_tree(number_of_strings,alphabet_size,-1);
}

static struct node *stree_edge_split(struct stree *tree,
				       struct node *node, int len)
{
  struct internal_node *newnode, *parent;
  if (node == stree_get_root(tree) ||
      len == 0 || node_get_edgelen(node) <= len)
    return NULL;

  newnode = new_internal_node(tree, node->edgestr, len);
  if (newnode == NULL)
    return NULL;

  parent = (struct internal_node*)node_get_parent(node);
  node_reconnect(tree, parent, node, &newnode->n);

  node->edgestr += len;
  node->edgelen -= len;

  if (node_connect(tree, &newnode->n, node) == NULL)
  {
    fprintf(stderr,"Error in reconnecting node at %s (%d), terminating program.\n",
	    __FILE__,__LINE__);
      
  }

  return &newnode->n;
}








#if 0
static void debugPrintTruncatedTree(TRUNC_SUFFIX_TREE tree)
{
	printf("%s l. %d: debugPrintTree\n",__FILE__,__LINE__);
	printf("\tnum_nodes:\t\t%d\n",tree->num_nodes);
	printf("\troot node:\n*****\n");
	debugPrintNode(tree->root,0);
	printf("*****\n");
	printf("\ttree_size:%d\n",tree->tree_size);
	printf("\tstrsize:%d\n",tree->strsize);
	printf("\tstrings entered in tree structure:\n");
	/*for (i=0;i<tree->nextslot;++i) {
	*	printf("\t\tstring %d: %s\n", i ,tree->strings[i]);
	}*/
	printf("\tnextslot: %d\n",tree->nextslot);
}

#endif

/*
 * debugPrintNode
 *
 * print out node values
 *************************/
#if 0
static void debugPrintNode(NODE node, int depth)
{
	INTLEAF ileaf;
	LEAF leaf;

	if (!node) return;

	if (depth  == 0)
	{
		printf("root\n");
	}
/*	printf("%s l. %d: debugPrintNode\n",__FILE__,__LINE__);
	printf("depth: %d,\t isaleaf:\t%d;\n",depth,node_isaleaf(node));
	printf("\tedgestr: %s\t\n",node->edgestr);
	printf("\tedgelen:%d\n",node->edgelen);*/
	if (!node_isaleaf(node)){
		if (node->leaves != NULL){
			printf("depth %d ,nodeedgelen: %d, char %c (ILEAF)",depth,node->edgelen,*(node->edgestr)+48);
		ileaf = node->leaves;
		while(ileaf!=NULL){
			printf(" id: %d, pos: %d\t", ileaf->strid ,ileaf->pos);
			ileaf = ileaf->next;
		}
		printf("\n");
		}
		   debugPrintNode(node_find_child(node, 0),depth+1);

		   debugPrintNode(node_find_child(node, 1),depth+1);

		   debugPrintNode(node_find_child(node, 2),depth+1);

		   debugPrintNode(node_find_child(node, 3),depth+1);


	}
	else{
		leaf = (LEAF) node;
		printf("depth: %d, strid: %d, pos: %d, char %c ,edgelen: %d (LEAF)\n",depth,leaf->strid,leaf->pos,*(leaf->edgestr)+48,leaf->edgelen);
	}
}

void FUNCTIONNAME(debugPrintTruncatedTreeRecursive)(TRUNC_SUFFIX_TREE tree)
{
	int depth = 0;
	printf("%s l. %d: debugPrintTreeRecursive\n",__FILE__,__LINE__);
	printf("\tnum_nodes:\t\t%d\n",tree->num_nodes);
	printf("number of seqs in the tree: %d\n",tree->nextslot);
	printf("Print out the whole tree:\n");
	debugPrintNode(tree->root,depth);
}

#endif


static int myDebugPrintNode(struct node *n)
{
	int i;
	char *p;
	printf("Node at %p\n",n);
	if (n->isaleaf)
		printf("\tNode is a leaf\n");
	else
		printf("\tNode is a nonleaf\n");
	printf("\tEdgelen=%d\tedgestr=",n->edgelen);
	p=n->edgestr;
	for (i=0;i<n->edgelen;++i) { putchar(*p); p++; }
	printf("\n");
	
	
}

static int kTstree_walk_to_leaf(struct stree *tree, struct node *node, int pos,
                           const char *T, int N, struct node **node_out, int *pos_out)
{
	int len, edgelen;
	char *edgestr;
	struct node *child;

	if (node_isaleaf(node))
	{
		*node_out = node;
		*pos_out = pos;
		return 0;
	}

	edgestr = node_get_edgestr(node);
	edgelen = node_get_edgelen(node);

	len = 0;
	while (1)
	{
		while (len < N && pos < edgelen && T[len] == edgestr[pos])
		{
			pos++;
			len++;
		}

		if (len == N || pos < edgelen || (child = node_find_child(node, T[len])) == NULL) {
		  break;
		}

		if (node_isaleaf(child))
		{
			*node_out = child;
			*pos_out = 0; 
			return len;
		}

		node = child;
		edgestr = node_get_edgestr(node);
		edgelen = node_get_edgelen(node);
		pos = 1;
		len++;
	}

	*node_out = node;
	*pos_out = pos;
	return len;
}

static int kTstree_walk(struct stree *tree, struct node *node,
			int pos, const char *T, int N,
               struct node **node_out, int *pos_out)
{
  int len, endpos, edgelen;
  char *edgestr;
  struct node *endnode;

  len = kTstree_walk_to_leaf(tree, node, pos, T, N, &endnode, &endpos);

  if (!node_isaleaf(endnode) || len == N)
  {
    *node_out = endnode;
    *pos_out = endpos;
    return len;
  }

  edgestr = node_get_edgestr(endnode);
  edgelen = node_get_edgelen(endnode);

  while (len < N && endpos < edgelen && T[len] == edgestr[endpos])
  {
    len++;
    endpos++;
  }

  *node_out = endnode;
  *pos_out = endpos;
  return len;
}

static int kTstree_match(struct stree *tree, const char *T, int N,
                struct node **node_out, int *pos_out)
{
  return kTstree_walk(tree, stree_get_root(tree), 0, T, N, node_out, pos_out);
}

static int kTstree_get_num_children(struct stree *tree, struct node *node)
{
#ifdef USE_LINKEDLIST
  int count;
  struct node *child;

  if (node_isaleaf(node)) return 0;

  for (child = ((struct internal_node*)node)->first,count=0;child != NULL;child = child->next)
	count++;

  return count;
#else
	int i, count;
	struct node **children;

	if (node_isaleaf(node))
		return 0;

	count = 0;
	children = ((struct internal_node*)node)->children;

    for (i=0; i < tree->alpha_size; i++)
    	if (children[i] != NULL)
			count++;

	return count;
#endif
}

#if 0
void kTstree_traverse(TRUNC_SUFFIX_TREE tree, int (*preorder_fn)(),
                    int (*postorder_fn)())
{
  kTstree_traverse_subtree(tree, kTstree_get_root(tree), preorder_fn,
                         postorder_fn);
}
#endif


static void kTstree_traverse_subtree(struct stree *tree, struct node *node,
                            int (*preorder_fn)(), int (*postorder_fn)())
{
  enum { START, FIRST, MIDDLE, DONE, DONELEAF } state;
  int i, num, childnum;
  struct node *root, *child;

  root = node;
  state = START;
  while (1) {
    if (state == START) {
      if (preorder_fn != NULL)
        if(!(*preorder_fn)(tree, node))
           return;

      num = kTstree_get_num_children(tree, node);
      if (num > 0)
        state = FIRST;
      else
        state = DONELEAF;
    }

    if (state == FIRST || state == MIDDLE) {
      if (state == FIRST)
        childnum = 0;
      else
        childnum = node->isaleaf;

#ifdef USE_LINKEDLIST
      for (child = ((struct internal_node*)node)->first,i=0;child != NULL && i<childnum;child = child->next,i++);
#else
	{
		struct node **children;

        children = ((struct internal_node*)node)->children;
        for (i=childnum; i < tree->alpha_size; i++) {
          if (children[i] != NULL)
            break;

        }
        child = (i < tree->alpha_size ? children[i] : NULL);
	}
#endif

      if (child == NULL)
        state = DONE;
      else {
        node->isaleaf = i + 1;
        node = child;

        state = START;
      }
    }

    if (state == DONE  || state == DONELEAF) {
      if (state == DONE)
        node->isaleaf = 0;

      if (postorder_fn != NULL)
        (*postorder_fn)(tree, node);

      if (node == root)
        break;

      node = node_get_parent(node);
      state = MIDDLE;
    }
  }
}

static int add_match(struct stree *tree, struct node *node)
{
  struct leaf_node *leaf;
  struct int_leaf *intleaf;
  int pos, id;

  if (node_isaleaf(node))
    {
      leaf = (struct leaf_node*)node;
      pos = leaf->pos;

      id = leaf_get_strid(leaf);
      tree->hits++;

      return tree->hitcallback(id, pos, tree->userdata);
    }  else
    {
      intleaf = ((struct internal_node*)node)->leaves;

      for (;intleaf != NULL;intleaf=intleaf->next)
	{
	  pos = intleaf->pos;
	  id = leaf_get_strid(intleaf);

	  tree->hits++;
	  if (!tree->hitcallback(id, pos, tree->userdata))
	    return 0;
	}
    }
  return 1;
}

int FUNCTIONNAME(truncatedSuffixTree_lookup)(struct stree *tree, const char *pattern, int pattern_length, int (*hitcallback)(int index, int pos, void *userdata), void *userdata)
{
	int pos,matchlen;
	struct node *node;
	char *mapped_pattern;

#ifndef USE_LINKEDLIST
	char buf[16];
#endif

	tree->hitcallback = hitcallback;
	tree->userdata = userdata;
	tree->hits = 0;

#ifndef USE_LINKEDLIST
	if (pattern_length < sizeof(buf)) mapped_pattern = buf;
	else mapped_pattern = gsuffix_malloc(pattern_length*sizeof(char));
	memcpy(mapped_pattern,pattern,pattern_length);
	map_string(mapped_pattern,pattern_length);
#else
	mapped_pattern = (char*)pattern;
#endif

	node = stree_get_root(tree);

	matchlen = kTstree_match(tree, mapped_pattern, pattern_length, &node, &pos);

	if (matchlen == pattern_length)
		kTstree_traverse_subtree(tree, node, add_match, (int (*)()) NULL);

#ifndef USE_LINKEDLIST
	if (pattern_length >= sizeof(buf))
		gsuffix_free(mapped_pattern);
#endif

	return tree->hits;
}



int FUNCTIONNAME(truncatedSuffixTree_lookup_exists)(struct stree *tree, const char *pattern, int pattern_length)
{
	int pos,matchlen;
	struct node *node;
	char *mapped_pattern;

#ifndef USE_LINKEDLIST
	char buf[16];
#endif

#ifndef USE_LINKEDLIST
	if (pattern_length < sizeof(buf)) mapped_pattern = buf;
	else mapped_pattern = gsuffix_malloc(pattern_length*sizeof(char));
	memcpy(mapped_pattern,pattern,pattern_length);
	map_string(mapped_pattern,pattern_length);
#else
	mapped_pattern = (char*)pattern;
#endif

	node = stree_get_root(tree);

	matchlen = kTstree_match(tree, mapped_pattern, pattern_length, &node, &pos);

   #ifndef USE_LINKEDLIST
	if (pattern_length >= sizeof(buf))
		gsuffix_free(mapped_pattern);
#endif

	if (matchlen == pattern_length)
		return 1;
	else
		return 0;

}

struct stree *FUNCTIONNAME(kTstree_create_from_char)
	(char **strings, int nstrings, int k)
{
	int i, len;
	char *str;
	struct stree *tree;

	init_alphabet(DNA);

	if (!(tree = kTstree_new_tree(nstrings, DNA_SZ, k)))
		return NULL;

	for (i = 0; i < nstrings; ++i)
	{
		len = strlen(strings[i]);

#ifndef USE_LINKEDLIST
		str = create_mapped_string(strings[i], len);
#else
		str = strings[i];
#endif

		if (kTstree_ukkonen_add_string(tree, str, len, i + 1) < 1)
		{
			FUNCTIONNAME(kTstree_delete_tree)(tree);
			return NULL;
		}
	}
	return tree;
}

struct stree *FUNCTIONNAME(stree_create_from_char) (char **strings, int nstrings)
{
	int i,len;
	char *str;
	struct stree *tree;

	init_alphabet(DNA);

	if (!(tree = stree_new_tree(nstrings,DNA_SZ)))
		return NULL;

	for (i = 0; i < nstrings; ++i)
	{
		len = strlen(strings[i]);

#ifndef USE_LINKEDLIST
		str = create_mapped_string(strings[i], len);
#else
		str = strings[i];
#endif
		if (stree_ukkonen_add_string(tree, str, len, i+1) < 1)
		{
			FUNCTIONNAME(stree_delete_tree)(tree);
			return NULL;
		}
	}
	return tree;
}

void FUNCTIONNAME(stree_delete_tree)(struct stree *tree)
{
	item_pool_delete(tree->node_pool);
	item_pool_delete(tree->int_leaf_pool);
	freeable_item_pool_delete(tree->leaf_node_pool);
#ifndef USE_LINKEDLIST
	{
		int i;
		for (i = 0; i < tree->strsize; i++)
			gsuffix_free(tree->strings[i]);
	}
#endif
	gsuffix_free(tree->strings);
	gsuffix_free(tree->lengths);
	gsuffix_free(tree->ids);
	gsuffix_free(tree);
}

void FUNCTIONNAME(kTstree_delete_tree)(struct stree *tree)
{
	FUNCTIONNAME(stree_delete_tree)(tree);
}

static int kTstree_enum_k_mers_recursive(struct stree *tree,
		struct node *node, int depth,
		int (*callback)(char *kmer, int id, int pos, void *userdata), void *userdata)
{
	struct internal_node *inode;
	struct node *n;

	depth += node->edgelen;

	if (node_isaleaf(node))
	{
		if (depth == tree->k)
		{
			struct leaf_node *leaf;
			leaf = (struct leaf_node*)node;
			callback(tree->strings[leaf->strid]+leaf->pos,leaf->strid,leaf->pos,userdata);
		}
		return 1;
	}

	inode = (struct internal_node*)node;
	if (!(n = inode->first))
	{
		if (depth == tree->k)
		{
			struct int_leaf *ileaf = inode->leaves;
			char *kmer = tree->strings[ileaf->strid]+ileaf->pos;

			while (ileaf)
			{
				callback(kmer,ileaf->strid,ileaf->pos,userdata);
				ileaf = ileaf->next;
			}
		}
		return 1;
	}

	while (n)
	{
		if (!kTstree_enum_k_mers_recursive(tree, n, depth, callback,userdata))
			return 0;
		n = n->next;
	}
	return 1;
}




int FUNCTIONNAME(stree_enum_k_mers)(struct stree *tree, int k,
				       int (*callback)(char *kmer, int id, int pos, void *userdata),
				       void *userdata)
{
#ifdef USE_LINKEDLIST
	if (k != tree->k)
		return 0;

	kTstree_enum_k_mers_recursive(tree, tree->root, 0, callback, userdata);

	callback(NULL,0,0,userdata);

	return 1;
#else
	fprintf(stderr,"enumeration search not implemented for array based suffix tree (%s,%d)!",
		__FILE__,__LINE__);
	exit(-1);
	return 0;
#endif
}

static int add_match3(struct stree *tree, struct node *node, int X)
{
  struct leaf_node *leaf;
  struct int_leaf *intleaf;
  int pos, id;

  if (node_isaleaf(node))
    {
      leaf = (struct leaf_node*)node;
      pos = leaf->pos;

      id = leaf_get_strid(leaf);
      tree->hits++;

      return tree->hitcallback3(id, pos, X,tree->userdata);
    }  else
    {
      intleaf = ((struct internal_node*)node)->leaves;

      for (;intleaf != NULL;intleaf=intleaf->next)
	{

	  pos = intleaf->pos;
	  id = leaf_get_strid(intleaf);

	  tree->hits++;
	  if (!tree->hitcallback3(id, pos, X,tree->userdata))
	    return 0;

	}
    }
  return 1;
}

static void kTstree_traverse_subtree3(struct stree *tree, struct node *node,
				      int (*preorder_fn)(), int X)
{
  enum { START, FIRST, MIDDLE, DONE, DONELEAF } state;
  int i, num, childnum;
  struct node *root, *child;

  root = node;
  state = START;
  while (1) {
    if (state == START) {
      if (preorder_fn != NULL)
        if(!(*preorder_fn)(tree, node,X))
           return;

      num = kTstree_get_num_children(tree, node);
      if (num > 0)
        state = FIRST;
      else
        state = DONELEAF;
    }

    if (state == FIRST || state == MIDDLE) {
      if (state == FIRST)
        childnum = 0;
      else
        childnum = node->isaleaf;

#ifdef USE_LINKEDLIST
      for (child = ((struct internal_node*)node)->first,i=0;
	   child != NULL && i<childnum;
	   child = child->next, i++)  
#else
	{
	  struct node **children;
	  children = ((struct internal_node*)node)->children;
	  for (i=childnum; i < tree->alpha_size; i++) {
	    if (children[i] != NULL)
	      break;
	  }
	  child = (i < tree->alpha_size ? children[i] : NULL);
	}
#endif

	if (child == NULL)
	  state = DONE;
      else {
        node->isaleaf = i + 1; 
        node = child;
        state = START;
      }
    }

    if (state == DONE  || state == DONELEAF) {
      if (state == DONE)
        node->isaleaf = 0;
      if (node == root)
        break;

      node = node_get_parent(node);
      state = MIDDLE;
    }
  }
}



int stree_get_length_of_longest_match(struct stree *tree,const char *T, int N)
{
	int len, endpos;
	struct node *endnode, *root;
	root=stree_get_root(tree);
	len = kTstree_walk_to_leaf(tree,root, 0, T,N,&endnode,&endpos);
	return len;
}

