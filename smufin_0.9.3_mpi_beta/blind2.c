#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "common.h"

extern UChar  *Text;                   // input string+ overshoot
extern Int32  Text_size;               // size of input string
extern UChar  *Upper_text_limit;       // Text+Text_size

typedef struct nodex {
  Int32 skip;
  UChar key;  
  struct nodex  *down;      // first child
  struct nodex *right;      // next brother
} node;


#define BUFSIZE 1000
#define FREESIZE 5000
void *freearr[FREESIZE];
static node *bufn;
static int bufn_num=0, free_num=0;
static Int32 *Aux, Aux_written;
node **Stack;
int Stack_size;

static int neg_integer_cmp(const void *, const void *);
static node *find_companion(node *head, UChar *s);
static void insert_suffix(node *h, Int32 suf, int n, UChar mmchar);
static void traverse_trie(node *h);
static Int32 compare_suffixes(Int32 suf1, Int32 suf2, Int32 depth);
static void free_node_mem();
static node *get_leaf(node *head);

void blind_ssort(Int32 *a, Int32 n, Int32 depth)
{
  Int32 i,j,aj,lcp;
  node nh, *root, *h;

  qsort(a,n, sizeof(Int32), neg_integer_cmp);
  
for(j=0;j<n;j++)
    if(a[j]+depth < Text_size)
      break;
  if(j>=n-1) return;  // everything is already sorted!

  Stack = (node **) malloc(n*sizeof(node *));
  if(Stack==NULL) {
    fprintf(stderr,"Out of memory! (blind_ssort)\n");
    exit(1);
  }

  nh.skip = -1;   nh.right = NULL; nh.down = (void *) a[j]; 
  root = &nh;

  for(i=j+1;i<n;i++) {
    h=find_companion(root, Text+a[i]);
    assert(h->skip==-1);
    assert(Stack_size<=i-j);
    aj=(Int32) h->down;
    assert(aj>a[i]);
    lcp = compare_suffixes(aj,a[i],depth);
    insert_suffix(root, a[i], lcp, Text[aj+lcp]);
  }

  Aux=a;  Aux_written = j;
  traverse_trie(root);
  assert(Aux_written==n);
 
  free_node_mem();
  free(Stack);
}

static node *find_companion(node *head, UChar *s)
{
  UChar c;
  node *p;
  int t;

  Stack_size = 0;                // init stack
  while(head->skip >= 0) {
    Stack[Stack_size++] = head;
    t = head->skip;
    if(s+t>=Upper_text_limit)    // s[t] does not exist: mismatch 
      return get_leaf(head);
    c = s[t]; p = head->down;
  repeat:
    if(c==p->key) {              // found branch corresponding to c
      head = p;
      continue;
    }
    else if(c<p->key)            // no branch corresponding to c: mismatch
      return get_leaf(head);
    if((p=(p->right))==NULL)     // no other branches: mismatch
      return get_leaf(head);
    goto repeat;                 // look at next branch
  }
  Stack[Stack_size++] = head;
  return head;
}

static node *get_leaf(node *head)
{
  assert(head->skip>=0);

  do {
    head = head->down;
  } while(head->skip>=0);
  return head;
}

__inline__ node *new_node__blind_ssort(void)
{
  if(bufn_num-- == 0) {
    bufn = (node *) malloc(BUFSIZE * sizeof(node));
    if(bufn==NULL) {
      fprintf(stderr,"Out of mem (new_node1)\n"); exit(1);}
    freearr[free_num++] = (void *) bufn; 
    if(free_num>=FREESIZE) {
      fprintf(stderr,"Out of mem (new_node2)\n"); exit(1);}
   bufn_num = BUFSIZE-1;
  }
  return bufn++;
}


static void insert_suffix(node *h, Int32 suf, int n, UChar mmchar)
{
  Int32 t;
  UChar c, *s;
  node *p, **pp;

  s = Text + suf;

#if 0
  // ---------- find the insertion point
  while( (t=h->skip) < n) {
    if( t < 0) break;          // insert "suf" just above *h
    c=s[t];  p=h->down;  
  repeat:
    if(c==p->key) {              // found branch corresponding to c
      h = p;                     // go down
      continue;
    }
    if(c>p->key && p->right!=NULL) {
      p=p->right;
      goto repeat;
    }
    fprintf(stderr,"Error in blind_sort (insert_string)\n");
    exit(1);
  }
#else
  for(t=0;t<Stack_size;t++) {
    h=Stack[t];
    if(h->skip<0 || h->skip>=n) break;
  }  
#endif
  
  assert(s[n]!=mmchar || h->skip==-1 || h->skip==n);

  if(h->skip!=n) {
    p = new_node__blind_ssort();     // create and init new node
    p->key = mmchar;
    p->skip = h->skip;  // p inherits skip and children of *h
    p->down = h->down;   
    p->right = NULL;
    h->skip = n;
    h->down = p;        // now *h has p as the only child 
  }
  assert(h->skip==n);

  c=s[n]; pp = &(h->down);
  while((*pp)!=NULL) {
    if((*pp)->key>=c) 
      break;
    pp = &((*pp)->right);
  }
  p = new_node__blind_ssort();
  p->skip = -1;
  p->key = c; 
  p->right = *pp; *pp = p;
  p->down = (void *) suf;
  return;
}

static void traverse_trie(node *h)
{
  node *p, *nextp;

  if(h->skip<0)
    Aux[Aux_written++] = (Int32) h->down;
  else {
    p = h->down;
    assert(p!=NULL);
    do {
      nextp = p->right;
      if(nextp!=NULL) {
	assert(nextp->key>=p->key);
	if(nextp->key==p->key) {
	  traverse_trie(nextp);
	  traverse_trie(p);
	  p = nextp->right;
	  continue;
	}
      }
      traverse_trie(p);
      p=nextp;
    } while(p!=NULL);
  }
}

static __inline__
Int32 get_lcp_unrolled(UChar *b1, UChar *b2, Int32 cmp_limit)
{
  Int32 cmp2do; 
  UChar c1, c2;
  assert(b1 != b2);

  cmp2do = cmp_limit;
  do {
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      break;}
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  1; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  2; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  3; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  4; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  5; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  6; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  7; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  8; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  9; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 10; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 11; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 12; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 13; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 14; break; }
    b1++; b2++; 
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 15; break; }
    b1++; b2++; 

    cmp2do -= 16;
  } while(cmp2do>0);


  if(cmp_limit - cmp2do < cmp_limit)
    return cmp_limit-cmp2do;

  return cmp_limit-1;
} 



static Int32 compare_suffixes(Int32 suf1, Int32 suf2, Int32 depth)
{
  int limit;
  UChar *s1, *s2;

  assert(suf1>suf2);
  s1  = Text + depth +suf1;
  s2  = Text + depth +suf2;
  limit = Text_size - suf1 - depth;
  return depth + get_lcp_unrolled(s1 ,s2, limit);
}

static int neg_integer_cmp(const void *a, const void *b)
{
  return *((Int32 *) b) -  *((Int32 *) a); 
}

static void free_node_mem()
{
  int i;

  for(i=free_num-1;i>=0;i--) {
    assert(freearr[i]!=NULL);
    free(freearr[i]);
  }
  bufn_num=free_num=0;
}

