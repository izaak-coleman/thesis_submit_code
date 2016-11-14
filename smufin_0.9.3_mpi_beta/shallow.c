#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "common.h"

extern UChar *Text;                // start of input string
extern UChar *Upper_text_limit;    // Text+Text_size
extern Int32 _ds_Word_size;        // # of bytes in a word in mkq
extern Int32 Mk_qs_thresh;         // recursion limit for mk quicksort: 
                                   // groups smaller than this are sorted 
                                   // using insertion sort
Int32 Shallow_limit;               // Max depth for shallow sorting
UChar *Shallow_text_limit;         // Text+Shallow_limit

#define UNROLL 1                   // if !=0 partially unroll shallow_mkq

void helped_sort(Int32 *a, Int32 n, Int32 depth);
static void shallow_inssort_lcp(Int32 *a, Int32 n, UChar *text_depth);

static void shallow_mkq(Int32 *a, int n, UChar *text_depth);
static void shallow_mkq16(Int32 *a, int n, UChar *text_depth);
static void shallow_mkq32(Int32 *a, int n, UChar *text_depth);

void shallow_sort(Int32 *a, int n, int shallow_limit) 
{ 
  Shallow_limit = shallow_limit;        
  Shallow_text_limit = Text + shallow_limit;
  switch(_ds_Word_size) {
  case(1): shallow_mkq(a, n, Text+2); break;
  case(2): shallow_mkq16(a, n, Text+2); break;
  case(4): shallow_mkq32(a, n, Text+2); break;
  default:
    fprintf(stderr,
	    "Invalid word size for mkqs (%d) (shallow_sort)\n",_ds_Word_size);
    exit(1);
  }     
}


__inline__ void vecswap2(Int32 *a, Int32 *b, int n)
{   while (n-- > 0) {
        Int32 t = *a;
        *a++ = *b;
        *b++ = t;
    }
}

#define swap2(a, b) { t = *(a); *(a) = *(b); *(b) = t; }
#define ptr2char(i) (*(*(i) + text_depth))

__inline__ Int32 *med3func(Int32 *a, Int32 *b, Int32 *c, UChar *text_depth)
{   int va, vb, vc;
    if ((va=ptr2char(a)) == (vb=ptr2char(b)))
        return a;
    if ((vc=ptr2char(c)) == va || vc == vb)
        return c;       
    return va < vb ?
          (vb < vc ? b : (va < vc ? c : a ) )
        : (vb > vc ? b : (va < vc ? a : c ) );
}
#define med3(a, b, c) med3func(a, b, c, text_depth)


static void shallow_mkq(Int32 *a, int n, UChar *text_depth)
{
  int d, r, partval;
  Int32 *pa, *pb, *pc, *pd, *pl, *pm, *pn, t;
  UChar *next_depth;

  if (n < Mk_qs_thresh) {
    shallow_inssort_lcp(a, n, text_depth);
    return;
  }

 repeat:
  pl = a;
  pm = a + (n/2);
  pn = a + (n-1);
  if (n > 30) { // On big arrays, pseudomedian of 9
    d = (n/8);
    pl = med3(pl, pl+d, pl+2*d);
    pm = med3(pm-d, pm, pm+d);
    pn = med3(pn-2*d, pn-d, pn);
  }
  pm = med3(pl, pm, pn);
  swap2(a, pm);
  partval = ptr2char(a);
  pa = pb = a + 1;
  pc = pd = a + n-1;
  // -------- partition -----------------
  for (;;) {
    while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) {
      if (r == 0) { swap2(pa, pb); pa++; }
      pb++;
    }
    while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) {
      if (r == 0) { swap2(pc, pd); pd--; }
      pc--;
    }
    if (pb > pc) break;
    swap2(pb, pc);
    pb++;
    pc--;
  }

#if UNROLL
  if(pa>pd) {
    if( (next_depth = text_depth+1) >= Shallow_text_limit) {
      helped_sort(a, n, next_depth-Text);
      return;
    }
    else {
      text_depth = next_depth;
      goto repeat;
    }
  }
#endif
  pn = a + n;
  r = min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
  r = min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
  if ((r = pb-pa) > 1)
    shallow_mkq(a, r, text_depth);
  if( (next_depth = text_depth+1) < Shallow_text_limit)
    shallow_mkq(a + r, pa-pd+n-1, next_depth);
  else 
    helped_sort(a + r, pa-pd+n-1, next_depth-Text);
  if ((r = pd-pc) > 1)
    shallow_mkq(a + n-r, r, text_depth);
}



#define ptr2char16(i) (getword16(*(i) + text_depth))
#define getword16(s) ((unsigned)((*(s) << 8) | *((s)+1)))

#if 0
__inline__ Int32 *med3func16(Int32 *a, Int32 *b, Int32 *c, UChar *text_depth)
{   int va, vb, vc;
    if ((va=ptr2char16(a)) == (vb=ptr2char16(b)))
        return a;
    if ((vc=ptr2char16(c)) == va || vc == vb)
        return c;       
    return va < vb ?
          (vb < vc ? b : (va < vc ? c : a ) )
        : (vb > vc ? b : (va < vc ? a : c ) );
}
#define med3_16(a, b, c) med3func16(a, b, c, text_depth)
#endif

static void shallow_mkq16(Int32 *a, int n, UChar *text_depth)
{
  int d, r, partval;
  Int32 *pa, *pb, *pc, *pd, *pl, *pm, *pn, t;
  UChar *next_depth;

  if (n < Mk_qs_thresh) {
    shallow_inssort_lcp(a, n, text_depth);
    return;
  }

 repeat:
  pl = a;
  pm = a + (n/2);
  pn = a + (n-1);
  if (n > 30) { // On big arrays, pseudomedian of 9
    d = (n/8);
    pl = med3(pl, pl+d, pl+2*d);
    pm = med3(pm-d, pm, pm+d);
    pn = med3(pn-2*d, pn-d, pn);
  }
  pm = med3(pl, pm, pn);
  swap2(a, pm);
  partval = ptr2char16(a);
  pa = pb = a + 1;
  pc = pd = a + n-1;
  for (;;) {
    while (pb <= pc && (r = ptr2char16(pb)-partval) <= 0) {
      if (r == 0) { swap2(pa, pb); pa++; }
      pb++;
    }
    while (pb <= pc && (r = ptr2char16(pc)-partval) >= 0) {
      if (r == 0) { swap2(pc, pd); pd--; }
      pc--;
    }
    if (pb > pc) break;
    swap2(pb, pc);
    pb++;
    pc--;
  }
#if UNROLL
  if(pa>pd) {
    if( (next_depth = text_depth+2) >= Shallow_text_limit) {
      helped_sort(a, n, next_depth-Text);
      return;
    }
    else {
      text_depth = next_depth;
      goto repeat;
    }
  }
#endif
  pn = a + n;
  r = min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
  r = min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
  if ((r = pb-pa) > 1)
    shallow_mkq16(a, r, text_depth);
  if( (next_depth = text_depth+2) < Shallow_text_limit)
    shallow_mkq16(a + r, pa-pd+n-1, next_depth);
  else 
    helped_sort(a + r, pa-pd+n-1, next_depth-Text);
  if ((r = pd-pc) > 1)
    shallow_mkq16(a + n-r, r, text_depth);
}


#define ptr2char32(i) (getword32(*(i) + text_depth))
#define getword32(s) ((unsigned)( (*(s) << 24) | ((*((s)+1)) << 16) \
                                  | ((*((s)+2)) << 8) | (*((s)+3)) ))
static void shallow_mkq32(Int32 *a, int n, UChar *text_depth)
{
  UInt32 partval, val;
  Int32 *pa, *pb, *pc, *pd, *pl, *pm, *pn, t, d, r;
  UChar *next_depth;

  if (n < Mk_qs_thresh) {
    shallow_inssort_lcp(a, n, text_depth);
    return;
  }

 repeat:
  pl = a;
  pm = a + (n/2);
  pn = a + (n-1);
  if (n > 30) { // On big arrays, pseudomedian of 9
    d = (n/8);
    pl = med3(pl, pl+d, pl+2*d);
    pm = med3(pm-d, pm, pm+d);
    pn = med3(pn-2*d, pn-d, pn);
  }
  pm = med3(pl, pm, pn);
  swap2(a, pm);
  partval = ptr2char32(a);
  pa = pb = a + 1;
  pc = pd = a + n-1;
  for (;;) {
    while (pb <= pc &&  (val=ptr2char32(pb)) <=  partval) {
      if (val == partval) { swap2(pa, pb); pa++; }
      pb++;
    }
    while (pb <= pc && (val=ptr2char32(pc)) >= partval) {
      if (val == partval) { swap2(pc, pd); pd--; }
      pc--;
    }
    if (pb > pc) break;
    swap2(pb, pc);
    pb++;
    pc--;
  }
#if UNROLL
  if(pa>pd) {
    if( (next_depth = text_depth+4) >= Shallow_text_limit) {
      helped_sort(a, n, next_depth-Text);
      return;
    }
    else {
      text_depth = next_depth;
      goto repeat;
    }
  }
#endif
  pn = a + n;
  r = min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
  r = min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
  if ((r = pb-pa) > 1)
    shallow_mkq32(a, r, text_depth);
  if( (next_depth = text_depth+4) < Shallow_text_limit)
    shallow_mkq32(a + r, pa-pd+n-1, next_depth);
  else 
    helped_sort(a + r, pa-pd+n-1, next_depth-Text);
  if ((r = pd-pc) > 1)
    shallow_mkq32(a + n-r, r, text_depth);
}

static Int32 Cmp_left;
__inline__ 
Int32 cmp_unrolled_shallow_lcp(UChar *b1, UChar *b2)
{

  UChar c1, c2;
  assert(b1 != b2);

  do {
    // 1
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 2
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  1; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 3
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  2; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 4
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  3; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 5
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  4; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 6
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  5; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 7
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  6; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 8
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  7; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 9
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  8; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 10
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -=  9; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 11
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -= 10; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 12
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -= 11; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 13
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -= 12; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 14
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -= 13; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 15
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -= 14; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    // 16
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      Cmp_left -= 15; return ((UInt32)c1 - (UInt32)c2); }
    b1++; b2++; 
    Cmp_left -= 16;
    if(Cmp_left<=0) return 0;
  } while(1);
  return b2 - b1;
} 

static int lcp_aux[1+Max_thresh];
static int *lcp=lcp_aux+1; 
static void shallow_inssort_lcp(Int32 *a, Int32 n, UChar *text_depth)
{   
  Int32 i, j, j1, lcp_new, r, ai,lcpi;
  Int32 cmp_from_limit;
  UChar *text_depth_ai;

  lcp_aux[0] = -1;               // set lcp[-1] = -1
  for(i=0;i<n;i++) lcp[i]=0;     // I think this loop is not necessary
  cmp_from_limit = Shallow_text_limit-text_depth;

  for (i = 1; i< n ; i++) {
    ai = a[i]; lcpi = 0;
    text_depth_ai = ai + text_depth;
    j=i; j1=j-1;                  // j1 is a shorhand for j-1
    while(1) {           

      Cmp_left = cmp_from_limit-lcpi;  
      r = cmp_unrolled_shallow_lcp(lcpi+a[j1]+text_depth,lcpi+text_depth_ai);
      lcp_new = cmp_from_limit - Cmp_left;       // lcp between ai and a[j1] 
      assert(r!=0 || lcp_new>= cmp_from_limit);

      if(r<=0) {         // we have a[j-1] <= ai
	lcp[j1]=lcp_new; // ai will be written in a[j]; update lcp[j-1]
	break;
      }

      lcpi = lcp_new;                
      do {
	a[j] = a[j1];               // move down a[j-1]
        lcp[j] = lcp[j1];           // move down lcp[j-1]
	j=j1; j1--;                 // update j and j1=j-1
      } while(lcpi<lcp[j1]);        // recall that lcp[-1]=-1

      if(lcpi>lcp[j1]) break;       // ai will be written in position j

    }     // end for(j=i ...
    a[j]=ai;
    lcp[j]=lcpi;
  }       // end for(i=1 ... 

  for(i=0;i<n-1;i=j+1) {
    for(j=i; j<n ;j++)
      if(lcp[j]<cmp_from_limit) break;
    if(j-i>0) 
      helped_sort(a+i,j-i+1,Shallow_limit); 
  }
}











