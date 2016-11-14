#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gmem.h"


#define TRACKMEM 0

#if TRACKMEM

static int memallocated;
static int maxmem;

void *gsuffix_malloc(int size)
{
	int *mem = malloc(size+sizeof(int));
	if (!mem)
	{
		printf("Could not allocate %d bytes of memory! Aborting\n",size+4);
		exit(20);
	}
	mem[0] = size;
	memallocated += size;
	if (memallocated > maxmem) maxmem = memallocated;
	return mem+1;
}

void gsuffix_free(void *mem)
{
	int *m = (int*)mem;
	if (!m) return;

	memallocated -= m[-1];
	free(m-1);
}

void *gsuffix_realloc(void *mem, int size)
{
	int *m = (int*)mem;
	m--;
	memallocated -= m[0];
	m = realloc(m,size+4);
	if (!m)
	{
		printf("Could not reallocate %d bytes of memory! Aborting\n",size+4);
		exit(20);
	}
	memallocated += size;
	m[0] = size;
	if (memallocated > maxmem) maxmem = memallocated;
	return m+1;
}

unsigned int gsuffix_memallocated(void)
{
	return memallocated;
}

unsigned int gsuffix_maxmem(void)
{
	return maxmem;
}

void gsuffix_resetmaxmem(void)
{
	maxmem = memallocated;
}

#else   

void *gsuffix_malloc(int size)
{
	void *mem = malloc(size);
	if (!mem)
	{
		printf("Could not allocate %d bytes of memory! Aborting\n",size+4);
		exit(20);
	}
	return mem;
}

void gsuffix_free(void *mem)
{
	free(mem);
}

void *gsuffix_realloc(void *mem, int size)
{
	mem = realloc(mem,size);
	if (!mem)
	{
		printf("Could not reallocate %d bytes of memory! Aborting\n",size+4);
		exit(20);
	}
	return mem;
}

unsigned int gsuffix_memallocated(void)
{

	fprintf(stderr,"Need to recompile gsuffix library with TRACKMEM defined "
		" in order to use memory tracking functions.\n");
	exit(1);
	return -1;
}

unsigned int gsuffix_maxmem(void)
{
	fprintf(stderr,"Need to recompile gsuffix library with TRACKMEM defined "
		" in order to use memory tracking functions.\n");
	exit(1);
	return -1;
}

void gsuffix_resetmaxmem(void)
{
	fprintf(stderr,"Need to recompile gsuffix library with TRACKMEM defined "
		" in order to use memory tracking functions.\n");
	exit(1);
}



#endif




struct mypage
{
	struct mypage *next;	/* embedded node structure */

	int size;				/* number of bytes covered by the page */
	char *mem;				/* address of the page covered by the page */

	char *next_item;		/* pointer to next item */
	int items_left;			/* how many items are left within this page? */
};

struct item_pool
{
	struct mypage *head;
	int item_size;
};

struct item_pool *item_pool_create(int item_size)
{
	struct item_pool *ip = gsuffix_malloc(sizeof(*ip));

	ip->head = NULL;
	ip->item_size = (item_size + 3)/4*4;

	return ip;
}

void item_pool_delete(struct item_pool *ip)
{
	struct mypage *p = ip->head;
	while (p)
	{
		struct mypage *next;

		next = p->next;
		gsuffix_free(p->mem);
		gsuffix_free(p);
		p = next;
	}
	gsuffix_free(ip);
}


static struct mypage *item_pool_alloc_page(int item_size, int page_size)
{
	struct mypage *p;

	p = gsuffix_malloc(sizeof(*p));
	p->mem = p->next_item = gsuffix_malloc(page_size);
	p->size = page_size;
	p->next = NULL;
	p->items_left = page_size / item_size;

	return p;
}

void *item_pool_alloc(struct item_pool *pool)
{
	char *item;

	if (!pool->head || pool->head->items_left == 0)
	{
		struct mypage *p = item_pool_alloc_page(pool->item_size,32768);
		p->next = pool->head;
		pool->head = p;
	}

	item = pool->head->next_item;
	pool->head->next_item += pool->item_size;
	pool->head->items_left--;
	return item;
}

struct mypage2
{
	struct mypage2 *next;	/* embedded node structure */

	int size;				/* number of bytes covered by the page */
	void *mem;				/* address of the page covered by the page */

	int nitems;
};

struct item
{
	struct item *next;
};

struct freeable_item_pool
{
	struct mypage2 *page_head;
	struct item *items_head;
	int item_size;
};

struct freeable_item_pool *freeable_item_pool_create(int item_size)
{
	struct freeable_item_pool *ip = gsuffix_malloc(sizeof(*ip));

	ip->page_head = NULL;
	ip->items_head = NULL;
	ip->item_size = (item_size + 7)/4*4;

	return ip;
}

void freeable_item_pool_delete(struct freeable_item_pool *ip)
{
	struct mypage2 *p = ip->page_head;
	while (p)
	{
		struct mypage2 *next;

		next = p->next;
		gsuffix_free(p->mem);
		gsuffix_free(p);
		p = next;
	}
	gsuffix_free(ip);
}

static struct mypage2 *freeable_item_pool_alloc_page(int item_size, int page_size)
{
	struct mypage2 *p;

	p = gsuffix_malloc(sizeof(*p));
	p->mem = gsuffix_malloc(page_size);
	p->size = page_size;
	p->next = NULL;
	p->nitems = page_size / item_size;
	return p;
}

void *freeable_item_pool_alloc(struct freeable_item_pool *pool)
{
	struct item *item;

	if (!pool->items_head)
	{
		struct mypage2 *p;
		int i;

		p = freeable_item_pool_alloc_page(pool->item_size,32768);
		item = pool->items_head = p->mem;

		for (i=0;i < p->nitems - 1;i++)
		{
			struct item *next_item = (struct item*)(((unsigned char*)item) + pool->item_size);
			item->next = next_item;
			item = next_item;
		}
		item->next = NULL;

		p->next = pool->page_head;
		pool->page_head = p;
	}

	item = pool->items_head;
	pool->items_head = item->next;

	return ((unsigned char*)item)+4;
}

void freeable_item_pool_free(struct freeable_item_pool *pool, void *mem)
{
	struct item *item = (struct item*)(((unsigned char*)mem) - 4);
	item->next = pool->items_head;
	pool->items_head = item;
}
