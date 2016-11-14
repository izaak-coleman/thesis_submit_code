#ifndef SUFFIX__GMEM_H
#define SUFFIX__GMEM_H

void *gsuffix_malloc(int size);
void gsuffix_free(void *);
void *gsuffix_realloc(void *mem, int size);


unsigned int gsuffix_memallocated(void);
unsigned int gsuffix_maxmem(void);
void gsuffix_resetmaxmem(void);

struct item_pool;

struct item_pool *item_pool_create(int item_size);
void item_pool_delete(struct item_pool *ip);
void *item_pool_alloc(struct item_pool *pool);
void item_pool_free(struct item_pool *pool, void *item);

struct freeable_item_pool;

struct freeable_item_pool *freeable_item_pool_create(int freeable_item_size);
void freeable_item_pool_delete(struct freeable_item_pool *ip);
void *freeable_item_pool_alloc(struct freeable_item_pool *pool);
void freeable_item_pool_free(struct freeable_item_pool *pool, void *item);

#endif
