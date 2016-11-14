#include <assert.h>
#include <stdlib.h>

#include "gmem.h"
#include "stack.h"

struct stack
{
	int top;

	int *values;
	
	int num_values;
};

struct stack *stack_new(void)
{
	struct stack *stk = gsuffix_malloc(sizeof(struct stack));

	stk->values = gsuffix_malloc(10*sizeof(stk->values[0]));
	stk->num_values = 10;
	stk->top = -1;

	return stk;
}

void stack_delete(struct stack *stk)
{
	gsuffix_free(stk->values);
	stk->values = NULL;
	gsuffix_free(stk);
}

void stack_push(struct stack *stk, int value)
{
	assert(stk->values != NULL);
	stk->top++;

	if (stk->top == stk->num_values)
	{
		stk->values = gsuffix_realloc(stk->values, stk->num_values * sizeof(stk->values[0])*2);
		if (!stk->values) exit(20);
		stk->num_values *= 2;
	}

	stk->values[stk->top] = value;
}

int stack_top(struct stack *stk)
{
	assert(stk->values != NULL);
	assert(stk->top < stk->num_values);
	assert(stk->top >= 0);

	return stk->values[stk->top];
}

int stack_isempty(struct stack *stk)
{
	return stk->top == -1;
}

int stack_pop(struct stack *stk)
{
	assert(stk->values != NULL);
	assert(stk->top < stk->num_values);
	assert(stk->top >= 0);

	return stk->values[stk->top--];
}

#ifdef TESTME

int main(void)
{
	int i,val;

	struct stack *stk;

	stk = stack_new();

	printf("Output should be: 1,1,5,4,3,4,3,2\n");
	printf("Output is:        ");

	stack_push(stk,1);
	val = stack_top(stk);
	printf("%d,",val);

	val = stack_pop(stk);
	printf("%d,",val);
	
	stack_push(stk,2);
	stack_push(stk,3);
	stack_push(stk,4);
	stack_push(stk,3);
	stack_push(stk,4);
	stack_push(stk,5);

	val = stack_pop(stk);
	printf("%d,",val);
	val = stack_pop(stk);
	printf("%d,",val);
	val = stack_pop(stk);
	printf("%d,",val);
	val = stack_pop(stk);
	printf("%d,",val);
	val = stack_pop(stk);
	printf("%d,",val);
	val = stack_pop(stk);
	printf("%d",val);

	printf("\n");

	printf("Stack is empty? %d (should be 1)\n",stack_isempty(stk));

/* must be even */
#define NUM_ELEMENTS 10

	printf("Pushing %d elements\n",NUM_ELEMENTS);

	for (i=0;i<NUM_ELEMENTS;i++)
		stack_push(stk,i);

	printf("Now popping %d elements\n",NUM_ELEMENTS/2);
	for (;i>NUM_ELEMENTS/2;i--)
	{
		int v;

		if ((v = stack_pop(stk)) != i - 1)
		{
			printf("Expected %d, but got %d\n",i,v);
		}
	}

	printf("Stack is empty? %d (should be 0)\n",stack_isempty(stk));

	printf("Poping rest (should be %d elements)\n",NUM_ELEMENTS/2);
	i = 0;

	while (!stack_isempty(stk))
	{
		stack_pop(stk);
		i++;
	}

	printf("Popped %d elememts\n",i);

	printf("Stack is empty? %d (should be 1)\n",stack_isempty(stk));

	stack_delete(stk);
}

#endif
