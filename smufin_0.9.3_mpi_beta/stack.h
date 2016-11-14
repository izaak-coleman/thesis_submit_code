#ifndef __STACK_H
#define __STACK_H

struct stack;

struct stack *stack_new(void);
void stack_delete(struct stack *);

void stack_push(struct stack *s, int value);
int stack_top(struct stack *s);
int stack_pop(struct stack *s);

#endif
