#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"

int find_interval(const int *posArray, int size, int val)
{
	int l,r;

	l = 0;
	r = size;

	while (l<=r)
	{
		int m = l + (r - l - 1) / 2;

		if (val < posArray[m])
			r = m - 1;
		else
		{
			if (val < posArray[m+1] || m+1 == size)
			{
				return m;
			}
			
			l = m + 1;
		}
	}

	assert(0 == 1);
	return -1;
}


