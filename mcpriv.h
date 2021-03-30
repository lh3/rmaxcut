#ifndef MCPRIV_H
#define MCPRIV_H

#include <stdlib.h>
#include "rmaxcut.h"

#define MC_MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define MC_CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MC_REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))
#define MC_BZERO(ptr, len) memset((ptr), 0, (len) * sizeof(*(ptr)))
#define MC_EXPAND(a, m) do { \
		(m) = (m)? (m) + ((m)>>1) : 16; \
		MC_REALLOC((a), (m)); \
	} while (0)

void radix_sort_mc64(uint64_t*, uint64_t*);

#endif
