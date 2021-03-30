#include <string.h>
#include <stdio.h>
#include "mcpriv.h"

void mc_opt_init(mc_opt_t *opt)
{
	memset(opt, 0, sizeof(mc_opt_t));
}

uint64_t *mc_find_cc_core(const mc_graph_t *g)
{
	uint32_t i, x, y, *flag, *stack = 0, ns = 0, ms = 0;
	uint64_t *cc;

	MC_MALLOC(flag, g->n_node);
	for (i = 0; i < g->n_node; ++i)
		flag[i] = (uint32_t)-1;

	// connected componets
	for (i = 0; i < g->n_node; ++i) {
		if (flag[i] != (uint32_t)-1) continue;
		if (ns == ms) MC_EXPAND(stack, ms);
		stack[ns++] = i;
		while (ns > 0) {
			uint32_t k, j, n, s;
			k = stack[--ns];
			flag[k] = i;
			n = (uint32_t)g->idx[k];
			s = g->idx[k] >> 32;
			for (j = 0; j < n; ++j) {
				uint32_t t = (uint32_t)g->edge[s + j].x;
				if (flag[t] != (uint32_t)-1) continue;
				if (ns == ms) MC_EXPAND(stack, ms);
				stack[ns++] = t;
			}
		}
	}
	free(stack);

	// precalculate the size of each cc
	MC_CALLOC(cc, g->n_node);
	for (i = 0; i < g->n_node; ++i)
		cc[i] = (uint64_t)flag[i] << 32 | i;
	radix_sort_mc64(cc, cc + g->n_node);
	for (i = 1, x = y = 0; i <= g->n_node; ++i) {
		if (i == g->n_node || cc[i]>>32 != cc[x]>>32) {
			uint32_t j;
			for (j = x; j < i; ++j)
				cc[j] = (uint64_t)y << 32 | (uint32_t)cc[j];
			++y, x = i;
		}
	}
	free(flag);
	if (mc_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] found %d connected components\n", __func__, mc_realtime(), y);
	return cc;
}

void mc_find_cc(mc_graph_t *g)
{
	g->cc = mc_find_cc_core(g);
}

void mc_solve(mc_graph_t *g)
{
	mc_find_cc(g);
}
