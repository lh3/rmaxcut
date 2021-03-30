#include <string.h>
#include <stdio.h>
#include "mcpriv.h"

static inline uint64_t kr_splitmix64(uint64_t x)
{
	uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

static inline double kr_drand_r(uint64_t *x)
{
    union { uint64_t i; double d; } u;
	*x = kr_splitmix64(*x);
    u.i = 0x3FFULL << 52 | (*x) >> 12;
    return u.d - 1.0;
}

static inline void ks_shuffle_uint32_t(size_t n, uint32_t a[], uint64_t *x)
{
	size_t i, j;
	for (i = n; i > 1; --i) {
		uint32_t tmp;
		j = (size_t)(kr_drand_r(x) * i);
		tmp = a[j]; a[j] = a[i-1]; a[i-1] = tmp;
	}
}

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

typedef struct {
	uint64_t x; // RNG
	uint32_t cc_off, cc_size;
	uint32_t n_cc_edge, m_cc_edge;
	uint32_t *cc_node;
	uint64_t *cc_edge;
	uint64_t *srt;
	int8_t *s, *s_tmp;
} mc_buf_t;

mc_buf_t *mc_buf_init(const mc_graph_t *g)
{
	uint32_t st, i, max_deg = 0, max_cc = 0;
	mc_buf_t *b;
	MC_CALLOC(b, 1);
	for (i = 0; i < g->n_node; ++i) {
		uint32_t n = (uint32_t)g->idx[i];
		max_deg = max_deg > n? max_deg : n;
	}
	MC_CALLOC(b->srt, max_deg);
	for (st = 0, i = 1; i <= g->n_node; ++i)
		if (i == g->n_node || g->cc[st]>>32 != g->cc[i]>>32)
			max_cc = max_cc > i - st? max_cc : i - st, st = i;
	MC_MALLOC(b->cc_node, max_cc);
	MC_CALLOC(b->s, g->n_node);
	MC_CALLOC(b->s_tmp, g->n_node);
	return b;
}

void mc_buf_destroy(mc_buf_t *b)
{
	free(b->cc_edge); free(b->cc_node);
	free(b->srt);
	free(b->s); free(b->s_tmp);
	free(b);
}

int64_t mc_score(const mc_graph_t *g, mc_buf_t *b, int32_t topn)
{
	uint32_t i;
	int64_t z = 0;
	for (i = 0; i < g->n_node; ++i) {
		uint32_t o = g->idx[i] >> 32;
		uint32_t j, n = (uint32_t)g->idx[i];
		if (topn > 0) { // only evaluate edges with large abs weight
			for (j = 0; j < n; ++j) {
				int32_t w = g->edge[o + j].w;
				w = w > 0? w : -w;
				b->srt[j] = (uint64_t)((uint32_t)-1 - w) << 32 | (o + j);
			}
			radix_sort_mc64(b->srt, b->srt + n);
			for (j = 0; j < n && j < topn; ++j) {
				const mc_edge_t *e = &g->edge[(uint32_t)b->srt[j]];
				z += g->node[e->x>>32].s * g->node[(uint32_t)e->x].s * e->w;
			}
		} else { // standard max-cut
			for (j = 0; j < n; ++j) {
				const mc_edge_t *e = &g->edge[o + j];
				z += g->node[e->x>>32].s * g->node[(uint32_t)e->x].s * e->w;
			}
		}
	}
	return z;
}

void mc_init_spin(const mc_graph_t *g, mc_buf_t *b)
{
	uint32_t i;
	b->n_cc_edge = 0;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)g->cc[b->cc_off + i];
		uint32_t o = g->idx[k] >> 32;
		uint32_t n = (uint32_t)g->idx[k], j;
		b->cc_node[i] = k;
		for (j = 0; j < n; ++j) {
			int32_t w = g->edge[o + j].w;
			w = w > 0? w : -w;
			if (b->n_cc_edge == b->m_cc_edge)
				MC_EXPAND(b->cc_edge, b->m_cc_edge);
			b->cc_edge[b->n_cc_edge++] = (uint64_t)((uint32_t)-1 - w) << 32 | (o + j);
		}
	}
	radix_sort_mc64(b->cc_edge, b->cc_edge + b->n_cc_edge);
	for (i = 0; i < b->n_cc_edge; ++i) { // from the strongest edge to the weakest
		const mc_edge_t *e = &g->edge[(uint32_t)b->cc_edge[i]];
		uint32_t n1 = e->x>>32, n2 = (uint32_t)e->x;
		if (b->s[n1] == 0 && b->s[n2] == 0) {
			b->x = kr_splitmix64(b->x);
			b->s[n1] = b->x&1? 1 : -1;
			b->s[n2] = e->w > 0? -b->s[n1] : b->s[n1];
		} else if (b->s[n1] == 0) {
			b->s[n1] = e->w > 0? -b->s[n1] : b->s[n1];
		} else if (b->s[n2] == 0) {
			b->s[n2] = -b->s[n1];
		}
	}
}

void mc_solve(mc_graph_t *g)
{
	mc_buf_t *b;
	mc_find_cc(g);
	b = mc_buf_init(g);
	mc_buf_destroy(b);
}
