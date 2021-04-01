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

static void ks_shuffle_uint32_t(size_t n, uint32_t a[], uint64_t *x)
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
	opt->topn_pos = opt->topn_neg = -1;
	opt->n_perturb = 2000;
	opt->f_perturb = 0.1;
	opt->max_iter = 1000;
	opt->seed = 11;
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
	uint32_t n_sub;
	uint32_t *cc_node;
	uint64_t *cc_edge;
	uint64_t *idx;
	mc_edge_t *sub;
	int8_t *s, *s_opt;
} mc_svaux_t;

mc_edge_t *mc_topn(const mc_graph_t *g, int32_t topn_pos, int32_t topn_neg, uint32_t *n_edge_)
{
	uint32_t i, max_deg = 0, k;
	uint64_t *srt_pos, *srt_neg;
	mc_edge_t *edge;

	MC_MALLOC(edge, g->n_edge);
	if (topn_pos <= 0 && topn_neg <= 0) {
		memcpy(edge, g->edge, g->n_edge * sizeof(mc_edge_t));
		*n_edge_ = g->n_edge;
		return edge;
	}

	*n_edge_ = 0;
	for (i = 0; i < g->n_node; ++i) {
		uint32_t n = (uint32_t)g->idx[i];
		max_deg = max_deg > n? max_deg : n;
	}
	MC_MALLOC(srt_pos, max_deg);
	MC_MALLOC(srt_neg, max_deg);
	for (i = 0, k = 0; i < g->n_node; ++i) {
		uint32_t o = g->idx[i] >> 32;
		uint32_t j, n = (uint32_t)g->idx[i];
		uint32_t n_pos = 0, n_neg = 0;
		for (j = 0; j < n; ++j) {
			int32_t w = g->edge[o + j].w;
			if (w > 0)
				srt_pos[n_pos++] = (uint64_t)((uint32_t)-1 - w) << 32 | (o + j);
			else if (w < 0)
				srt_neg[n_neg++] = (uint64_t)((uint32_t)-1 + w) << 32 | (o + j);
		}
		radix_sort_mc64(srt_pos, srt_pos + n_pos);
		radix_sort_mc64(srt_neg, srt_neg + n_neg);
		for (j = 0; j < n_pos && j < topn_pos; ++j)
			edge[k++] = g->edge[(uint32_t)srt_pos[j]];
		for (j = 0; j < n_neg && j < topn_neg; ++j)
			edge[k++] = g->edge[(uint32_t)srt_neg[j]];
	}
	free(srt_neg); free(srt_pos);
	*n_edge_ = k;
	if (mc_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] reduced to %u edges\n", __func__, mc_realtime(), k);
	return edge;
}

mc_svaux_t *mc_svaux_init(const mc_graph_t *g, uint64_t x, int32_t topn_pos, int32_t topn_neg)
{
	uint32_t st, i, max_cc = 0;
	mc_svaux_t *b;
	MC_CALLOC(b, 1);
	b->x = x;
	for (st = 0, i = 1; i <= g->n_node; ++i)
		if (i == g->n_node || g->cc[st]>>32 != g->cc[i]>>32)
			max_cc = max_cc > i - st? max_cc : i - st, st = i;
	MC_MALLOC(b->cc_node, max_cc);
	MC_CALLOC(b->s, g->n_node);
	MC_CALLOC(b->s_opt, g->n_node);
	b->sub = mc_topn(g, topn_pos, topn_neg, &b->n_sub);
	b->idx = mc_index_core(b->sub, b->n_sub, g->n_node);
	return b;
}

void mc_svaux_destroy(mc_svaux_t *b)
{
	free(b->cc_edge); free(b->cc_node);
	free(b->s); free(b->s_opt);
	free(b->idx); free(b->sub);
	free(b);
}

int64_t mc_score(const mc_graph_t *g, mc_svaux_t *b)
{
	uint32_t i;
	int64_t z = 0;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)g->cc[b->cc_off + i];
		uint32_t o = b->idx[k] >> 32;
		uint32_t j, n = (uint32_t)b->idx[k];
		for (j = 0; j < n; ++j) {
			const mc_edge_t *e = &b->sub[o + j];
			z += -b->s[e->x>>32] * b->s[(uint32_t)e->x] * e->w;
		}
	}
	return z;
}

int64_t mc_init_spin(const mc_opt_t *opt, const mc_graph_t *g, mc_svaux_t *b)
{
	uint32_t i;
	b->n_cc_edge = 0;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)g->cc[b->cc_off + i];
		uint32_t o = b->idx[k] >> 32;
		uint32_t n = (uint32_t)b->idx[k], j;
		b->cc_node[i] = k;
		for (j = 0; j < n; ++j) {
			int32_t w = b->sub[o + j].w;
			w = w > 0? w : -w;
			if (b->n_cc_edge == b->m_cc_edge)
				MC_EXPAND(b->cc_edge, b->m_cc_edge);
			b->cc_edge[b->n_cc_edge++] = (uint64_t)((uint32_t)-1 - w) << 32 | (o + j);
		}
	}
	radix_sort_mc64(b->cc_edge, b->cc_edge + b->n_cc_edge);
	for (i = 0; i < b->n_cc_edge; ++i) { // from the strongest edge to the weakest
		const mc_edge_t *e = &b->sub[(uint32_t)b->cc_edge[i]];
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
	return mc_score(g, b);
}

static void mc_perturb(const mc_opt_t *opt, const mc_graph_t *g, mc_svaux_t *b)
{
	uint32_t i;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)g->cc[b->cc_off + i];
		double y;
		y = kr_drand_r(&b->x);
		if (y < opt->f_perturb)
			b->s[k] = -b->s[k];
	}
}

static int64_t mc_optimize_local(const mc_opt_t *opt, const mc_graph_t *g, mc_svaux_t *b, uint32_t *n_iter)
{
	int32_t n_iter_local = 0;
	while (n_iter_local < opt->max_iter) {
		uint32_t i, n_flip = 0;
		++(*n_iter);
		ks_shuffle_uint32_t(b->cc_size, b->cc_node, &b->x);
		for (i = 0; i < b->cc_size; ++i) {
			uint32_t k = b->cc_node[i];
			uint32_t o = b->idx[k] >> 32;
			uint32_t n = (uint32_t)b->idx[k], j;
			int64_t z[2];
			int8_t s;
			z[0] = z[1] = 0;
			for (j = 0; j < n; ++j) {
				const mc_edge_t *e = &b->sub[o + j];
				if (b->s[(uint32_t)e->x] > 0) z[0] += e->w;
				else if (b->s[(uint32_t)e->x] < 0) z[1] += e->w;
			}
			if (z[0] == z[1]) continue;
			s = z[0] > z[1]? -1 : 1;
			if (b->s[k] != s)
				b->s[k] = s, ++n_flip;
		}
		++n_iter_local;
		if (n_flip == 0) break;
	}
	return mc_score(g, b);
}

uint32_t mc_solve_cc(const mc_opt_t *opt, const mc_graph_t *g, mc_svaux_t *b, uint32_t cc_off, uint32_t cc_size)
{
	uint32_t j, k, n_iter = 0;
	int64_t sc_ori, sc_opt = -(1<<30), sc;

	b->cc_off = cc_off, b->cc_size = cc_size;
	if (b->cc_size < 2) return 0;

	// first guess
	sc_ori = mc_init_spin(opt, g, b);
	if (b->cc_size == 2) return 0;

	// optimize
	sc_opt = mc_optimize_local(opt, g, b, &n_iter);
	for (j = 0; j < b->cc_size; ++j)
		b->s_opt[b->cc_node[j]] = b->s[b->cc_node[j]];
	for (k = 0; k < opt->n_perturb; ++k) {
		mc_perturb(opt, g, b);
		sc = mc_optimize_local(opt, g, b, &n_iter);
		if (sc > sc_opt) {
			for (j = 0; j < b->cc_size; ++j)
				b->s_opt[b->cc_node[j]] = b->s[b->cc_node[j]];
			sc_opt = sc;
		} else {
			for (j = 0; j < b->cc_size; ++j)
				b->s[b->cc_node[j]] = b->s_opt[b->cc_node[j]];
		}
	}
	for (j = 0; j < b->cc_size; ++j)
		b->s[b->cc_node[j]] = b->s_opt[b->cc_node[j]];
	fprintf(stderr, "[%s::%.3f] group:%d, size:%d, #edges:%d, #iter:%d, sc_ori:%ld, sc_opt:%ld\n", __func__, mc_realtime(),
			(uint32_t)(g->cc[b->cc_off]>>32), b->cc_size, b->n_cc_edge, n_iter, (long)sc_ori, (long)sc_opt);
	return n_iter;
}

void mc_write_info(mc_graph_t *g, mc_svaux_t *b)
{
	uint32_t i;
	for (i = 0; i < g->n_node; ++i) {
		mc_node_t *t = &g->node[i];
		uint32_t o = g->idx[i] >> 32;
		uint32_t j, n = (uint32_t)g->idx[i];
		t->s = b->s[i];
		t->w[0][0] = t->w[0][1] = t->w[1][0] = t->w[1][1] = 0;
		for (j = 0; j < n; ++j) {
			mc_edge_t *e = &g->edge[o + j];
			if (b->s[(uint32_t)e->x] > 0) {
				if (e->w > 0) t->w[0][0] += e->w;
				else if (e->w < 0) t->w[0][1] += e->w;
			} else if (b->s[(uint32_t)e->x] < 0) {
				if (e->w > 0) t->w[1][0] += e->w;
				else if (e->w < 0) t->w[1][1] += e->w;
			}
		}
	}
}

void mc_solve(const mc_opt_t *opt, mc_graph_t *g)
{
	uint32_t st, i;
	mc_svaux_t *b;
	mc_find_cc(g);
	b = mc_svaux_init(g, opt->seed, opt->topn_pos, opt->topn_neg);
	for (st = 0, i = 1; i <= g->n_node; ++i) {
		if (i == g->n_node || g->cc[st]>>32 != g->cc[i]>>32) {
			mc_solve_cc(opt, g, b, st, i - st);
			st = i;
		}
	}
	mc_write_info(g, b);
	mc_svaux_destroy(b);
}

void mc_print_cut(FILE *fp, const mc_graph_t *g)
{
	uint32_t i;
	for (i = 0; i < g->n_node; ++i) {
		const mc_node_t *t = &g->node[i];
		fprintf(fp, "N\t%s\t%d\t%d\t%d\t%d\t%d\n", t->name, t->s,
				t->w[0][0], t->w[0][1], t->w[1][0], t->w[1][1]);
	}
	for (i = 0; i < g->n_edge; ++i) {
		const mc_edge_t *e = &g->edge[i];
		const mc_node_t *t1 = &g->node[e->x>>32];
		const mc_node_t *t2 = &g->node[(uint32_t)e->x];
		fprintf(fp, "E\t%s\t%s\t%d\t%d\t%d\n", t1->name, t2->name, e->w, t1->s, t2->s);
	}
}
