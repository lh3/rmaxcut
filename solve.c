#include <string.h>
#include <assert.h>
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
	opt->n_perturb = 20000;
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
	int64_t z[2];
} mc_pairsc_t;

typedef struct {
	uint64_t x; // RNG
	uint32_t cc_off, cc_size;
	uint32_t n_cc_edge, m_cc_edge;
	uint32_t n_sub;
	uint32_t *cc_node;
	uint64_t *cc_edge;
	uint32_t *bfs, *bfs_mark;
	mc_pairsc_t *z, *z_opt;
	int8_t *s, *s_opt;
} mc_svaux_t;

mc_svaux_t *mc_svaux_init(const mc_graph_t *g, uint64_t x)
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
	MC_MALLOC(b->bfs, g->n_node);
	MC_MALLOC(b->bfs_mark, g->n_node);
	MC_CALLOC(b->z, g->n_node);
	MC_CALLOC(b->z_opt, g->n_node);
	for (i = 0; i < g->n_node; ++i) b->bfs_mark[i] = (uint32_t)-1;
	return b;
}

void mc_svaux_destroy(mc_svaux_t *b)
{
	free(b->cc_edge); free(b->cc_node);
	free(b->s); free(b->s_opt);
	free(b->z); free(b->z_opt);
	free(b->bfs); free(b->bfs_mark);
	free(b);
}

void mc_reset_z(const mc_graph_t *g, mc_svaux_t *b)
{
	uint32_t i;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)g->cc[b->cc_off + i];
		uint32_t o = g->idx[k] >> 32;
		uint32_t j, n = (uint32_t)g->idx[k];
		b->z[k].z[0] = b->z[k].z[1] = 0;
		for (j = 0; j < n; ++j) {
			const mc_edge_t *e = &g->edge[o + j];
			uint32_t t = (uint32_t)e->x;
			if (b->s[t] > 0) b->z[k].z[0] += e->w;
			else if (b->s[t] < 0) b->z[k].z[1] += e->w;
		}
	}
}

static void mc_set_spin(const mc_graph_t *g, mc_svaux_t *b, uint32_t k, int8_t s)
{
	uint32_t o, j, n;
	int8_t s0 = b->s[k];
	if (s0 == s) return;
	o = g->idx[k] >> 32;
	n = (uint32_t)g->idx[k];
	for (j = 0; j < n; ++j) {
		const mc_edge_t *e = &g->edge[o + j];
		uint32_t t = (uint32_t)e->x;
		if (s0 != 0) b->z[t].z[(s0 < 0)] -= e->w;
		if (s  != 0) b->z[t].z[(s  < 0)] += e->w;
	}
	b->s[k] = s;
}

int64_t mc_score(const mc_graph_t *g, mc_svaux_t *b)
{
	uint32_t i;
	int64_t z = 0;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)g->cc[b->cc_off + i];
		z += -b->s[k] * (b->z[k].z[0] - b->z[k].z[1]);
	}
	return z;
}

int64_t mc_init_spin(const mc_opt_t *opt, const mc_graph_t *g, mc_svaux_t *b)
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
		}
	}
	mc_reset_z(g, b);
	return mc_score(g, b);
}

static uint32_t mc_bfs(const mc_graph_t *g, mc_svaux_t *b, uint32_t k0, int32_t bfs_round, uint32_t max_size)
{
	uint32_t i, n_bfs = 0, st, en, r;
	b->bfs[n_bfs++] = k0, b->bfs_mark[k0] = k0;
	st = 0, en = n_bfs;
	for (r = 0; r < bfs_round; ++r) {
		for (i = st; i < en; ++i) {
			uint32_t k = b->bfs[i];
			uint32_t o = g->idx[k] >> 32;
			uint32_t n = (uint32_t)g->idx[k], j;
			for (j = 0; j < n; ++j) {
				uint32_t t = (uint32_t)g->edge[o + j].x;
				if (b->bfs_mark[t] != k0)
					b->bfs[n_bfs++] = t, b->bfs_mark[t] = k0;
			}
		}
		st = en, en = n_bfs;
		if (max_size > 0 && n_bfs > max_size) break;
	}
	return n_bfs;
}

static void mc_perturb_node(const mc_opt_t *opt, const mc_graph_t *g, mc_svaux_t *b, int32_t bfs_round)
{
	uint32_t i, k, n_bfs = 0;
	k = (uint32_t)(kr_drand_r(&b->x) * b->cc_size + .499);
	k = (uint32_t)g->cc[b->cc_off + k];
	n_bfs = mc_bfs(g, b, k, bfs_round, (int32_t)(b->cc_size * opt->f_perturb));
	for (i = 0; i < n_bfs; ++i)
		mc_set_spin(g, b, b->bfs[i], -b->s[b->bfs[i]]);
}

static void mc_perturb(const mc_opt_t *opt, const mc_graph_t *g, mc_svaux_t *b)
{
	uint32_t i;
	for (i = 0; i < b->cc_size; ++i) {
		uint32_t k = (uint32_t)g->cc[b->cc_off + i];
		double y;
		y = kr_drand_r(&b->x);
		if (y < opt->f_perturb)
			mc_set_spin(g, b, k, -b->s[k]);
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
			int8_t s;
			if (b->z[k].z[0] == b->z[k].z[1]) continue;
			s = b->z[k].z[0] > b->z[k].z[1]? -1 : 1;
			if (b->s[k] != s) {
				mc_set_spin(g, b, k, s);
				++n_flip;
			}
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
	for (j = 0; j < b->cc_size; ++j) {
		b->s_opt[b->cc_node[j]] = b->s[b->cc_node[j]];
		b->z_opt[b->cc_node[j]] = b->z[b->cc_node[j]];
	}
	for (k = 0; k < opt->n_perturb; ++k) {
		if (k&1) mc_perturb(opt, g, b);
		else mc_perturb_node(opt, g, b, 3);
		sc = mc_optimize_local(opt, g, b, &n_iter);
		if (sc > sc_opt) {
			for (j = 0; j < b->cc_size; ++j) {
				b->s_opt[b->cc_node[j]] = b->s[b->cc_node[j]];
				b->z_opt[b->cc_node[j]] = b->z[b->cc_node[j]];
			}
			sc_opt = sc;
		} else {
			for (j = 0; j < b->cc_size; ++j) {
				b->s[b->cc_node[j]] = b->s_opt[b->cc_node[j]];
				b->z[b->cc_node[j]] = b->z_opt[b->cc_node[j]];
			}
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
	b = mc_svaux_init(g, opt->seed);
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
