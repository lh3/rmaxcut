#include <zlib.h>
#include <stdio.h>
#include "mcpriv.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

#include "khashl.h"
KHASHL_MAP_INIT(, mc_s2i_t, mc_s2i, const char*, int32_t, kh_hash_str, kh_eq_str)

#include "ksort.h"
#define mc_edge_key(e) ((e).x)
KRADIX_SORT_INIT(mce, mc_edge_t, mc_edge_key, 8)

int mc_verbose = 3;

char *mc_strdup(const char *src)
{
	int32_t len;
	char *dst;
	len = strlen(src);
	MC_MALLOC(dst, len + 1);
	memcpy(dst, src, len + 1);
	return dst;
}

uint32_t mc_node_add(mc_graph_t *g, const char *name)
{
	mc_s2i_t *h = (mc_s2i_t*)g->h_name;
	int absent;
	khint_t k;
	k = mc_s2i_put(h, name, &absent);
	if (absent) {
		mc_node_t *t;
		if (g->n_node == g->m_node) MC_EXPAND(g->node, g->m_node);
		t = &g->node[g->n_node++];
		kh_key(h, k) = t->name = mc_strdup(name);
		kh_val(h, k) = g->n_node - 1;
	}
	return kh_val(h, k);
}

void mc_sort(mc_graph_t *g)
{
	radix_sort_mce(g->edge, g->edge + g->n_edge);
}

void mc_merge(mc_graph_t *g) // MUST BE sorted
{
	uint32_t i, j, k, st;
	for (st = 0, i = 1, k = 0; i <= g->n_edge; ++i) {
		if (i == g->n_edge || g->edge[i].x != g->edge[st].x) {
			int64_t w = 0;
			if (i - st > 1) {
				for (j = st; j < i; ++j)
					w += g->edge[j].w;
				w = (int64_t)((double)w / (i - st) + .499);
			} else w = g->edge[j].w;
			g->edge[k] = g->edge[st];
			g->edge[k++].w = w;
			st = i;
		}
	}
	g->n_edge = g->m_edge = k;
	MC_REALLOC(g->edge, g->m_edge);
}

void mc_index(mc_graph_t *g) // MUST BE sorted
{
	uint32_t i, st;
	free(g->idx);
	MC_CALLOC(g->idx, g->n_node);
	for (st = 0, i = 1; i <= g->n_edge; ++i) {
		if (i == g->n_edge || g->edge[i].x>>32 != g->edge[st].x>>32) {
			uint32_t t = g->edge[st].x>>32;
			g->idx[t] = (uint64_t)st << 32 | (i - st);
			st = i;
		}
	}
}

mc_graph_t *mc_read(const char *fn)
{
	mc_graph_t *g;
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int dret;

	fp = gzopen(fn, "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);

	MC_CALLOC(g, 1);
	g->h_name = mc_s2i_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		char *p, *q, *u1 = 0, *u2 = 0;
		int32_t i, w;
		for (i = 0, p = q = str.s;; ++p) {
			if (*p == '\t' || *p == 0) {
				int c = *p;
				*p = 0;
				if (i == 0) u1 = q;
				else if (i == 1) u2 = q;
				else if (i == 2) w = atoi(q);
				++i, q = p + 1;
				if (c == 0 || i == 3) break;
			}
		}
		if (i == 3) {
			uint32_t nid1, nid2;
			nid1 = mc_node_add(g, u1);
			nid2 = mc_node_add(g, u2);
			if (g->n_edge == g->m_edge) MC_EXPAND(g->edge, g->m_edge);
			g->edge[g->n_edge].x = (uint64_t)nid1 << 32 | nid2;
			g->edge[g->n_edge++].w = w;
			if (g->n_edge == g->m_edge) MC_EXPAND(g->edge, g->m_edge);
			g->edge[g->n_edge].x = (uint64_t)nid2 << 32 | nid1;
			g->edge[g->n_edge++].w = w;
		}
	}
	ks_destroy(ks);
	gzclose(fp);
	mc_sort(g);
	mc_merge(g);
	if (mc_verbose >= 3)
		fprintf(stderr, "[%s] read %u nodes and %u edges\n", __func__, g->n_node, g->n_edge);
	mc_index(g);
	return g;
}

void mc_destroy(mc_graph_t *g)
{
	free(g->node); free(g->edge); free(g->idx);
	mc_s2i_destroy(g->h_name);
	free(g);
}
