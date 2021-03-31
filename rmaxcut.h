#ifndef RMAXCUT_H
#define RMAXCUT_H

#include <stdint.h>
#include <stdio.h>

typedef struct {
	char *name;
	int32_t w[2][2];
	int8_t s;
} mc_node_t;

typedef struct {
	uint64_t x;
	int32_t w;
} mc_edge_t;

typedef struct {
	uint32_t n_node, m_node;
	mc_node_t *node;
	uint32_t n_edge, m_edge;
	mc_edge_t *edge;
	uint64_t *idx; // index to find edges
	uint64_t *cc; // connected components
	void *h_name;
} mc_graph_t;

typedef struct {
	int32_t topn;
	int32_t n_perturb;
	double f_perturb;
	uint64_t seed;
} mc_opt_t;

extern int mc_verbose;

mc_graph_t *mc_read(const char *fn);
void mc_destroy(mc_graph_t *g);
void mc_print_cut(FILE *fp, const mc_graph_t *g);

void mc_opt_init(mc_opt_t *opt);
void mc_solve(const mc_opt_t *opt, mc_graph_t *g);

double mc_realtime(void);

#endif
