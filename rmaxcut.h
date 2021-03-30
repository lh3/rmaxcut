#ifndef RMAXCUT_H
#define RMAXCUT_H

#include <stdint.h>

typedef struct {
	char *name;
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
	uint64_t *idx;
	void *h_name;
} mc_graph_t;

extern int mc_verbose;

mc_graph_t *mc_read(const char *fn);
void mc_destroy(mc_graph_t *g);

#endif
