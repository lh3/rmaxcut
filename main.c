#include <stdio.h>
#include <stdlib.h>
#include "ketopt.h"
#include "rmaxcut.h"

static void print_usage(FILE *fp, const mc_opt_t *opt)
{
	fprintf(fp, "Usage: rmaxcut [options] <in.txt>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -s INT      RNG seed [%lu]\n", (unsigned long)opt->seed);
}

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	mc_graph_t *g;
	mc_opt_t opt;
	int c;

	mc_realtime();
	mc_opt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "s:", 0)) >= 0) {
		if (c == 's') opt.seed = atol(o.arg);
	}
	if (o.ind == argc) {
		print_usage(stdout, &opt);
		return 1;
	}
	g = mc_read(argv[o.ind]);
	mc_solve(&opt, g);
	mc_destroy(g);
	return 0;
}
