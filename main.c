#include <stdio.h>
#include "ketopt.h"
#include "rmaxcut.h"

static void print_usage(FILE *fp)
{
	fprintf(fp, "Usage: rmaxcut [options] <in.txt>\n");
}

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	mc_graph_t *g;
	int c;
	mc_realtime();
	while ((c = ketopt(&o, argc, argv, 1, "", 0)) >= 0) {
	}
	if (o.ind == argc) {
		print_usage(stdout);
		return 1;
	}
	g = mc_read(argv[o.ind]);
	mc_solve(g);
	mc_destroy(g);
	return 0;
}
