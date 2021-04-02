#include <stdio.h>
#include <stdlib.h>
#include "ketopt.h"
#include "rmaxcut.h"

static void print_usage(FILE *fp, const mc_opt_t *opt)
{
	fprintf(fp, "Usage: rmaxcut [options] <in.txt>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -s INT      RNG seed [%lu]\n", (unsigned long)opt->seed);
	fprintf(fp, "  -r INT      rounds of perturbation [%d]\n", opt->n_perturb);
	fprintf(fp, "  -f FLOAT    fraction to flip for perturbation [%.3g]\n", opt->f_perturb);
}

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	mc_graph_t *g;
	mc_opt_t opt;
	int c;

	mc_realtime();
	mc_opt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "s:p:f:r:", 0)) >= 0) {
		if (c == 's') opt.seed = atol(o.arg);
		else if (c == 'r') opt.n_perturb = atoi(o.arg);
		else if (c == 'f') opt.f_perturb = atof(o.arg);
	}
	if (o.ind == argc) {
		print_usage(stdout, &opt);
		return 1;
	}
	g = mc_read(argv[o.ind]);
	mc_solve(&opt, g);
	mc_print_cut(stdout, g);
	mc_destroy(g);
	if (mc_verbose >= 3) {
		int i;
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MC_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, mc_realtime(), mc_cputime(), mc_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
