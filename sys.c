#include <sys/resource.h>
#include <sys/time.h>
#include "mcpriv.h"

double mc_realtime0 = -1.0;

double mc_cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long mc_peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

static double mc_realtime_core(void)
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

double mc_realtime(void)
{
	if (mc_realtime0 < 0) mc_realtime0 = mc_realtime_core();
	return mc_realtime_core() - mc_realtime0;
}
