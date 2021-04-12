## Introduction

Rmaxcut finds an *approximate* solution to a weighted [max-cut][max-cut]
problem via random perturbation. Each line in an input file consists of
the first nodeID, the second nodeID and an integer weight. The output gives
```txt
N  nodeID   spin     w++  w+-  w-+  w--
E  nodeID1  nodeID2  weight  spin1  spin2
```
where `spin` is either 1 or -1, indicating the partition of `nodeID`. On an
N-line, `w++` is the sum of positive weights of neighbors with positive spins;
other `w**` numbers are similar.  These numbers are mostly for debugging
purposes.

To try rmaxcut, you may acquire `x-all.txt.gz` from the [download
page][download] and run
```sh
rmaxcut -r20000 x-all.txt.gz > test.out 2> test.err
```
It will take a couple of minutes on a single thread. Increasing option `-r`,
the number of iterations, often leads to a better solution for the largest few
connected components. Rmaxcut emits good enough solutions to problems at my
hand. I haven't compared it to other more sophisticated max-cut solvers. Use
with caution.

## Algorithms

The weighted [max-cut][max-cut] problem is to find a bipartition in an
undirected graph such that the sum of edge weights between the two partitions 
are maximized. Max-cut is [NP-complete][np-comp]. In physics, solving max-cut
is equivalent to finding optimal states in an [Ising model][ising] without
external fields. We often take ideas in physics to find *approximate* solutions
to large-scale max-cut problems. This repo provides such an implementation.
It uses a fairly naive algorithm:

1. Find connected component. Do the following in each component.
2. Randomly partition nodes, or equivalently, randomly assign a spin {-1,1} to
   each node.
3. For each node, flip its spin if doing that increases the weight. Repeat
   until the total weight can't be improved this way. This is a local maximum.
4. Randomly flip 10% of spins, or flip the spins of nodes within a certain
   distance from a random node (done by [BFS][bfs]). Do step 3 again. If the
   new local maximum is better, do step 4 on the new state; otherwise, do step
   4 on the old state.
5. Repeat 3-4 for many times. Then move to the next component.

The algorithm is akin to [simulated annealing][sa] but without proper annealing.
When I have time, I need to try more sophisticated procedures in vein of
Monte Carlo methods. This may speed up heuristic search.

[max-cut]: https://en.wikipedia.org/wiki/Maximum_cut
[np-comp]: https://en.wikipedia.org/wiki/NP-completeness
[ising]: https://en.wikipedia.org/wiki/Ising_model#Connection_to_graph_maximum_cut
[sa]: https://en.wikipedia.org/wiki/Simulated_annealing
[download]: https://github.com/lh3/rmaxcut/releases/tag/data1
[bfs]: https://en.wikipedia.org/wiki/Breadth-first_search
