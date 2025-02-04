# Random-Weakly-Connected-DAG-Complex

Create a random weakly-connected Directed Acyclic Graph (DAG) on n nodes with m edges with the additional constraint of x start-nodes (source nodes) and y exit-nodes (sink nodes).

Start-nodes are nodes with no in-edges, and exit nodes are nodes with no out-edges, this also gaurantees that all nodes which are not start or exit nodes will have both in and out edges. Much like how an Erdos-Renyi graph can be created in a G(n, m) form, where n represents the number of nodes and m the number of edges, this algorithm assembles a simple weakly-connected DAG in much the same way.

This algorithm starts by creating a complete DAG on n nodes, that is the upper triangle of the adjacency matrix. Next it removes all in edges from the first x nodes, and all out edges of the last y nodes. Next the algorithm starts to randomly removes edges from that semi-complete DAG until it arrives at the final configuration. Because it is a destructive process which is randomly removing edges it can sometimes trap itself in topologies that cannot be reduced to the desired number of edges, so it also randomly adds edges if there are 2 generations in a row where an edge cannot be removed.

Note: Unlike my [previous methods](https://github.com/samrosean/Random-Weakly-Connected-DAG-Matlab) for creating random weakly-connected DAGs, this is an algorithmic method. It frequently gets stuck in incorrect topologies, so the heauristic simply makes back steps to avoid getting stuck, but these back steps are often not large enough to walk out of topological sinks. This algorithm could be optimized to avoid this, but it works for my current requirements.
