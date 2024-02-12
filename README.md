# cloud_evolution_networks
Pipeline to generate a directed graph of the time evolution of all star-forming giant molecular clouds in a simulated galaxy, and extract their lifetimes and merger rates from this graph.

## Class to 'walk through' the network and output the lifetimes, mass evolution and merger history of clouds
Once the evolution of all clouds in a galaxy has been 'tracked' to produce the cloud evolution graph, there are two obvious ways to define the evolution of an individual cloud:
1. One cloud is a one weakly-connected component of the graph, including all its mergers and splits. I.e. clouds by definition do not interact with each other. If two objects merge, they are defined as part of the same cloud.
2. One cloud is a straight component of the graph, containing no mergers or splits. I.e. if a cloud splits in two, it generates a new cloud that starts with a lifetime of 0. If a cloud merges, it consumes a cloud and ends its lifetime.

Definition (1) is more physically-meaningful, but definition (2) is a better representation of the observable cloud population at a snapshot in time. That is, we don't know whether two clouds that appear separate in an observation will merge at a later time.

It's therefore useful to characterize cloud evolution in both of the two ways above.

### (1) Straight components or 'trajectories' of the graph
The basic implementation of this algorithm is described in the Appendix of [Jeffreson et al. 2021a]([URL](https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.1678J/abstract)https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.1678J/abstract). There is a random choice associated with the path to take through a merger or a split, so this is a Monte Carlo 'walk' through the graph, repeated to obtain a representative set of cloud evolution 'trajectories' and their associated lifetime distribution.

One Monte Carlo iteration of the algorithm used to return the lifetime and trajectories is given by XXX function.

A non-MC version that follows the most massive cloud through mergers and splits is given by YYY function.
