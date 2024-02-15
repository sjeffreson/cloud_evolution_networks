# cloud_evolution_networks
Pipeline to generate a directed graph of the time evolution of all star-forming giant molecular clouds in a simulated galaxy, and extract their lifetimes and merger rates from this graph.

As of February 2024, all code has been tested using Python 3.10.

This is the pipeline used for molecular cloud identification, tracking and analysis in the following works:
- [The role of galactic dynamics in shaping the physical properties of giant molecular clouds in Milky Way-like galaxies](https://ui.adsabs.harvard.edu/abs/2020MNRAS.498..385J/abstract)
- [A scaling relation for the molecular cloud lifetime in Milky Way-like galaxies](https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.1678J/abstract)
- [Momentum feedback from marginally resolved H II regions in isolated disc galaxies](https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.3470J/abstract)
- [Clouds of Theseus: long-lived molecular clouds are composed of short-lived H2 molecules](https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.7093J/abstract)

## `cloud_id.py`: Cloud identification

### Dependencies
This script uses the [h5py](https://www.h5py.org/) package to read the galaxy snapshots (hdf5 format). It also uses the [astrodendro](https://github.com/dendrograms/astrodendro/tree/stable) package. Note that Astrodendro has not been updated for 8 years and so contains some deprecated features that might show up when the package is imported, but which are easily fixed.

Note that Astrodendro is used because of the information it generates about the hierarchical structure of each cloud. This information is not currently analyzed or stored in this repo, but could easily be tracked by adding additional properties to the dictionaries defined in `cloud_id.py`. See the Astrodendro docs for all possibilities.

### Inputs
- A set of hdf5 snapshots of an isolated galaxy simulation, containing a mesh of gas particles with centroid positions, velocities etc.
- A corresponding set of 2D gas column density maps
- A corresponding set of 2D mean CO luminosity maps
In Arepo, the maps can be computed via the native ray-tracing prescription, see the [Arepo documentation](https://gitlab.mpcdf.mpg.de/vrs/arepo). These are used to identify the molecular clouds via a threshold on their CO-luminous gas column density. The properties of the corresponding clouds are then extracted from the full 3D gas cell mesh in the snapshot, via binning.

### Running
```
python3 cloud_id.py $INPUT_DIR -s $SEARCH_STRING -b $FIRST_SNAP_NO -e $LAST_SNAP_NO -w $MAP_DIAMETER -r $MAP_RESOLUTION -p $CLOUDDIR/$CLOUDSTRING -lmv $LOG_MIN_DENS -mp $MIN_PX -N $NO_CPUS
```

Where:
- Snapshots and density maps are found at `$INPUT_DIR`
- The search string for the total gas column density maps is `$SEARCH_STRING`, e.g. `density_proj_*`
- The first snapshot to analyze (time-ordered) is `$FIRST_SNAP_NO`
- The last snapshot to analyze (time-ordered) is `$LAST_SNAP_NO`
- The diameter of the maps in kiloparsec is `$MAP_DIAMETER`
- The resolution of the maps is `$MAP_RESOLUTION`
- The directory where the cloud dictionary will be stored is `$CLOUDDIR`, and the common string for the filenames is `$CLOUDSTRING`
- The threshold for cloud identification in the CO-luminous molecular gas column density is `$LOG_MIN_DENS`
- The minimum number of pixels in an identified cloud is `$MIN_PX`
- And finally, the number of CPUs for the parallel computation of cloud properties is `$NO_CPUS`
These definitions and default values can all be accessed via `python3 cloud_id.py --help`

### Checking outputs
Key properties of the output dictionaries can be checked `check_cloud-ID.ipynb`.

## `cloud_track.py`: Cloud tracking to produce cloud evolution network

### Dependencies
In addition to h5py, this script uses the [NetworkX](https://networkx.org/documentation/stable/index.html) network/graph analysis software.

### Inputs
- A set of hdf5 snapshots of an isolated galaxy simulation, containing a mesh of gas particles with centroid positions, velocities etc. Note that these _must_ form a continuous set with a fixed time-step.
- The set of cloud dictionaries generated via `cloud_id.py`

### Running
```python3 cloud_track.py $INPUT_DIR -s $SEARCH_STRING -p $CLOUDDIR/$CLOUDSTRING -b $FIRST_SNAP_NO -e $LAST_SNAP_NO -w $MAP_DIAMETER -r $MAP_RESOLUTION```

Where:
- Snapshots are found at `$INPUT_DIR`
- The search string for the snapshots is `$SEARCH_STRING`, e.g. `snap_*`
- The directory where the cloud dictionary will be stored is `$CLOUDDIR`, and the common string for the filenames is `$CLOUDSTRING`
- The first snapshot to analyze (time-ordered) is `$FIRST_SNAP_NO`
- The last snapshot to analyze (time-ordered) is `$LAST_SNAP_NO`
- The diameter of the maps in kiloparsec is `$MAP_DIAMETER`
- The resolution of the maps is `$MAP_RESOLUTION`
These definitions and default values can all be accessed via `python3 cloud_id.py --help`

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
