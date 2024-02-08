import sys, os, glob, re
regex = re.compile(r"\d+")
import numpy as np
import pickle
import argparse, logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

import h5py
from scipy.stats import binned_statistic_2d
import networkx as nx
import astro_helper as ah

def sort_by_time(nodes, times):
    sidx = np.argsort(times)
    split_idx = np.flatnonzero(np.diff(np.take(times,sidx))>0)+1
    out = np.split(np.take(nodes,sidx,axis=0), split_idx)
    out = [list(elem) for elem in out]
    return out

def main(input_path, snap_list, cloudfile_list, width_kpc, rsln_px):

    # initialize empty graph
    G = nx.DiGraph()
    min_cloud_id = 0
    prevtime = 0.

    for snap, cloudfile in zip(snap_list, cloudfile_list):
        snap_data = ah.get_gas_info_from_snap(snap, width_kpc, rsln_px)
        with open(cloudfile, 'rb') as f:
            cloud_dict_of_lists = pickle.load(f)
        cloud_attr_dicts = [dict(zip(cloud_dict_of_lists, v)) for v in zip(*cloud_dict_of_lists.values())]
        cloud_lens = [len(c['particles']) for c in cloud_attr_dicts]
        cloud_ids = [str(i) for i in range(min_cloud_id, len(cloud_attr_dicts))]

        # add nodes to graph
        G.add_nodes_from(cloud_ids)
        for id, attrs, in zip(cloud_ids, cloud_attr_dicts):
            nx.set_node_attributes(G, {id: {'time': snap_data['time']}})
            nx.set_node_attributes(G, {id: {'inttime': int(np.round(snap_data['time']*ah.Gyr_to_Myr))}})
            nx.set_node_attributes(G, {id: attrs})

        if(snap_list.index(snap) == 0):
            min_cloud_id = cloud_ids[-1]+1
            prevtime = snap_data['time']
            previnttime = int(np.round(snap_data['time']*ah.Gyr_to_Myr))
            continue

        # predict positions of these clouds in the last snapshot, using the velocities
        # of their particles, and bin for comparison to current masks
        deltat = (snap_data['time'] - prevtime) * ah.Gyr_to_s

        x_bin_edges = np.linspace(-width_kpc/2.*kpc_to_cm, width_kpc/2.*kpc_to_cm, rsln_px+1)
        y_bin_edges = np.linspace(width_kpc/2.*kpc_to_cm, -width_kpc/2.*kpc_to_cm, rsln_px+1)

        pixel_ids_proj = {}
        for id, attrs, in zip(cloud_ids, cloud_attr_dicts):
            xproj = snap_data['x_coords'][attrs['particles']] - deltat * snap_data['velxs'][attrs['particles']]
            yproj = snap_data['y_coords'][attrs['particles']] - deltat * snap_data['velys'][attrs['particles']]

            cnd = (np.fabs(xproj) < width*kpc_to_cm/2.) & (np.fabs(yproj) < width*kpc_to_cm/2.)
            xproj_bin_idcs = np.digitize(xproj[cnd], bins=x_bin_edges)-1
            yproj_bin_idcs = np.digitize(yproj[cnd], bins=y_bin_edges)-1
            pxsproj = yproj_bin_idcs * rsln_px + xproj_bin_idcs
            for pxproj in pxsproj:
                if pxproj in pixel_ids_proj:
                    pixel_ids_proj[pxproj].append(id) # pixels which the projected clouds overlap
                else:
                    pixel_ids_proj[pxproj] = [id]

        # similar pixel dictionary from the real positions of the previous clouds
        G_inttimes = nx.get_node_attributes(G, 'inttime')
        G_masks = nx.get_node_attributes(G, 'mask')
        G_masks_prev = {id: mask for id, mask in G_masks.items() if G_inttimes[id] == previnttime}

        pixel_ids_prev = {}
        for id, mask in G_masks_prev.items():
            xprev_bin_idcs, yprev_bin_idcs = mask
            pxsprev = yprev_bin_idcs * rsln_px + xprev_bin_idcs
            for pxprev in pxsprev:
                pixel_ids_prev[pxprev] = id # there's only one previous cloud on each pixel (true position)

        # compare the pixel dictionaries and link cloud ids that share a pixel
        for px, proj_ids in pixel_ids_proj.items():
            prev_id = pixel_ids_prev[px]
            for proj_id in proj_ids:
                G.add_edge(prev_id, proj_id)

        # update the previous time
        min_cloud_id = cloud_ids[-1]+1
        prevtime = snap_data['time']
        previnttime = int(np.round(snap_data['time']*ah.Gyr_to_Myr))

    # save merger tree
    nx.write_gpickle(G, os.path.join(input_path, "merger_tree.pkl"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Cloud tracking",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input",
        help="Input folder containing galaxy snapshots to be analyzed"
    )
    parser.add_argument(
        "-s", "--search",
        default="*",
        help="Search pattern to find galaxy snapshots in input dir"
    )
    parser.add_argument(
        "-p", "--prefix",
        default="clouds",
        help="Prefix for cloud filenames in input dir"
    )
    parser.add_argument(
        "-b", "--begsnapno",
        default="300",
        help="Snapshot at which to begin cloud tracking"
    )
    parser.add_argument(
        "-e", "--endsnapno",
        default="600",
        help="Snapshot at which to end cloud tracking"
    )
    parser.add_argument(
        "-w", "--width",
        default="30",
        help="Diameter of galaxy in kpc"
    )
    parser.add_argument(
        "-r", "--rsln_px",
        default="5000",
        help="Resolution of the galaxy snapshot in pixels (one side)"
    )
    args = parser.parse_args()

    if not os.path.exists(args.input):
        logger.critical(
            "Input folder cannot be found: {:s}".format(args.input)
        )
        exit()

    # list of galaxy snapshot hdf5 snaps
    snapnames = [snapname for snapname in sorted(glob.glob(os.path.join(args.input, args.search)))]
    cloudnames = [cloudname for cloudname in sorted(glob.glob(os.path.join(args.input, args.prefix+"*")))]

    # cut out the image names that are not within the specified range
    snapnames_cut = [
        snapname for snapname in snapnames 
        if(int(regex.findall(snapname.replace(".hdf5", ""))[-1]) >= int(args.begsnapno)) &
        (int(regex.findall(snapname.replace(".hdf5", ""))[-1]) <= int(args.endsnapno))
    ]
    cloudnames_cut = [
        cloudname for cloudname in cloudnames 
        if(int(regex.findall(cloudname)[-1]) >= int(args.begsnapno)) &
        (int(regex.findall(cloudname)[-1]) <= int(args.endsnapno))
    ]

    # check that all snapshots and all cloud snaps exist within this range,
    # at intervals of 1 Myr
    snapnos_cut = [regex.findall(snapname.replace(".hdf5", ""))[-1] for snapname in snapnames_cut]
    cloudnos_cut = [regex.findall(cloudname)[-1] for cloudname in cloudnames_cut]
    if not np.array_equal(np.array(snapnos_cut, dtype=int), np.arange(int(args.begsnapno), int(args.endsnapno)+1)):
        logger.critical(
            "Missing snapshots! {:s}".format(str(np.setdiff1d(np.arange(int(args.begsnapno), int(args.endsnapno)+1), np.array(snapnos_cut, dtype=int))))
        )
        exit()
    if not np.array_equal(np.array(cloudnos_cut, dtype=int), np.arange(int(args.begsnapno), int(args.endsnapno)+1)):
        logger.critical(
            "Missing cloud snaps! {:s}".format(str(np.setdiff1d(np.arange(int(args.begsnapno), int(args.endsnapno)+1), np.array(cloudnos_cut, dtype=int))))
        )
        exit()

    logger.info(
        "Performing cloud tracking for snapshots/clouds in the interval {:s} to {:s}".format(
            snapnos_cut[0], snapnos_cut[-1]
        )
    )

    main(
        args.input,
        snapnames_cut,
        cloudnames_cut,
        float(args.width),
        int(args.rsln_px),
    )
