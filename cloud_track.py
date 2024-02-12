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

def main(input_path, snap_list, cloudfile_list, width, rsln_px):

    # initialize empty graph
    G = nx.DiGraph()
    min_cloud_id = 0
    prevtime = 0.

    for snap, cloudfile in zip(snap_list, cloudfile_list):
        snap_data = ah.get_gas_info_from_snap(snap, width, rsln_px)
        with open(cloudfile, 'rb') as f:
            cloud_dict_of_lists = pickle.load(f)
        cloud_attr_dicts = [dict(zip(cloud_dict_of_lists, v)) for v in zip(*cloud_dict_of_lists.values())]
        cloud_lens = [len(c['particles']) for c in cloud_attr_dicts]
        cloud_ids = [i for i in range(min_cloud_id, min_cloud_id+len(cloud_attr_dicts))]

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

        x_bin_edges = np.linspace(-width*ah.kpc_to_cm/2., width*ah.kpc_to_cm/2., rsln_px+1)
        y_bin_edges = np.linspace(width*ah.kpc_to_cm/2., -width*ah.kpc_to_cm/2., rsln_px+1)

        pixel_ids_proj = {} # pixels that contain part of a projected cloud
        for id, attrs, in zip(cloud_ids, cloud_attr_dicts):
            xproj = snap_data['x_coords'][attrs['particles']] - deltat * snap_data['velxs'][attrs['particles']]
            yproj = snap_data['y_coords'][attrs['particles']] - deltat * snap_data['velys'][attrs['particles']]

            cnd = (np.fabs(xproj) < width*ah.kpc_to_cm/2.) & (np.fabs(yproj) < width*ah.kpc_to_cm/2.)
            xproj_bin_idcs = np.digitize(xproj[cnd], bins=x_bin_edges)
            yproj_bin_idcs = np.digitize(yproj[cnd], bins=y_bin_edges)
            pxsproj = (xproj_bin_idcs-1) * rsln_px + (yproj_bin_idcs-1) # recall 90-degree rotation of masks
            for pxproj in pxsproj:
                if pxproj in pixel_ids_proj:
                    pixel_ids_proj[pxproj].append(id)
                else:
                    pixel_ids_proj[pxproj] = [id]

        # now see which pixels are occupied by clouds in the previous snapshot, for comparison
        # NOTE: only one previous cloud can occupy a pixel
        G_inttimes = nx.get_node_attributes(G, 'inttime')
        G_masks = nx.get_node_attributes(G, 'mask')
        G_masks_prev = {id: mask for id, mask in G_masks.items() if G_inttimes[id] == previnttime}

        pixel_ids_prev = {} # pixels that contain part of a real cloud during the previous time-step
        for id, mask in G_masks_prev.items():
            xprev_bin_idcs, yprev_bin_idcs = mask
            pxsprev = yprev_bin_idcs * rsln_px + xprev_bin_idcs
            for pxprev in pxsprev:
                pixel_ids_prev[pxprev] = id

        # compare the pixel dictionaries and link cloud ids that share a pixel
        new_edges = []
        for px, proj_ids in pixel_ids_proj.items():
            if px in pixel_ids_prev:
                prev_id = pixel_ids_prev[px]
                for proj_id in proj_ids:
                    new_edges.append((prev_id, proj_id))
        new_edges = set(new_edges)
        G.add_edges_from(new_edges)

        logger.info(
            "Added {:d} edges between clouds in snapshot {:s} and snapshot {:s}".format(
                len(new_edges), str(int(previnttime)), str(int(np.round(snap_data['time']*ah.Gyr_to_Myr)))
            )
        )

        # update the previous time
        min_cloud_id = cloud_ids[-1]+1
        prevtime = snap_data['time']
        previnttime = int(np.round(snap_data['time']*ah.Gyr_to_Myr))

    # save merger tree
    with open(os.path.join(input_path, "merger_tree_"+str(args.begsnapno)+"-"+str(args.endsnapno)+".pkl"), "wb") as f:
        pickle.dump(G, f, pickle.HIGHEST_PROTOCOL)
    logger.info(
        "Saved merger tree to {:s}".format(os.path.join(input_path, "merger_tree_"+str(args.begsnapno)+"-"+str(args.endsnapno)+".pkl"))
    )

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
