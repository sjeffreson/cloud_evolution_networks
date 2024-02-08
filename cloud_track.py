import sys, os, glob, re
regex = re.compile(r"\d+")
import numpy as np
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

def main(input_path, snap_list, cloudsnap_list, width_kpc, rsln_px):
    # initialize empty graph
    G = nx.DiGraph()

    # get info for first snapshot
    snap_data = ah.get_gas_info_from_snap(snap_list[0], width_kpc, rsln_px) # remember that time for this script is in Myr not Gyr
    cloud_dict_of_lists = np.load(cloudsnap_list[0], allow_pickle=True)
    print(type(cloud_dict_of_lists))
    cloud_dicts = [dict(zip(cloud_dict_of_lists, v)) for v in zip(*cloud_dict_of_lists.values())]
    cloud_lens = [len(c['particles']) for c in cloud_dicts]
    cloud_ids = [str(i) for i in range(len(cloud_dicts))]
    max_cloud_id = cloud_ids[-1]

    # add first level of nodes to graph
    # maybe we can do this all at once without looping through?
    for ID, cloud in zip(clouds_ids, clouds):            
        G.add_node(
            ID,
            time=snap_data['time'],
            mask=cloud['mask'],
            #pno=len(cloud['particles']),
            centroid=cloud['centroid'],
            vcentroid=cloud['vcentroid'],
            mass=cloud['mass'],
            div=cloud['div'],
            veldispz=cloud['veldispz'],
            veldisp3D=cloud['veldisp3D'],
            angmomz=cloud['angmomz'],
            angmomR=cloud['angmomR'],
            temp=cloud['temp'],
            SFR=cloud['SFR'],
        )
    print("first tree nodes added")

    # G_times = nx.get_node_attributes(G, 'time')
    # G_cents = nx.get_node_attributes(G, 'centroid')
    # G_colors = nx.get_node_attributes(G, 'color')
    # node_groups = sort_by_time(list(G_times.keys()), list(G_times.values()))
    # current_nodes = node_groups[len(node_groups)-1]

    # # project positions for the next snapshot
    # nextsnap = snap_list[1]
    # nextdata = h5py.snap(nextsnap, "r")
    # nextheader = nextdata['Header']

    # outpath = cloudsnapname.rsplit('/',1)[0]
    
    # # initialise loop
    # I = 1
    # while I < len(snap_list):

    #     # get the info from the next snapshot
    #     snapname_n = snap_list[I]
    #     data_n = h5py.snap(snapname_n, "r")
    #     gas_n = data_n['PartType0']
    #     header_n = data_n['Header']
    #     print("header/gas info extracted for snapshot " + snap_list[I].rsplit('snap-DESPOTIC_',1)[1].rsplit('.hdf5',1)[0])
        
    #     # get gas info
    #     x_n, y_n, vx_n, vy_n, masses_n, time_n = get_gas_info(header_n, gas_n, width_kpc)

    #     # get info from cloud snap
    #     cloudsnapname_n = cloudsnap_list[I]
    #     cloudsnap_n = np.load(cloudsnapname_n, allow_pickle=True)

    #     clouds_n, clouds_n_len = [], []
    #     for cloud_dict_of_lists in cloudsnap_n:
    #         cloud_list = [dict(zip(cloud_dict_of_lists, v)) for v in zip(*cloud_dict_of_lists.values())]
    #         clouds_n.append(cloud_list)
    #         clouds_n_len.append([len(c['particles']) for c in cloud_list])
    #     clouds_n = np.array(flatten_list(clouds_n))
    #     clouds_n_len = np.array(flatten_list(clouds_n_len))

    #     # remove clouds with < 1 particle
    #     clouds_n = clouds_n[clouds_n_len>1]

    #     # assign cloud ids
    #     clouds_n_ids = []
    #     for i in range(len(clouds_n)):
    #         clouds_n_ids += [str(int(max_clouds_id) + 1 + i)]
    #     max_clouds_id = clouds_n_ids[-1]
            
    #     # add new nodes to graph without any edges yet
    #     for name, cloud in zip(clouds_n_ids, clouds_n):
    #         G.add_node(
    #             name,
    #             time=time_n,
    #             pno=len(cloud['particles']),
    #             TrIDs=cloud['TracerIDs'],
    #             Tridcs=cloud['Traceridcs'],
    #             centroid=cloud['centroid'],
    #             H2centroid=cloud['H2centroid'],
    #             vcentroid=cloud['vcentroid'],
    #             H2vcentroid=cloud['H2vcentroid'],
    #             mass=cloud['mass'],
    #             H2mass=cloud['H2mass'],
    #             area=cloud['area'],
    #             div=cloud['div'],
    #             H2div=cloud['H2div'],
    #             veldisp=cloud['veldispz'],
    #             H2veldisp=cloud['H2veldispz'],
    #             veldisp3D=cloud['veldisp3D'],
    #             H2veldisp3D=cloud['H2veldisp3D'],
    #             Hz=cloud['Hz'],
    #             H2Hz=cloud['H2Hz'],
    #             HR=cloud['HR'],
    #             Htheta=cloud['Htheta'],
    #             temp=cloud['temp'],
    #             H2temp=cloud['H2temp'],
    #             SFR=cloud['SFR'],
    #             youngstarIDs=cloud['youngstarIDs'],
    #             starIDs=cloud['starIDs'],
    #             ID=name,
    #             mask=cloud['mask'],
    #             count=np.sum(cloud['mask']),
    #             color=np.random.rand(3,)
    #         )

    #     # add edges between two times
    #     # generate spatial bins
    #     x_edges = np.linspace(-width_kpc/2.*kpc_to_cm, width_kpc/2.*kpc_to_cm, rsln_px+1)
    #     y_edges = np.linspace(width_kpc/2.*kpc_to_cm, -width_kpc/2.*kpc_to_cm, rsln_px+1)

    #     # bin the projected clouds
    #     deltat = (time_n - time) * Myr_to_s
    #     x_bin_idcs, y_bin_idcs, x_vals, y_vals = {}, {}, {}, {}
    #     for name, cloud in zip(clouds_ids, clouds):
    #         part = cloud['particles']
    #         cent = cloud['centroid']
    #         x_proj = x[part] + deltat * vx[part]
    #         y_proj = y[part] + deltat * vy[part]
    #         x_bin_idcs[name] = np.digitize(x_proj, bins=x_edges, right=False)
    #         y_bin_idcs[name] = np.digitize(y_proj, bins=y_edges, right=False)

    #         cnd = ((x_bin_idcs[name]!=0) & (y_bin_idcs[name]!=0) &
    #             (x_bin_idcs[name]!=rsln_px+1) & (y_bin_idcs[name]!=rsln_px+1))
    #         x_bin_idcs[name] = x_bin_idcs[name][cnd] # remove values outside bin range
    #         y_bin_idcs[name] = y_bin_idcs[name][cnd]

    #         x_vals[name] = cent[0]
    #         y_vals[name] = cent[1]
    #     print("projected clouds binned")

    #     # bin the clouds for the next snapshot
    #     xn_bin_idcs, yn_bin_idcs, xn_vals, yn_vals = {}, {}, {}, {}
    #     for name, cloud in zip(clouds_n_ids, clouds_n):
    #         part = cloud['particles']
    #         cent = cloud['centroid']
    #         xn_bin_idcs[name] = np.digitize(x_n[part], bins=x_edges, right=False)
    #         yn_bin_idcs[name] = np.digitize(y_n[part], bins=y_edges, right=False)

    #         cnd = ((xn_bin_idcs[name]!=0) & (yn_bin_idcs[name]!=0) &
    #             (xn_bin_idcs[name]!=rsln_px+1) & (yn_bin_idcs[name]!=rsln_px+1))
    #         xn_bin_idcs[name] = xn_bin_idcs[name][cnd] # remove values outside bin range
    #         yn_bin_idcs[name] = yn_bin_idcs[name][cnd]

    #         xn_vals[name] = cent[0]
    #         yn_vals[name] = cent[1]
    #     print("next clouds binned")

    #     # initialise array to store the lists of clouds
    #     cloud_array_1D = np.empty((rsln_px,), dtype=object)
    #     for i in range(rsln_px):
    #         cloud_array_1D[i] = ''
    #     cloud_array = np.tile(cloud_array_1D, (rsln_px, 1))
    #     cloud_array_n = np.tile(cloud_array_1D, (rsln_px, 1))

    #     # store the list of projected clouds for each pixel
    #     for (k, x_idcs_), (k, y_idcs_) in zip(x_bin_idcs.items(), y_bin_idcs.items()):
    #         cloud_array[x_idcs_-1, y_idcs_-1] += k + ','
            
    #     # and now add the list of new clouds
    #     for (k, x_idcs_), (k, y_idcs_) in zip(xn_bin_idcs.items(), yn_bin_idcs.items()):
    #         cloud_array_n[x_idcs_-1, y_idcs_-1] += k + ','

    #     # binned masses of particles in next snap
    #     BoxSize = header_n.attrs['BoxSize'] * header.attrs['UnitLength_in_cm']
    #     mass_array_n, x_edges_bs, y_edges_bs, binnumber = binned_statistic_2d(
    #         (x_n+0.5*BoxSize)/kpc_to_cm, (y_n+0.5*BoxSize)/kpc_to_cm,
    #         masses_n,
    #         bins=((x_edges+0.5*BoxSize)/kpc_to_cm, (y_edges[::-1]+0.5*BoxSize)/kpc_to_cm),
    #         statistic='sum')
    #     mass_array_n = mass_array_n[:,::-1]

    #     # for each pixel that contains a non-empty string in
    #     # BOTH the previous and new arrays, link all real
    #     # current clouds to previous projected clouds (note
    #     # that each pixel can only belong to one current cloud,
    #     # although it can belong to multiple projected ones)
    #     cnd = (cloud_array!='') & (cloud_array_n!='')
    #     cloud_array_full = cloud_array[cnd]
    #     cloud_array_full_n = cloud_array_n[cnd]
    #     mass_array_full_n = mass_array_n[cnd]

    #     # pairs stores the proj/next cloud pairs that sit on a
    #     # shared pixel, pair_masses stores the CO-mass masked by that
    #     # one pixel
    #     pairs, pair_masses = [], []
    #     for elem, elem_n, mass in zip(cloud_array_full, cloud_array_full_n, mass_array_full_n):
    #         cloud_list = elem.split(",")
    #         cloud_list = list(filter(None, cloud_list))
    #         cloud_list_n = elem_n.split(",")
    #         cloud_list_n = list(filter(None, cloud_list_n))

    #         len_cloud_list = len(cloud_list)
    #         len_cloud_list_n = len(cloud_list_n)

    #         # add all projected clouds to the current cloud
    #         for i in range(len_cloud_list):
    #             pairs.append(np.array([cloud_list[i], cloud_list_n[0]]))
    #             pair_masses.append(mass)
    #     pairs = np.array(pairs)
    #     pair_masses = np.array(pair_masses)
    #     print("calculated non-unique pairs: there are "+str(len(pairs))+" pairs")

    #     # find the unique pairs, not worrying about permutations
    #     # because the pairs will already be ordered: proj, second
    #     unique_pairs, unique_counts = np.unique(pairs, axis=0, return_counts=True)

    #     # get the overlap mass/particle count shared by each pair
    #     overlap_masses, overlap_counts = [], []
    #     for unique_pair in unique_pairs:
    #         overlap_masses.append(np.sum(pair_masses[np.where((pairs[:,0]==unique_pair[0]) & (pairs[:,1]==unique_pair[1]))]))
    #     overlap_masses = np.array(overlap_masses)

    #     for pair, olmass, count in zip(unique_pairs, overlap_masses, unique_counts):
    #         pair = np.array(pair)
    #         if(count >= 1): # if the pixel overlap is 1 or more
    #             G.add_edge(pair[0], pair[-1], count=count, overlapmass=olmass)
    #     print("connected nodes found")

    #     # get next nodes
    #     G_times = nx.get_node_attributes(G, 'time')
    #     G_cents = nx.get_node_attributes(G, 'centroid')
    #     G_ids = nx.get_node_attributes(G, 'ID')
    #     node_groups = sort_by_time(list(G_times.keys()), list(G_times.values()))
    #     next_nodes = node_groups[len(node_groups)-1]

    #     # take colour of most massive parent
    #     for node in next_nodes:
    #         parents = [G.nodes[p] for p in G.predecessors(G_ids[node])]

    #         # sort, if the masses of the parents are different
    #         if(len(np.unique([parent['mass'] for parent in parents])) != len(parents)):
    #             parents_sorted = parents
    #         else:
    #             parents_sorted = [p for _,p in sorted(list(zip([parent['mass'] for parent in parents], parents)))]
    #         if len(parents) > 0:
    #             parent_color = parents_sorted[-1]['color']
    #             G.add_node(G_ids[node], color=parent_color)

    #     outpath = cloudsnapname_n.rsplit('/',1)[0]

    #     # save merger tree
    #     if((I==299) | (I==300)):
    #         #mergname = snapname_n.rsplit('/',1)[1].replace('.hdf5', '_mergertree-IH2thresh.gpickle')
    #         mergname = 'mergertree-IH2thresh-sIDs.gpickle'
    #         nx.write_gpickle(G, outpath + '/' + mergname)
    #         print("merger tree for snap " + snap_list[I].rsplit('snap-DESPOTIC_')[1].rsplit('.hdf5')[0] + ": " + outpath + '/' + mergname)

    #     # iterate
    #     current_nodes = next_nodes
    #     clouds = clouds_n
    #     clouds_ids = clouds_n_ids
    #     x = x_n
    #     y = y_n
    #     vx = vx_n
    #     vy = vy_n
    #     masses = masses_n
    #     time = time_n

    #     I += 1

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
            "Missing snapshots! {:s}".format(str(np.setdiff1d(np.arange(int(args.begsnapno), int(args.endsnapno)+1), np.array(snapnos, dtype=int))))
        )
        exit()
    if not np.array_equal(np.array(cloudnos_cut, dtype=int), np.arange(int(args.begsnapno), int(args.endsnapno)+1)):
        logger.critical(
            "Missing cloud snaps! {:s}".format(str(np.setdiff1d(np.arange(int(args.begsnapno), int(args.endsnapno)+1), np.array(cloudnos, dtype=int))))
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
