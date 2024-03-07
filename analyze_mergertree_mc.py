import h5py
import numpy as np
import glob, os, re, sys
import pickle

import astro_helper as ah
from MergertreesProps import MergertreeProps

num_mc_iter = int(sys.argv[1])
num_mc_workers = int(sys.argv[2])

# MergertreeProps instance
mergertree_name = 'merger_tree_601-799.pkl'
props = MergertreeProps(mergertree_name, num_mc_iter=num_mc_iter, num_mc_workers=num_mc_workers)

props.cut_wcs(time_cut=True, radius_cut=True, time_rsln_cut=True)

props.save_cloud_evol_mc()