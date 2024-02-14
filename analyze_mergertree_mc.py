import h5py
import numpy as np
import glob, os, re, sys
import pickle

import astro_helper as ah
from MergertreesProps import MergertreeProps

# MergertreeProps instance
mergertree_name = 'merger_tree_500-799.pkl'
props = MergertreeProps(mergertree_name, num_mc_iter=2, num_mc_workers=2)

print(props.get_timestep())
print(props.get_width()/ah.kpc_to_cm)

print(props.num_self_loops())

props.cut_wcs(time_cut=True, radius_cut=True, time_rsln_cut=True)

print(props.num_self_loops())

props.save_cloud_evol_mc()