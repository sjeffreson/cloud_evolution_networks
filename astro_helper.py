# Helper functions and unit conversions that come in handy for astronomy projects
import numpy as np
import h5py

# standard units
Msol_to_g = 1.99E33
kpc_to_cm = 3.086E21
pc_to_cm = 3.086E18
LCO_to_H2mass = 2.7E28
gamma = 5./3.
kB_cgs = 1.38E-16
mp_cgs = 1.67E-24
mu = 1.4

# Arepo units
ArepoMass_to_g = 1.989E43
ArepoLength_to_cm = 3.085678E21

def get_image_data(filename):
	with open(filename, "rb") as f:
		xpix = np.fromfile(f, dtype=np.int32, count=1)[0]
		ypix = np.fromfile(f, dtype=np.int32, count=1)[0]
		img = np.fromfile(f, dtype=np.float32, count=xpix*ypix)
	img = np.reshape(img, (xpix, ypix))
	img = np.rot90(img)
	return img

def flatten_list(lst):
    if type(lst[0])==list:
        return [item for sublist in lst for item in sublist]
    else:
        return lst

def get_gas_info_from_snap(snapname, width, rsln_px):
	snap = h5py.File(snapname, "r")
	header = snap["Header"]
	gas = snap["PartType0"]

	x_coords = (gas['Coordinates'][:,0] - 0.5 * header.attrs['BoxSize']) * gas['Coordinates'].attrs['to_cgs']
	y_coords = (gas['Coordinates'][:,1] - 0.5 * header.attrs['BoxSize']) * gas['Coordinates'].attrs['to_cgs']
	R_coords = np.sqrt(x_coords*x_coords + y_coords*y_coords)
	z_coords = (gas['Coordinates'][:,2] - 0.5 * header.attrs['BoxSize']) * gas['Coordinates'].attrs['to_cgs']
	cnd = (np.fabs(x_coords) < width*kpc_to_cm/2.) & (np.fabs(y_coords) < width*kpc_to_cm/2.)

	snap_data = {}
	snap_data["time"] = header.attrs["Time"]
	snap_data["x_coords"] = x_coords[cnd]
	snap_data["y_coords"] = y_coords[cnd]
	snap_data["z_coords"] = z_coords[cnd]
	snap_data["R_coords"] = R_coords[cnd]
	snap_data["velxs"] = gas['Velocities'][:,0][cnd] * gas['Velocities'].attrs['to_cgs']
	snap_data["velys"] = gas['Velocities'][:,1][cnd] * gas['Velocities'].attrs['to_cgs']
	snap_data["velzs"] = gas['Velocities'][:,2][cnd] * gas['Velocities'].attrs['to_cgs']
	snap_data["masses"] = gas['Masses'][:][cnd] * gas['Masses'].attrs['to_cgs']
	snap_data["IH2s"] = gas['ChemicalAbundances'][:,0][cnd]
	snap_data["voldenses"] = gas['Density'][:][cnd] * gas['Density'].attrs['to_cgs']
	snap_data["Us"] = gas['InternalEnergy'][:][cnd] * gas['InternalEnergy'].attrs['to_cgs']
	snap_data["temps"] = (gamma - 1.) * snap_data["Us"]/kB_cgs * mu * mp_cgs
	snap_data["SFRs"] = gas['StarFormationRate'][:][cnd] * gas['StarFormationRate'].attrs['to_cgs']

	# digitize onto grid of size rsln_px x rsln_px
	x_bin_edges = np.linspace(-width*kpc_to_cm/2., width*kpc_to_cm/2., rsln_px+1)
	y_bin_edges = np.linspace(width*kpc_to_cm/2., -width*kpc_to_cm/2., rsln_px+1)
	snap_data["x_bin_idx"] = np.digitize(snap_data["x_coords"], x_bin_edges)-1
	snap_data["y_bin_idx"] = np.digitize(snap_data["y_coords"], y_bin_edges)-1

	return snap_data