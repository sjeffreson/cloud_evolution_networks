# Helper functions and unit conversions that come in handy for astronomy projects
import numpy as np

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

