import sys, os, glob, re
regex = re.compile(r"\d+")
import numpy as np
import argparse, logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

import multiprocessing
from functools import partial
from contextlib import contextmanager

import pickle

import h5py
import astrodendro
from astrodendro import Dendrogram
import astro_helper as ah

def get_galaxy_dendrogram(LCOimgname, imgname, min_value = -1.5, min_npix=1):
    LCOimg = ah.get_image_data(LCOimgname)
    img = ah.get_image_data(imgname)
    LCOmass = (
        img * LCOimg * ah.LCO_to_H2mass * ah.ArepoMass_to_g
        / ah.Msol_to_g / (ah.ArepoLength_to_cm / ah.pc_to_cm)**2
    )

    d = Dendrogram.compute(np.log10(LCOmass), min_value=min_value, min_npix=min_npix)

    return d

@contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

def clouds_for_one_img(mask_idcs_group, snap_data=None, masks_where=None, rsln_px=None):
    logger.info(
        "Clouds {:d} to {:d}".format(mask_idcs_group[0], mask_idcs_group[-1])
    )

    clouds_dict = {
        'particles': [],
        'mask': [],
        'centroid': [],
        'vcentroid': [],
        'mass': [],
        'veldispz': [],
        'veldisp3D': [],
        'divergence': [],
        'angmomz': [],
        'angmomR': [],
        'temp': [],
        'starformrate': []
    }
    for idx in mask_idcs_group:
        mask_where = masks_where[idx]
        cloud_mask = np.zeros((rsln_px, rsln_px), dtype=bool)
        cloud_mask[mask_where] = True

        # remember that the images are rotated by 90 degrees relative to the gas particle co-ordinates
        cloud = cloud_mask[snap_data["y_bin_idx"]-1, snap_data["x_bin_idx"]-1].nonzero()[0] # gas particles

        # cloud masses
        masses_cl = snap_data["masses"][cloud] * snap_data["IH2s"][cloud] * 2. / ah.mu

        # ensure against clouds with no mass, can happen if the min px number per cloud is small
        if(np.sum(masses_cl) > 0.):
            # center of mass and velocity
            xC_cl = np.average(snap_data["x_coords"][cloud], weights=masses_cl)
            yC_cl = np.average(snap_data["y_coords"][cloud], weights=masses_cl)
            zC_cl = np.average(snap_data["z_coords"][cloud], weights=masses_cl)

            vxC_cl = np.average(snap_data["velxs"][cloud], weights=masses_cl)
            vyC_cl = np.average(snap_data["velys"][cloud], weights=masses_cl)
            vzC_cl = np.average(snap_data["velzs"][cloud], weights=masses_cl)

            # co-ordinates and velocities relative to cloud center of mass
            xrelCs_cl = snap_data["x_coords"][cloud] - xC_cl
            yrelCs_cl = snap_data["y_coords"][cloud] - yC_cl
            zrelCs_cl = snap_data["z_coords"][cloud] - zC_cl

            vxrelCs_cl = snap_data["velxs"][cloud] - vxC_cl
            vyrelCs_cl = snap_data["velys"][cloud] - vyC_cl
            vzrelCs_cl = snap_data["velzs"][cloud] - vzC_cl

            # velocity dispersion in galactic z-direction and in 3D
            veldispz_cl = np.sqrt(np.average(vzrelCs_cl**2, weights=masses_cl))
            veldisp3D_cl = np.sqrt(
                np.average(vxrelCs_cl**2, weights=masses_cl)+
                np.average(vyrelCs_cl**2, weights=masses_cl)+
                np.average(vzrelCs_cl**2, weights=masses_cl)
            )

            # divergence of the velocity field inside clouds
            div_cl = np.average(
                ((vxrelCs_cl * xrelCs_cl + vyrelCs_cl * yrelCs_cl + vzrelCs_cl * zrelCs_cl) /
                    np.sqrt(xrelCs_cl**2 + yrelCs_cl**2 + zrelCs_cl**2)
                ), weights=masses_cl
            )

            # angular momentum
            angmomx_cl = np.average(yrelCs_cl * vzrelCs_cl - zrelCs_cl * vyrelCs_cl, weights=masses_cl)
            angmomy_cl = np.average(zrelCs_cl * vxrelCs_cl - xrelCs_cl * vzrelCs_cl, weights=masses_cl)
            angmomz_cl = np.average(xrelCs_cl * vyrelCs_cl - yrelCs_cl * vxrelCs_cl, weights=masses_cl)
            angmomR_cl = (-angmomx_cl * xrelCs_cl - angmomy_cl * yrelCs_cl) / np.sqrt(xrelCs_cl**2 + yrelCs_cl**2)

            # temperature
            temp_cl = np.average(snap_data["temps"][cloud], weights=masses_cl)

            # star formation rate
            SFR_cl = np.sum(snap_data["SFRs"][cloud])

            # save to dictionary
            clouds_dict["particles"].append(np.array(cloud))
            clouds_dict["mask"].append(cloud_mask)
            clouds_dict["centroid"].append(np.array([xC_cl, yC_cl, zC_cl]))
            clouds_dict["vcentroid"].append(np.array([vxC_cl, vyC_cl, vzC_cl]))
            clouds_dict["mass"].append(np.sum(masses_cl))
            clouds_dict["veldispz"].append(veldispz_cl)
            clouds_dict["veldisp3D"].append(veldisp3D_cl)
            clouds_dict["divergence"].append(div_cl)
            clouds_dict["angmomz"].append(angmomz_cl)
            clouds_dict["angmomR"].append(angmomR_cl)
            clouds_dict["temp"].append(temp_cl)
            clouds_dict["starformrate"].append(SFR_cl)

    return clouds_dict

def main(imgnames, width, rsln_px, prefix, logminvalue, minpix, Ncores):
    # Loop over time in the galaxy simulation (one time == one snapshot)
    for imgname in imgnames:
        LCOimgname = imgname.replace("density_proj", "LCO_proj")
        d = get_galaxy_dendrogram(LCOimgname, imgname, min_value = float(args.logminvalue), min_npix = int(args.minpix))
        
        # the trunk of the dendrogram is the list of clouds
        trunk = d.trunk
        logger.info(
            "Number of clouds in dendrogram is {:d}".format(len(trunk))
        )

        # one branch of the dendrogram is one cloud
        masks_where = []
        for branch in trunk:
            cloud_mask = branch.get_mask()
            masks_where.append(np.where(cloud_mask==1)) # pixels
        mask_idcs = list(range(len(masks_where)))
        print("done with mask computation")

        # the snapshot contains all the gas particle information
        snap_data = ah.get_gas_info_from_snap(imgname.replace("density_proj", "snap-DESPOTIC")+".hdf5", width, rsln_px)

        # divide clouds (branches) into groups for the number of cores
        N_pergroup = int(np.ceil(len(masks_where)/Ncores))
        mask_idcs_grouped = [
            mask_idcs[i*N_pergroup:i*N_pergroup+N_pergroup]
            for i in range(Ncores) if i*N_pergroup<len(mask_idcs)
        ]

        # compute cloud properties in parallel
        clouds_for_one_img_data = partial(clouds_for_one_img, snap_data=snap_data, masks_where=masks_where, rsln_px=rsln_px)

        with poolcontext(processes=Ncores) as pool:
            clouds_dict = pool.map(clouds_for_one_img_data, mask_idcs_grouped)

        # print the whole clouds_dict with the logger
        logger.info(len(clouds_dict))

        # combine the results from the different cores
        clouds_dict = {
            key: ah.flatten_list([clouds_dict[i][key] for i in range(len(clouds_dict))]) for key in clouds_dict[0].keys()
        }

        with open(args.prefix + "_" + regex.findall(imgname)[-1] + ".pickle", 'wb') as f:
            pickle.dump(clouds_dict, f)
    
        logger.info(
            "Cloud properties for snapshot saved to file {:s}".format(args.prefix + "_" + regex.findall(imgname)[-1] + ".pickle")
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Cloud tracking",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input",
        help="Input folder containing galaxy snapshot images to be analyzed"
    )
    parser.add_argument(
        "-s", "--search",
        default="*",
        help="Search pattern to find galaxy snapshot images in input dir"
    )
    parser.add_argument(
        "-b", "--begsnapno",
        default="300",
        help="Will not analyze galaxy snapshots before this number"
    )
    parser.add_argument(
        "-e", "--endsnapno",
        default="600",
        help="Will not analyze galaxy snapshots after this number"
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
    parser.add_argument(
        "-p", "--prefix",
        default="clouds",
        help="Prefix for cloud filenames"
    )
    parser.add_argument(
        "-lmv", "--logminvalue",
        default="-1.5",
        help="The base 10 logarithm of the minimum H2 surface density to be considered a cloud, in Msol/pc^2"
    )
    parser.add_argument(
        "-mp", "--minpix",
        default="9",
        help="The minimum number of pixels to be considered a cloud"
    )
    parser.add_argument(
        "-N", "--Ncores",
        default="48",
        help="Number of cores for the parallel computation of cloud properties"
    )
    args = parser.parse_args()

    if not os.path.exists(args.input):
        logger.critical(
            "Input folder cannot be found: {:s}".format(args.input)
        )
        exit()

    # list of galaxy snapshot images
    imgnames = [imgname for imgname in sorted(glob.glob(os.path.join(args.input, args.search)))]

    # cut out the image names that are not within the specified range, or
    # for which cloud files already exist
    imgnos = [regex.findall(imgname)[-1] for imgname in imgnames]
    exists = [os.path.exists(args.prefix + "_" + imgno + ".npy") for imgno in imgnos]
    imgnames_cut = [
        imgname for imgname, imgno, exist in zip(imgnames, imgnos, exists)
        if not exist and int(imgno) >= int(args.begsnapno) and int(imgno) <= int(args.endsnapno)
    ]
    if len(imgnames_cut) == 0:
        logger.info("No images to analyze")
        exit()

    logger.info(
        "Analyzing {:d} snapshots in the interval {:s} to {:s}".format(
            len(imgnames_cut), imgnames_cut[0], imgnames_cut[-1]
        )
    )

    main(
        imgnames_cut,
        float(args.width),
        int(args.rsln_px),
        args.prefix,
        float(args.logminvalue),
        int(args.minpix),
        int(args.Ncores),
    )