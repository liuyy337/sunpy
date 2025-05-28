#!/home/potato/miniforge3/envs/sunpy/bin/python

import astropy.units as u
import gc
import matplotlib.pyplot as plt
import numpy as np
import os
import sunpy.map
import time
from aiapy.calibrate import register, update_pointing
from aiapy.calibrate.util import get_pointing_table
from astropy.coordinates import SkyCoord
from astropy.io import fits
# from datetime import datetime as dtdt
from multiprocessing import Pool, cpu_count
from sunpy.coordinates import propagate_with_solar_surface

directory = {"/media/potato/solar/data/20240511/aia/0094/": "cutted/0094/",
             "/media/potato/solar/data/20240511/aia/0131/": "cutted/0131/",
             "/media/potato/solar/data/20240511/aia/0171/": "cutted/0171/",
             "/media/potato/solar/data/20240511/aia/0193/": "cutted/0193/",
             "/media/potato/solar/data/20240511/aia/0211/": "cutted/0211/",
             "/media/potato/solar/data/20240511/aia/0304/": "cutted/0304/",
             "/media/potato/solar/data/20240511/aia/0335/": "cutted/0335/",
             "/media/potato/solar/data/20240511/aia/1600/": "cutted/1600/",
             "/media/potato/solar/data/20240511/aia/1700/": "cutted/1700/"}
x0, y0 = 715, -370
dx, dy = 150, 150

def create_header(original_header, reprojected_header):
    new_header = fits.Header()

    new_header['SIMPLE'] = original_header['SIMPLE']
    new_header['BITPIX'] = original_header['BITPIX']
    new_header['NAXIS'] = 2
    new_header['EXPTIME'] = original_header['EXPTIME']
    new_header['T_REC'] = original_header['T_REC']
    new_header['T_OBS'] = original_header['T_OBS']
    new_header['DATE-OBS'] = original_header['DATE-OBS']
    
    # WCS
    new_header['CTYPE1'] = original_header['CTYPE1']
    new_header['CTYPE2'] = original_header['CTYPE2']
    new_header['CRVAL1'] = reprojected_header['CRVAL1']
    new_header['CRVAL2'] = reprojected_header['CRVAL2']
    new_header['CRPIX1'] = reprojected_header['CRPIX1']
    new_header['CRPIX2'] = reprojected_header['CRPIX2']
    new_header['CDELT1'] = original_header['CDELT1']
    new_header['CDELT2'] = original_header['CDELT2']
    new_header['CUNIT1'] = original_header['CUNIT1']
    new_header['CUNIT2'] = original_header['CUNIT2']
    
    new_header['TELESCOP'] = original_header['TELESCOP']
    new_header['INSTRUME'] = original_header['INSTRUME']
    new_header['WAVELNTH'] = original_header['WAVELNTH']
    new_header['WAVEUNIT'] = original_header['WAVEUNIT']
    
    new_header['DSUN_OBS'] = original_header['DSUN_OBS']
    new_header['RSUN_OBS'] = original_header['RSUN_OBS']
    new_header['HGLT_OBS'] = original_header['HGLT_OBS']
    new_header['HGLN_OBS'] = original_header['HGLN_OBS']
    
    return new_header

def cut_aia(args):
    # paramenters
    input_dir, output_dir, file, target_wcs, pointing_table = args
    os.makedirs(output_dir, exist_ok = True)
    input_path = os.path.join(input_dir, file)
    output_path = os.path.join(output_dir, file.replace("_lev1.fits", ".fits"))

    # aia_prep: registering and aligning level 1 data
    aia_map = sunpy.map.Map(input_path)
    aia_map_updated_pointing = update_pointing(aia_map, pointing_table=pointing_table)
    aia_registered = register(aia_map_updated_pointing)

    with propagate_with_solar_surface():
        reprojected_map = aia_registered.reproject_to(target_wcs)
    new_header = create_header(aia_registered.meta, reprojected_map.meta)
    data_int16 = reprojected_map.data.astype(np.int16)

    final_map = sunpy.map.Map(data_int16, new_header)
    final_map.save(output_path, overwrite=True)

def main():
    # prepare reference map
    ref_path = "/media/potato/solar/data/20240511/aia/0335/aia.lev1_euv_12s.2024-05-11T105926Z.335.image_lev1.fits"
    ref_map = sunpy.map.Map(ref_path)
    pointing_table = get_pointing_table("JSOC", time_range=(ref_map.date - 12 * u.h, ref_map.date + 12 * u.h))
    ref_map_updated_pointing = update_pointing(ref_map, pointing_table=pointing_table)
    ref_registered = register(ref_map_updated_pointing)
    
    corner = SkyCoord(x0 * u.arcsec, y0 * u.arcsec, frame=ref_registered.coordinate_frame)
    cutout_map = ref_registered.submap(corner, width = dx * u.arcsec, height = dy * u.arcsec)
    print(cutout_map.data.shape)
    
    # multiprocessing
    for input_dir, output_dir in directory.items():
        files = sorted([f for f in os.listdir(input_dir) if f.endswith('.fits')])
        args = [(input_dir, output_dir, file, cutout_map.wcs, pointing_table) for file in files]
        with Pool(cpu_count()) as pool:
            pool.map(cut_aia, args)
        gc.collect()

if __name__ == "__main__":
    start_time = time.time()
    main()
    print(f"Total processing time: {time.time()-start_time:.2f} seconds")

