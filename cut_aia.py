#!/home/potato/miniforge3/envs/sunpy/bin/python

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import os
import sunpy.map
from astropy.coordinates import SkyCoord

directories = {
    "94": "cut_94",
    "131": "cut_131",
    "171": "cut_171",
    "193": "cut_193",
    "211": "cut_211",
    "304": "cut_304",
    "335": "cut_335"
}

for input_dir, output_dir in directories.items():
    
    os.makedirs(output_dir, exist_ok=True)
    
    fits_files = [f for f in os.listdir(input_dir) if f.endswith(".fits")]

    for fits_file in fits_files:
        # read the fits
        input_file_path = os.path.join(input_dir, fits_file) 
        swap_map = sunpy.map.Map(input_file_path)
    
        # cut the fits
        top_right = SkyCoord(1200 * u.arcsec, -100 * u.arcsec, frame=swap_map.coordinate_frame)
        bottom_left = SkyCoord(700 * u.arcsec, -600 * u.arcsec, frame=swap_map.coordinate_frame)
        swap_submap = swap_map.submap(bottom_left, top_right=top_right)
    
        # save new fits
        output_file_path = os.path.join(output_dir, fits_file.replace(".fits", "_cut.fits"))
        swap_submap.save(output_file_path)
        print(f"Processed and saved: {output_file_path}")

