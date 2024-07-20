#!/home/potato/miniforge3/envs/sunpy/bin/python

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import os
import sunpy.map
from astropy.coordinates import SkyCoord

input_dir = "cut_171/"
output_dir = "fig_171/"
fits_files = [f for f in os.listdir(input_dir) if f.endswith(".fits")]

for fits_file in fits_files:
    input_file_path = os.path.join(input_dir, fits_file) 
    a = sunpy.map.Map(input_file_path)
    a.plot(autoalign=True)
    output_file_path = os.path.join(output_dir, fits_file.replace(".fits", ".png"))
    plt.savefig(output_file_path)
    print(f"Processed and saved: {output_file_path}")
    plt.clf()

