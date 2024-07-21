#!/home/potato/miniforge3/envs/sunpy/bin/python

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import os
import sunpy.map
from astropy.coordinates import SkyCoord

directories = {
    "cut_94": "fig_94",
#     "cut_131": "fig_131",
    "cut_171": "fig_171",
    "cut_193": "fig_193",
    "cut_211": "fig_211",
    "cut_304": "fig_304",
    "cut_335": "fig_335"
}

for input_dir, output_dir in directories.items():
	
	os.makedirs(output_dir, exist_ok=True)

	fits_files = [f for f in os.listdir(input_dir) if f.endswith(".fits")]

	for fits_file in fits_files:
		#read the fits
		input_file_path = os.path.join(input_dir, fits_file) 
		a = sunpy.map.Map(input_file_path)
		
		#plot the fits
		a.plot(autoalign=True)
		output_file_path = os.path.join(output_dir, fits_file.replace(".fits", ".png"))
		plt.savefig(output_file_path)
		print(f"Processed and saved: {output_file_path}")
		plt.clf()

