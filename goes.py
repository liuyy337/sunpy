#!/home/potato/miniforge3/envs/sunpy/bin/python

import matplotlib.pyplot as plt
import numpy as np
import os

from astropy.visualization import time_support
from sunpy import timeseries as ts
from sunpy.net import Fido
from sunpy.net import attrs as a

tstart = "2024-05-05 00:00"
tend = "2024-05-15 00:00"

# download goes data
result_goes16 = Fido.search(a.Time(tstart, tend), a.Instrument("XRS"), a.goes.SatelliteNumber(16), a.Resolution("flx1s"))
Fido.fetch(result_goes16)

# plot goes data
input_dir = '/home/potato/sunpy/data'
output_dir = '/home/potato/sunpy/data'
files = os.listdir(input_dir)
nc_files = [file for file in files if file.endswith('.nc')]

for file in nc_files:
	#read nc files
	input_file_path = os.path.join(input_dir, file) 
	goes_16 = ts.TimeSeries(input_file_path, concatenate=True)
	df = goes_16.to_dataframe()
	df = df[(df["xrsa_quality"] == 0) & (df["xrsb_quality"] == 0)]
	goes_16 = ts.TimeSeries(df, goes_16.meta, goes_16.units)	
	fig, ax = plt.subplots()
	goes_16.plot(axes=ax, columns=["xrsb"])

	# save png
	output_file_path = os.path.join(output_dir, file.replace(".nc", ".png"))
	plt.savefig(output_file_path, bbox_inches='tight')
	print(f"Processed and saved: {output_file_path}")
	plt.clf()

