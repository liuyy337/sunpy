#!/home/potato/miniforge3/envs/sunpy/bin/python

import astropy.units as u
import gc
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import sunpy.map
import time
from astropy.coordinates import SkyCoord
from astropy.io import fits
from multiprocessing import Pool, cpu_count
from sunpy.timeseries import TimeSeries

dir_0094 = "/media/potato/solar/flare/20240511/mgn/aia_mgn_small/0094/"
dir_0131 = "/media/potato/solar/flare/20240511/mgn/aia_mgn_small/0131/"
dir_0171 = "/media/potato/solar/flare/20240511/mgn/aia_mgn_small/0171/"
dir_0193 = "/media/potato/solar/flare/20240511/mgn/aia_mgn_small/0193/"
dir_0211 = "/media/potato/solar/flare/20240511/mgn/aia_mgn_small/0211/"
dir_0304 = "/media/potato/solar/flare/20240511/mgn/aia_mgn_small/0304/"
dir_0335 = "/media/potato/solar/flare/20240511/mgn/aia_mgn_small/0335/"
dir_1600 = "/media/potato/solar/flare/20240511/cut_aia/cutted_aia_small/1600/"
output_dir = "figure/combined_aia_small/"

def setup(ax, time_str):
    ax.set_title('')
    ax.text(0.02, 0.98, time_str, transform=ax.transAxes, 
         color='white', fontsize=10, verticalalignment='top')
    ax.tick_params(axis='both', which='major', labelsize=10)
    lon, lat = ax.coords
    lon.set_axislabel('X (arcsec)', fontsize=10)
    lat.set_axislabel('Y (arcsec)', fontsize=10)
    lon.set_ticks(spacing=50 * u.arcsec, color='white', size=6)
    lat.set_ticks(spacing=50 * u.arcsec, color='white', size=6)
    lon.set_format_unit(u.arcsec, show_decimal_unit=False)
    lat.set_format_unit(u.arcsec, show_decimal_unit=False)
    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    lat.set_minor_frequency(5)
    ax.coords.grid(False)

def combined_plot(args):
    i, file_0094, file_0131, file_0171, file_0193, file_0211, file_0304, file_0335, file_1600, output_path = args
    
    # read aia maps
    map_0094 = sunpy.map.Map(file_0094)
    map_0131 = sunpy.map.Map(file_0131)
    map_0171 = sunpy.map.Map(file_0171)
    map_0193 = sunpy.map.Map(file_0193)
    map_0211 = sunpy.map.Map(file_0211)
    map_0304 = sunpy.map.Map(file_0304)
    map_0335 = sunpy.map.Map(file_0335)
    map_1600 = sunpy.map.Map(file_1600)

    fig = plt.figure(figsize=(16, 8), dpi=100, facecolor='black')
    gs = gridspec.GridSpec(60, 120, wspace=0, hspace=1)

    # axes for aia
    ax1 = fig.add_subplot(gs[0:30, 0:30], projection=map_0131)
    time_str = "AIA 131 \u00C5 " + map_0131.date.strftime('%Y-%m-%d %H:%M:%S')
    vmin = np.max(map_0131.data) - 1.05
    vmax = np.max(map_0131.data) - 0.15
    map_0131.plot(axes=ax1, autoalign=True, norm=colors.Normalize(vmin, vmax))
    setup(ax1, time_str)
    ax1.coords[0].set_ticklabel_visible(False)
   
    ax2 = fig.add_subplot(gs[0:30, 30:60], projection=map_0094)
    time_str = "AIA 94 \u00C5 " + map_0094.date.strftime('%Y-%m-%d %H:%M:%S')
    vmin = np.max(map_0094.data) - 1.05
    vmax = np.max(map_0094.data) - 0.1
    map_0094.plot(axes=ax2, autoalign=True, norm=colors.Normalize(vmin, vmax))
    setup(ax2, time_str)
    ax2.coords[0].set_ticklabel_visible(False)
    ax2.coords[1].set_ticklabel_visible(False)
    ax2.set_xlabel('')
    ax2.set_ylabel('')

    ax3 = fig.add_subplot(gs[0:30, 60:90], projection=map_0335)
    time_str = "AIA 335 \u00C5 " + map_0335.date.strftime('%Y-%m-%d %H:%M:%S')
    vmin = np.max(map_0335.data) - 1.05
    vmax = np.max(map_0335.data) - 0.15
    map_0335.plot(axes=ax3, autoalign=True, norm=colors.Normalize(vmin, vmax))
    setup(ax3, time_str)
    ax3.coords[0].set_ticklabel_visible(False)
    ax3.coords[1].set_ticklabel_visible(False)
    ax3.set_xlabel('')
    ax3.set_ylabel('')

    ax4 = fig.add_subplot(gs[0:30, 90:120], projection=map_0211)
    time_str = "AIA 211 \u00C5 " + map_0211.date.strftime('%Y-%m-%d %H:%M:%S')
    vmin = np.max(map_0211.data) - 1.05
    vmax = np.max(map_0211.data) - 0.1
    map_0211.plot(axes=ax4, autoalign=True, norm=colors.Normalize(vmin, vmax))
    setup(ax4, time_str)
    ax4.coords[0].set_ticklabel_visible(False)
    ax4.coords[1].set_ticklabel_visible(False)
    ax4.set_xlabel('')
    ax4.set_ylabel('')

    ax5 = fig.add_subplot(gs[30:60, 0:30], projection=map_0193)
    time_str = "AIA 193 \u00C5 " + map_0193.date.strftime('%Y-%m-%d %H:%M:%S')
    vmin = np.max(map_0193.data) - 1.05
    vmax = np.max(map_0193.data) - 0.1
    map_0193.plot(axes=ax5, autoalign=True, norm=colors.Normalize(vmin, vmax))
    setup(ax5, time_str)

    ax6 = fig.add_subplot(gs[30:60, 30:60], projection=map_0171)
    time_str = "AIA 171 \u00C5 " + map_0211.date.strftime('%Y-%m-%d %H:%M:%S')
    vmin = np.max(map_0171.data) - 1.05
    vmax = np.max(map_0171.data) - 0.1
    map_0171.plot(axes=ax6, autoalign=True, norm=colors.Normalize(vmin, vmax))
    setup(ax6, time_str)
    ax6.coords[1].set_ticklabel_visible(False)
    ax6.set_ylabel('')
   
    ax7 = fig.add_subplot(gs[30:60, 60:90], projection=map_0304)
    time_str = "AIA 304 \u00C5 " + map_0304.date.strftime('%Y-%m-%d %H:%M:%S')
    vmin = np.max(map_0304.data) - 1.05
    vmax = np.max(map_0304.data) - 0.1
    map_0304.plot(axes=ax7, autoalign=True, norm=colors.Normalize(vmin, vmax))
    setup(ax7, time_str)
    ax7.coords[1].set_ticklabel_visible(False)
    ax7.set_ylabel('')
    
    ax8 = fig.add_subplot(gs[30:60, 90:120], projection=map_1600)
    map_1600 = map_1600 / map_1600.exposure_time
    time_str = "AIA 1600 \u00C5 " + map_1600.date.strftime('%Y-%m-%d %H:%M:%S')
    # print(np.max(map_1600.data), np.max(map_1600.data))
    vmin = 10
    vmax = 1000
    map_1600.plot(axes=ax8, autoalign=True, norm=colors.LogNorm(vmin, vmax))
    setup(ax8, time_str)
    ax8.coords[1].set_ticklabel_visible(False)
    ax8.set_ylabel('')
 
    # save the png
    plt.savefig(output_path, transparent=False, facecolor='black', bbox_inches='tight')
    plt.close()

def main():
    plt.style.use('dark_background')
    os.makedirs(output_dir, exist_ok=True)

    # read all files
    files_0094 = sorted([dir_0094 + f for f in os.listdir(dir_0094) if f.endswith('.fits')])
    files_0131 = sorted([dir_0131 + f for f in os.listdir(dir_0131) if f.endswith('.fits')])
    files_0171 = sorted([dir_0171 + f for f in os.listdir(dir_0171) if f.endswith('.fits')])
    files_0193 = sorted([dir_0193 + f for f in os.listdir(dir_0193) if f.endswith('.fits')])
    files_0211 = sorted([dir_0211 + f for f in os.listdir(dir_0211) if f.endswith('.fits')])
    files_0304 = sorted([dir_0304 + f for f in os.listdir(dir_0304) if f.endswith('.fits')])
    files_0335 = sorted([dir_0335 + f for f in os.listdir(dir_0335) if f.endswith('.fits')])
    files_1600 = sorted([dir_1600 + f for f in os.listdir(dir_1600) if f.endswith('.fits')])

    args = []
    for i in range(600):
        args.append((i,
            files_0094[min(2 * i, len(files_0094)-1)],
            files_0131[min(2 * i, len(files_0131)-1)],
            files_0171[min(2 * i, len(files_0171)-1)],
            files_0193[min(2 * i, len(files_0193)-1)],
            files_0211[min(2 * i, len(files_0211)-1)],
            files_0304[min(2 * i, len(files_0304)-1)],
            files_0335[min(2 * i, len(files_0335)-1)],
            files_1600[min(i, len(files_1600)-1)],
            os.path.join(output_dir, f"{i+1:04}.png")
        ))
    with Pool(cpu_count()) as pool:
        pool.map(combined_plot, args)
    gc.collect()

if __name__ == "__main__":
    start_time = time.time()
    main()
    print(f"Total processing time: {time.time()-start_time:.2f} seconds")


