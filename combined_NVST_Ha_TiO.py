# import re
# import logging
# from datetime import datetime, timedelta
# import gc
# import matplotlib.colors as colors
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# import numpy as np
# import os
# import time
# import astropy.units as u
# from astropy.io import fits
# from multiprocessing import Pool, cpu_count, freeze_support
# import sunpy.map
# from sunpy.coordinates import frames
# from astropy.coordinates import SkyCoord

# # 目录配置 - 使用绝对路径
# BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# dir_ha_c = os.path.join(BASE_DIR, "paded/ha_c/")
# dir_ha_b = os.path.join(BASE_DIR, "paded/ha_b/")
# dir_ha_r = os.path.join(BASE_DIR, "paded/ha_r/")
# dir_tio = os.path.join(BASE_DIR, "paded/tio/")
# output_dir = os.path.join(BASE_DIR, "figure/NVST/Ha_TiO_combined/")
# ref_path = os.path.join(BASE_DIR, "paded/ha_c/Ha_13468_2023-10-23T02_30_02.552_level1+_paded.fits")

# # 配置日志记录
# logging.basicConfig(
#     level=logging.INFO,
#     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
#     handlers=[
#         logging.FileHandler('combined_processing.log', mode='w'),
#         logging.StreamHandler()
#     ]
# )
# logger = logging.getLogger(__name__)

# # 归一化参数配置
# NORM_PARAMS = {
#     'ha': {'vmin': 150000, 'vmax': 950000},
#     'tio': {'vmin': 15000000, 'vmax': 58000000}
# }

# # 波长配置 (单位: 埃)
# WAVELENGTHS = {
#     'ha': 6563,
#     'tio': 7058
# }

# def extract_timestamp_from_filename(filename):
#     """从文件名中提取时间戳"""
#     try:
#         basename = os.path.basename(filename)
#         pattern = r"(\d{4})-(\d{2})-(\d{2})[T _](\d{2})[_-](\d{2})[_-](\d{2})\.?(\d*)"
#         time_match = re.search(pattern, basename)
        
#         if time_match:
#             year, month, day = time_match.group(1), time_match.group(2), time_match.group(3)
#             hour, minute, second = time_match.group(4), time_match.group(5), time_match.group(6)
#             millisec = time_match.group(7).ljust(6, '0')[:6] if time_match.group(7) else "000000"
            
#             time_str = f"{year}-{month}-{day} {hour}:{minute}:{second}.{millisec}"
#             return datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S.%f")
        
#         logger.warning(f"文件名格式不匹配: {basename}")
#         return None
#     except Exception as e:
#         logger.error(f"提取时间戳错误: {filename} - {e}")
#         return None

# def setup(ax, title):
#     """设置图表坐标轴和标签"""
#     ax.set_title('')
#     ax.text(0.02, 0.98, title, transform=ax.transAxes, 
#             color='white', fontsize=8, verticalalignment='top')
#     ax.tick_params(axis='both', which='major', labelsize=8)
    
#     if hasattr(ax, 'coords'):
#         lon, lat = ax.coords
#         lon.set_axislabel('X (arcsec)', fontsize=8)
#         lat.set_axislabel('Y (arcsec)', fontsize=8)
#         lon.set_ticks(spacing=50 * u.arcsec, color='white', size=3)
#         lat.set_ticks(spacing=50 * u.arcsec, color='white', size=3)
#         lon.set_format_unit(u.arcsec, show_decimal_unit=False)
#         lat.set_format_unit(u.arcsec, show_decimal_unit=False)
#         lon.display_minor_ticks(True)
#         lat.display_minor_ticks(True)
#         lat.set_minor_frequency(5)
#         ax.coords.grid(False)

# def create_sunmap(file_path, is_tio=False):
#     """创建地图对象并修复元数据"""
#     try:
#         with fits.open(file_path, memmap=False) as hdul:
#             header = hdul[0].header.copy()
#             data = hdul[0].data.astype(np.float32)
            
#             # 读取曝光时间(单位:ms)并转换为秒
#             exposure_time = header.get('EXPTIME', 20 if not is_tio else 0.9)  # 默认值单位ms
#             if isinstance(exposure_time, str):
#                 try:
#                     exposure_time = float(exposure_time)
#                 except ValueError:
#                     exposure_time = 20 if not is_tio else 0.9
            
#             # 将毫秒转换为秒
#             exposure_time_sec = exposure_time / 1000.0
            
#             # 设置必要的头信息
#             header['CTYPE1'] = 'HPLN-TAN'
#             header['CTYPE2'] = 'HPLT-TAN'
#             header['CUNIT1'] = 'arcsec'
#             header['CUNIT2'] = 'arcsec'
            
#             # 设置关键参数
#             header['CDELT1'] = 0.165 if not is_tio else 0.052  # arcsec/pixel
#             header['CDELT2'] = 0.165 if not is_tio else 0.052  # arcsec/pixel
#             header['CRPIX1'] = data.shape[1] // 2
#             header['CRPIX2'] = data.shape[0] // 2
#             header['CRVAL1'] = 0.0
#             header['CRVAL2'] = 0.0
#             header['EXPTIME'] = exposure_time_sec  # 存储为秒
#             header['RSUN_OBS'] = 6.96e8
#             header['DSUN_OBS'] = 1.496e11
#             header['HGLN_OBS'] = 0.0
#             header['HGLT_OBS'] = 0.0
#             header['CRLN_OBS'] = 0.0
#             header['CRLT_OBS'] = 0.0
#             header['TELESCOP'] = 'NVST'
#             header['INSTRUME'] = 'HRIS'
#             header['WAVELNTH'] = WAVELENGTHS['tio'] if is_tio else WAVELENGTHS['ha']
#             header['WAVEUNIT'] = 'angstrom'
#             header['PASSBAND'] = 10.0
#             header['OFFBAND'] = 0.0
            
#             return sunpy.map.Map(data, header)
#     except Exception as e:
#         logger.error(f"创建地图失败: {file_path} - {e}")
#         return None

# def apply_exposure_correction(smap):
#     """应用曝光时间校正"""
#     if smap is None:
#         return None
    
#     try:
#         exposure_time_sec = smap.meta.get('EXPTIME', 0.02)  # 单位已经是秒
#         if exposure_time_sec <= 0:
#             logger.warning(f"无效的曝光时间: {exposure_time_sec}s, 使用默认值0.02s")
#             exposure_time_sec = 0.02
        
#         corrected_data = smap.data / exposure_time_sec
#         new_meta = smap.meta.copy()
#         new_meta['EXPTIME'] = 1.0  # 校正后的曝光时间设为1.0秒
        
#         logger.info(f"应用曝光校正: 曝光时间={exposure_time_sec*1000:.1f}ms (={exposure_time_sec:.4f}s), "
#                    f"最大值校正前={np.max(smap.data):.1f}, 校正后={np.max(corrected_data):.1f}")
#         return sunpy.map.Map(corrected_data, new_meta)
#     except Exception as e:
#         logger.error(f"曝光校正失败: {e}")
#         return smap

# def simple_match_files(files_c, files_b, files_r, files_tio):
#     """简化文件匹配方法"""
#     matched = []
#     tolerance = 30  # 30秒的最大时间差
    
#     times_c = [extract_timestamp_from_filename(f) for f in files_c]
#     times_b = [extract_timestamp_from_filename(f) for f in files_b]
#     times_r = [extract_timestamp_from_filename(f) for f in files_r]
#     times_tio = [extract_timestamp_from_filename(f) for f in files_tio]
    
#     for i, ts_c in enumerate(times_c):
#         if ts_c is None:
#             continue
            
#         closest_b = min(times_b, key=lambda ts: abs((ts or datetime.max) - ts_c)) if times_b else None
#         if closest_b is None or abs((closest_b - ts_c).total_seconds()) > tolerance:
#             continue
            
#         closest_r = min(times_r, key=lambda ts: abs((ts or datetime.max) - ts_c)) if times_r else None
#         if closest_r is None or abs((closest_r - ts_c).total_seconds()) > tolerance:
#             continue
            
#         closest_tio = min(times_tio, key=lambda ts: abs((ts or datetime.max) - ts_c)) if times_tio else None
#         if closest_tio is None or abs((closest_tio - ts_c).total_seconds()) > tolerance:
#             continue
            
#         file_b = files_b[times_b.index(closest_b)]
#         file_r = files_r[times_r.index(closest_r)]
#         file_tio = files_tio[times_tio.index(closest_tio)]
        
#         matched.append((files_c[i], file_b, file_r, file_tio))
    
#     matched.sort(key=lambda x: extract_timestamp_from_filename(x[0]) or datetime.min)
#     return matched

# def combined_plot(args):
#     """组合并绘制多通道图像"""
#     i, file_ha_c, file_ha_b, file_ha_r, file_tio, output_path = args
    
#     fig = None
#     try:
#         logger.info(f"开始处理 #{i}: {os.path.basename(file_ha_c)}")
        
#         ts_c = extract_timestamp_from_filename(file_ha_c)
#         time_str = ts_c.strftime('%H:%M:%S') if ts_c else "Unknown"
        
#         ref_map = None
#         if os.path.exists(ref_path):
#             ref_map = create_sunmap(ref_path)
        
#         # 创建地图并应用曝光校正
#         map_ha_c = apply_exposure_correction(create_sunmap(file_ha_c))
#         map_ha_b = apply_exposure_correction(create_sunmap(file_ha_b))
#         map_ha_r = apply_exposure_correction(create_sunmap(file_ha_r))
#         map_tio = apply_exposure_correction(create_sunmap(file_tio, is_tio=True))
        
#         if None in [map_ha_c, map_ha_b, map_ha_r, map_tio]:
#             logger.error("无法加载一个或多个通道地图")
#             return False
            
#         if ref_map:
#             try:
#                 def align_to_ref(smap):
#                     for key in ['CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2']:
#                         if key in ref_map.meta:
#                             smap.meta[key] = ref_map.meta[key]
#                     return smap
                
#                 map_ha_c = align_to_ref(map_ha_c)
#                 map_ha_b = align_to_ref(map_ha_b)
#                 map_ha_r = align_to_ref(map_ha_r)
#                 map_tio = align_to_ref(map_tio)
#             except Exception as e:
#                 logger.error(f"坐标系对齐失败: {e}")
        
#         fig = plt.figure(figsize=(10, 10), facecolor="black", dpi=200)
#         gs = gridspec.GridSpec(2, 2, wspace=0.1, hspace=0.15)
        
#         # Hα-0.4Å
#         ax1 = fig.add_subplot(gs[0, 0], projection=map_ha_b)
#         map_ha_b.plot(axes=ax1, norm=colors.Normalize(**NORM_PARAMS['ha']))
#         setup(ax1, f"Hα-0.4Å {time_str}")
        
#         # Hα
#         ax2 = fig.add_subplot(gs[0, 1], projection=map_ha_c)
#         map_ha_c.plot(axes=ax2, norm=colors.Normalize(**NORM_PARAMS['ha']))
#         setup(ax2, f"Hα {time_str}")
#         ax2.coords[1].set_ticklabel_visible(False)
        
#         # Hα+0.4Å
#         ax3 = fig.add_subplot(gs[1, 0], projection=map_ha_r)
#         map_ha_r.plot(axes=ax3, norm=colors.Normalize(**NORM_PARAMS['ha']))
#         setup(ax3, f"Hα+0.4Å {time_str}")
        
#         # TiO
#         ax4 = fig.add_subplot(gs[1, 1], projection=map_tio)
#         map_tio.plot(axes=ax4, norm=colors.Normalize(**NORM_PARAMS['tio']))
#         setup(ax4, f"TiO {time_str}")
#         ax4.coords[1].set_ticklabel_visible(False)
        
#         plt.savefig(output_path, dpi=200, facecolor='black', bbox_inches='tight')
#         logger.info(f"成功保存: {output_path}")
#         return True
    
#     except Exception as e:
#         logger.error(f"处理文件失败: {os.path.basename(file_ha_c)} - {str(e)}")
#         return False
#     finally:
#         if fig:
#             plt.close(fig)
#         gc.collect()

# def validate_paths():
#     """验证所有路径是否存在"""
#     paths_to_check = [dir_ha_c, dir_ha_b, dir_ha_r, dir_tio, os.path.dirname(output_dir)]
    
#     for path in paths_to_check:
#         if not os.path.exists(path):
#             try:
#                 if path == os.path.dirname(output_dir):
#                     os.makedirs(path, exist_ok=True)
#                 else:
#                     logger.error(f"路径不存在: {path}")
#                     return False
#             except Exception as e:
#                 logger.error(f"无法创建输出目录: {output_dir} - {e}")
#                 return False
    
#     logger.info("所有路径验证通过")
#     return True

# def main():
#     """主函数"""
#     try:
#         start_time = time.perf_counter()
#         plt.style.use('dark_background')
        
#         logger.info("程序启动")
        
#         if not validate_paths():
#             return
            
#         os.makedirs(output_dir, exist_ok=True)
        
#         def get_fits_files(directory):
#             return [os.path.join(directory, f) for f in os.listdir(directory) 
#                     if f.endswith('.fits') or f.endswith('.fits.gz')]
        
#         files_ha_c = sorted(get_fits_files(dir_ha_c))
#         files_ha_b = sorted(get_fits_files(dir_ha_b))
#         files_ha_r = sorted(get_fits_files(dir_ha_r))
#         files_tio = sorted(get_fits_files(dir_tio))
        
#         logger.info(f"中心线文件: {len(files_ha_c)}")
#         logger.info(f"蓝移文件: {len(files_ha_b)}")
#         logger.info(f"红移文件: {len(files_ha_r)}")
#         logger.info(f"TiO文件: {len(files_tio)}")
        
#         matched_files = simple_match_files(files_ha_c, files_ha_b, files_ha_r, files_tio)
#         logger.info(f"匹配文件组: {len(matched_files)}")
        
#         if not matched_files:
#             logger.error("错误：找不到匹配的文件组")
#             return
        
#         args = []
#         for i, (file_c, file_b, file_r, file_tio) in enumerate(matched_files):
#             ts = extract_timestamp_from_filename(file_c) or datetime.now()
#             filename = f"NVST_Combined_{ts.strftime('%Y%m%d_%H%M%S')}.png"
#             output_path = os.path.join(output_dir, filename)
#             args.append((i, file_c, file_b, file_r, file_tio, output_path))
        
#         if len(args) == 1:
#             result = combined_plot(args[0])
#             logger.info(f"处理完成，结果: {'成功' if result else '失败'}")
#         else:
#             num_cores = min(cpu_count(), 30, len(args))
#             logger.info(f"使用 {num_cores} 个核心处理 {len(args)} 个任务")
            
#             try:
#                 if os.name == 'nt':
#                     freeze_support()
                    
#                 with Pool(processes=num_cores) as pool:
#                     results = pool.imap(combined_plot, args)
#                     success_count = sum(1 for result in results if result)
#                     logger.info(f"成功处理 {success_count}/{len(args)} 个图像")
#             except Exception as e:
#                 logger.error(f"并行处理失败: {e}")
        
#         end_time = time.perf_counter()
#         total_time = end_time - start_time
#         logger.info(f"总处理时间: {total_time:.2f} 秒")
#         print(f"总处理时间: {total_time:.2f} 秒")
    
#     except Exception as e:
#         logger.error(f"主函数错误: {e}")
#     finally:
#         logger.info("程序结束")

# if __name__ == "__main__":
#     main()


import re
import logging
from datetime import datetime, timedelta
import gc
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import time
import astropy.units as u
from astropy.io import fits
from multiprocessing import Pool, cpu_count, freeze_support
import sunpy.map
from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord

# 坐标轴范围配置
X_MIN = 94.248  # 角秒
X_MAX = 263.088  # 角秒
Y_MIN = -347.584  # 角秒
Y_MAX = -178.744  # 角秒

# 目录配置 - 使用绝对路径
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
dir_ha_c = os.path.join(BASE_DIR, "paded/ha_c/")
dir_ha_b = os.path.join(BASE_DIR, "paded/ha_b/")
dir_ha_r = os.path.join(BASE_DIR, "paded/ha_r/")
dir_tio = os.path.join(BASE_DIR, "paded/tio/")
output_dir = os.path.join(BASE_DIR, "figure/NVST/Ha_TiO_combined/")
ref_path = os.path.join(BASE_DIR, "paded/ha_c/Ha_13468_2023-10-23T02_30_02.552_level1+_paded.fits")

# 配置日志记录
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('combined_processing.log', mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# 归一化参数配置
NORM_PARAMS = {
    'ha': {'vmin': 150000, 'vmax': 950000},
    'tio': {'vmin': 15000000, 'vmax': 58000000}
}

# 波长配置 (单位: 埃)
WAVELENGTHS = {
    'ha': 6563,
    'tio': 7058
}

def extract_timestamp_from_filename(filename):
    """从文件名中提取时间戳"""
    try:
        basename = os.path.basename(filename)
        pattern = r"(\d{4})-(\d{2})-(\d{2})[T _](\d{2})[_-](\d{2})[_-](\d{2})\.?(\d*)"
        time_match = re.search(pattern, basename)
        
        if time_match:
            year, month, day = time_match.group(1), time_match.group(2), time_match.group(3)
            hour, minute, second = time_match.group(4), time_match.group(5), time_match.group(6)
            millisec = time_match.group(7).ljust(3, '0')[:3] if time_match.group(7) else "000"
            
            time_str = f"{year}-{month}-{day} {hour}:{minute}:{second}.{millisec}"
            return datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S.%f")
        
        logger.warning(f"文件名格式不匹配: {basename}")
        return None
    except Exception as e:
        logger.error(f"提取时间戳错误: {filename} - {e}")
        return None

def setup(ax, title):
    """设置图表坐标轴和标签"""
    ax.set_title('')
    ax.text(0.02, 0.98, title, transform=ax.transAxes, 
            color='white', fontsize=8, verticalalignment='top')
    ax.tick_params(axis='both', which='major', labelsize=8)
    
    if hasattr(ax, 'coords'):
        lon, lat = ax.coords
        lon.set_axislabel('X (arcsec)', fontsize=8)
        lat.set_axislabel('Y (arcsec)', fontsize=8)
        
        # 设置坐标轴范围
        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(Y_MIN, Y_MAX)
        
        # 设置刻度间隔
        lon.set_ticks(spacing=50 * u.arcsec, color='white', size=3)
        lat.set_ticks(spacing=50 * u.arcsec, color='white', size=3)
        lon.set_format_unit(u.arcsec, show_decimal_unit=False)
        lat.set_format_unit(u.arcsec, show_decimal_unit=False)
        lon.display_minor_ticks(True)
        lat.display_minor_ticks(True)
        lat.set_minor_frequency(5)
        ax.coords.grid(False)

def create_sunmap(file_path, is_tio=False):
    """创建地图对象并修复元数据"""
    try:
        with fits.open(file_path, memmap=False) as hdul:
            header = hdul[0].header.copy()
            data = hdul[0].data.astype(np.float32)
            
            # 检查数据是否为空
            if data.size == 0:
                logger.error(f"空数据文件: {file_path}")
                return None
            
            # 从文件名中提取时间戳并添加到header中
            timestamp = extract_timestamp_from_filename(file_path)
            if timestamp:
                header['DATE-OBS'] = timestamp.strftime("%Y-%m-%dT%H:%M:%S.%f")
            
            # 读取曝光时间(单位:ms)并转换为秒
            exposure_time = header.get('EXPTIME', 20 if not is_tio else 0.9)  # 默认值单位ms
            if isinstance(exposure_time, str):
                try:
                    exposure_time = float(exposure_time)
                except ValueError:
                    exposure_time = 20 if not is_tio else 0.9
            
            # 将毫秒转换为秒
            exposure_time_sec = exposure_time / 1000.0
            
            # 设置必要的头信息
            header['CTYPE1'] = 'HPLN-TAN'
            header['CTYPE2'] = 'HPLT-TAN'
            header['CUNIT1'] = 'arcsec'
            header['CUNIT2'] = 'arcsec'
            
            # 设置关键参数
            header['CDELT1'] = 0.165 if not is_tio else 0.052  # arcsec/pixel
            header['CDELT2'] = 0.165 if not is_tio else 0.052  # arcsec/pixel
            header['CRPIX1'] = data.shape[1] // 2
            header['CRPIX2'] = data.shape[0] // 2
            header['CRVAL1'] = 0.0
            header['CRVAL2'] = 0.0
            header['EXPTIME'] = exposure_time_sec  # 存储为秒
            header['RSUN_OBS'] = 6.96e8
            header['DSUN_OBS'] = 1.496e11
            header['HGLN_OBS'] = 0.0
            header['HGLT_OBS'] = 0.0
            header['CRLN_OBS'] = 0.0
            header['CRLT_OBS'] = 0.0
            header['TELESCOP'] = 'NVST'
            header['INSTRUME'] = 'HRIS'
            header['WAVELNTH'] = WAVELENGTHS['tio'] if is_tio else WAVELENGTHS['ha']
            header['WAVEUNIT'] = 'angstrom'
            header['PASSBAND'] = 10.0
            header['OFFBAND'] = 0.0
            
            return sunpy.map.Map(data, header)
    except Exception as e:
        logger.error(f"创建地图失败: {file_path} - {e}")
        return None

def apply_exposure_correction(smap):
    """应用曝光时间校正"""
    if smap is None:
        return None
    
    try:
        exposure_time_sec = smap.meta.get('EXPTIME', 0.02)  # 单位已经是秒
        if exposure_time_sec <= 0:
            logger.warning(f"无效的曝光时间: {exposure_time_sec}s, 使用默认值0.02s")
            exposure_time_sec = 0.02
        
        corrected_data = smap.data / exposure_time_sec
        new_meta = smap.meta.copy()
        new_meta['EXPTIME'] = 1.0  # 校正后的曝光时间设为1.0秒
        
        logger.info(f"应用曝光校正: 曝光时间={exposure_time_sec*1000:.1f}ms (={exposure_time_sec:.4f}s), "
                   f"最大值校正前={np.max(smap.data):.1f}, 校正后={np.max(corrected_data):.1f}")
        return sunpy.map.Map(corrected_data, new_meta)
    except Exception as e:
        logger.error(f"曝光校正失败: {e}")
        return smap

def simple_match_files(files_c, files_b, files_r, files_tio):
    """简化文件匹配方法"""
    matched = []
    tolerance = 30  # 30秒的最大时间差
    
    times_c = [extract_timestamp_from_filename(f) for f in files_c]
    times_b = [extract_timestamp_from_filename(f) for f in files_b]
    times_r = [extract_timestamp_from_filename(f) for f in files_r]
    times_tio = [extract_timestamp_from_filename(f) for f in files_tio]
    
    for i, ts_c in enumerate(times_c):
        if ts_c is None:
            continue
            
        closest_b = min(times_b, key=lambda ts: abs((ts or datetime.max) - ts_c)) if times_b else None
        if closest_b is None or abs((closest_b - ts_c).total_seconds()) > tolerance:
            continue
            
        closest_r = min(times_r, key=lambda ts: abs((ts or datetime.max) - ts_c)) if times_r else None
        if closest_r is None or abs((closest_r - ts_c).total_seconds()) > tolerance:
            continue
            
        closest_tio = min(times_tio, key=lambda ts: abs((ts or datetime.max) - ts_c)) if times_tio else None
        if closest_tio is None or abs((closest_tio - ts_c).total_seconds()) > tolerance:
            continue
            
        file_b = files_b[times_b.index(closest_b)]
        file_r = files_r[times_r.index(closest_r)]
        file_tio = files_tio[times_tio.index(closest_tio)]
        
        matched.append((files_c[i], file_b, file_r, file_tio))
    
    matched.sort(key=lambda x: extract_timestamp_from_filename(x[0]) or datetime.min)
    return matched

def combined_plot(args):
    """组合并绘制多通道图像"""
    i, file_ha_c, file_ha_b, file_ha_r, file_tio, output_path = args
    
    fig = None
    try:
        logger.info(f"开始处理 #{i}: {os.path.basename(file_ha_c)}")
        
        ts_c = extract_timestamp_from_filename(file_ha_c)
        time_str = ts_c.strftime('%H:%M:%S') if ts_c else "Unknown"
        
        ref_map = None
        if os.path.exists(ref_path):
            ref_map = create_sunmap(ref_path)
        
        # 创建地图并应用曝光校正
        map_ha_c = apply_exposure_correction(create_sunmap(file_ha_c))
        map_ha_b = apply_exposure_correction(create_sunmap(file_ha_b))
        map_ha_r = apply_exposure_correction(create_sunmap(file_ha_r))
        map_tio = apply_exposure_correction(create_sunmap(file_tio, is_tio=True))
        
        # 检查地图是否有效
        if None in [map_ha_c, map_ha_b, map_ha_r, map_tio]:
            logger.error("无法加载一个或多个通道地图")
            return False
            
        # 检查数据是否为空
        for name, smap in [('Hα中心', map_ha_c), ('Hα蓝移', map_ha_b), 
                          ('Hα红移', map_ha_r), ('TiO', map_tio)]:
            if smap.data.size == 0:
                logger.error(f"{name}通道数据为空")
                return False
        
        if ref_map:
            try:
                def align_to_ref(smap):
                    for key in ['CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2']:
                        if key in ref_map.meta:
                            smap.meta[key] = ref_map.meta[key]
                    return smap
                
                map_ha_c = align_to_ref(map_ha_c)
                map_ha_b = align_to_ref(map_ha_b)
                map_ha_r = align_to_ref(map_ha_r)
                map_tio = align_to_ref(map_tio)
            except Exception as e:
                logger.error(f"坐标系对齐失败: {e}")
        
        fig = plt.figure(figsize=(10, 10), facecolor="black", dpi=200)
        gs = gridspec.GridSpec(2, 2, wspace=0.1, hspace=0.15)
        
        # Hα-0.4Å
        ax1 = fig.add_subplot(gs[0, 0], projection=map_ha_b)
        map_ha_b.plot(axes=ax1, norm=colors.Normalize(**NORM_PARAMS['ha']))
        setup(ax1, f"Hα-0.4Å {time_str}")
        
        # Hα
        ax2 = fig.add_subplot(gs[0, 1], projection=map_ha_c)
        map_ha_c.plot(axes=ax2, norm=colors.Normalize(**NORM_PARAMS['ha']))
        setup(ax2, f"Hα {time_str}")
        ax2.coords[1].set_ticklabel_visible(False)
        
        # Hα+0.4Å
        ax3 = fig.add_subplot(gs[1, 0], projection=map_ha_r)
        map_ha_r.plot(axes=ax3, norm=colors.Normalize(**NORM_PARAMS['ha']))
        setup(ax3, f"Hα+0.4Å {time_str}")
        
        # TiO
        ax4 = fig.add_subplot(gs[1, 1], projection=map_tio)
        map_tio.plot(axes=ax4, norm=colors.Normalize(**NORM_PARAMS['tio']))
        setup(ax4, f"TiO {time_str}")
        ax4.coords[1].set_ticklabel_visible(False)
        
        plt.savefig(output_path, dpi=200, facecolor='black', bbox_inches='tight')
        logger.info(f"成功保存: {output_path}")
        return True
    
    except Exception as e:
        logger.error(f"处理文件失败: {os.path.basename(file_ha_c)} - {str(e)}")
        return False
    finally:
        if fig:
            plt.close(fig)
        gc.collect()

def validate_paths():
    """验证所有路径是否存在"""
    paths_to_check = [dir_ha_c, dir_ha_b, dir_ha_r, dir_tio, os.path.dirname(output_dir)]
    
    for path in paths_to_check:
        if not os.path.exists(path):
            try:
                if path == os.path.dirname(output_dir):
                    os.makedirs(path, exist_ok=True)
                else:
                    logger.error(f"路径不存在: {path}")
                    return False
            except Exception as e:
                logger.error(f"无法创建输出目录: {output_dir} - {e}")
                return False
    
    logger.info("所有路径验证通过")
    return True

def main():
    """主函数"""
    try:
        start_time = time.perf_counter()
        plt.style.use('dark_background')
        
        logger.info("程序启动")
        
        if not validate_paths():
            return
            
        os.makedirs(output_dir, exist_ok=True)
        
        def get_fits_files(directory):
            return [os.path.join(directory, f) for f in os.listdir(directory) 
                    if f.endswith('.fits') or f.endswith('.fits.gz')]
        
        files_ha_c = sorted(get_fits_files(dir_ha_c))
        files_ha_b = sorted(get_fits_files(dir_ha_b))
        files_ha_r = sorted(get_fits_files(dir_ha_r))
        files_tio = sorted(get_fits_files(dir_tio))
        
        logger.info(f"中心线文件: {len(files_ha_c)}")
        logger.info(f"蓝移文件: {len(files_ha_b)}")
        logger.info(f"红移文件: {len(files_ha_r)}")
        logger.info(f"TiO文件: {len(files_tio)}")
        
        matched_files = simple_match_files(files_ha_c, files_ha_b, files_ha_r, files_tio)
        logger.info(f"匹配文件组: {len(matched_files)}")
        
        if not matched_files:
            logger.error("错误：找不到匹配的文件组")
            return
        
        args = []
        for i, (file_c, file_b, file_r, file_tio) in enumerate(matched_files):
            ts = extract_timestamp_from_filename(file_c) or datetime.now()
            filename = f"NVST_Combined_{ts.strftime('%Y%m%d_%H%M%S')}.png"
            output_path = os.path.join(output_dir, filename)
            args.append((i, file_c, file_b, file_r, file_tio, output_path))
        
        if len(args) == 1:
            result = combined_plot(args[0])
            logger.info(f"处理完成，结果: {'成功' if result else '失败'}")
        else:
            num_cores = min(cpu_count(), 30, len(args))
            logger.info(f"使用 {num_cores} 个核心处理 {len(args)} 个任务")
            
            try:
                if os.name == 'nt':
                    freeze_support()
                    
                with Pool(processes=num_cores) as pool:
                    results = pool.imap(combined_plot, args)
                    success_count = sum(1 for result in results if result)
                    logger.info(f"成功处理 {success_count}/{len(args)} 个图像")
            except Exception as e:
                logger.error(f"并行处理失败: {e}")
        
        end_time = time.perf_counter()
        total_time = end_time - start_time
        logger.info(f"总处理时间: {total_time:.2f} 秒")
        print(f"总处理时间: {total_time:.2f} 秒")
    
    except Exception as e:
        logger.error(f"主函数错误: {e}")
    finally:
        logger.info("程序结束")

if __name__ == "__main__":
    main()