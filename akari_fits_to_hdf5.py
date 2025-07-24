"""
This script simply rewrites the AKARI fits files into more readable hdf5 files
for easier later processing
"""

import fitsio
import h5py

extension_info = {
        1: {'name': 'FIS_OBS'},
        2: {'name': 'FIS_HK'},
        3: {'name': 'IRC_HK'},
        4: {'name': 'HK_2'},
        5: {'name': 'AOCU'},
        6: {'name': 'GADS'},
        7: {'name': 'PR'},
        8: {'name': 'SE'}
    }


def generate_filelist():
    infilelist = ['/home/eirik/data/akari_analysis/FIS_LW_20060429130000_gb.fits',
                  '/home/eirik/data/akari_analysis/FIS_LW_20060923000000_gb.fits']
    outfilelist = ['/home/eirik/data/akari_analysis/FIS_LW_20060429130000_gb_gzip9compression.hdf5',
                   '/home/eirik/data/akari_analysis/FIS_LW_20060923000000_gb_gzip9compression.hdf5']

    return infilelist, outfilelist


infilelist, outfilelist = generate_filelist()

for infile, outfile in zip(infilelist, outfilelist):
    fits = fitsio.FITS(infile)
    hdf_out = h5py.File(outfile, 'w')
    for extension, info in extension_info.items():
        curr_extension_group = hdf_out.create_group(info['name'])
        for colname in fits[extension].get_colnames():
            curr_extension_group.create_dataset(colname, compression="gzip",
                                                compression_opts=9,
                                                data=fits[extension][colname].read())
    hdf_out.close()
    fits.close()
