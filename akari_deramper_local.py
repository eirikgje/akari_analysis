import h5py
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import healpy
import glob
from multiprocessing import Pool
import matplotlib.pyplot as plt
from akari_data_processing import process_file
from akari_defaults import num_processes, npix, instrument, detector, nside


if __name__ == '__main__':
    #    dirs = ['2006_04', '2006_05', '2006_06', '2006_07', '2006_08', '2006_09',
#            '2006_10', '2006_11']

    tot_map = np.zeros(npix)
    tot_hitmap = np.zeros(npix)
    data_file_list = glob.glob(f'/home/eirik/data/akari_analysis/akari_TSD/2006_04/FIS_{instrument}_20060415*hdf5')

    total_datalen = 0
    used_datalen = 0
    with Pool(processes=num_processes) as pool:
        reses = pool.map(process_file, data_file_list)
        for res in reses:
            cds, coadd, base_data = res
            total_datalen += len(base_data['postprocess_flags'])
            used_datalen += np.sum(base_data['postprocess_flags'])
            if len(cds['adu']) > 0:
                tot_map += np.bincount(cds['pixs_postprocess'], weights=cds['adu_postprocess'], minlength=npix)
                tot_hitmap += np.bincount(cds['pixs_postprocess'], minlength=npix)
            if len(coadd['adu']) > 0:
                tot_map += np.bincount(coadd['pixs_postprocess'], weights=coadd['deramped_data_postprocess'], minlength=npix)
                tot_hitmap += np.bincount(coadd['pixs_postprocess'], minlength=npix)
    print('Total length of data:', total_datalen)
    print('Used length of data:', used_datalen)
    print('Used fraction:', used_datalen/total_datalen)
    nonzero_mask = np.where(tot_hitmap != 0)
    tot_map[nonzero_mask] = tot_map[nonzero_mask] / tot_hitmap[nonzero_mask]
    healpy.write_map(f'freq_map_{instrument}_det{detector+1}_nside{nside}.fits', tot_map, overwrite=True)
    healpy.write_map(f'hitmap_{instrument}_det{detector+1}_nside{nside}.fits', tot_hitmap, overwrite=True)
    healpy.mollview(tot_map, min=0, max=200)
    plt.savefig(f'freq_map_{instrument}_det{detector}_nside{nside}.png', bbox_inches='tight')
    plt.clf()
    healpy.mollview(tot_hitmap)
    plt.savefig(f'hitmap_{instrument}_det{detector}_nside{nside}.png', bbox_inches='tight')
