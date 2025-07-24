import h5py
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import healpy
import glob
from multiprocessing import Pool
import matplotlib.pyplot as plt

nside = 256

instrument = 'SW'
detector = 95
#data_files = [f'/home/eirik/data/akari_analysis/akari_TSD/2006_04/FIS_{instrument}_20060415090000_gb.hdf5']
#data_files = glob.glob(f'/home/eirik/data/akari_analysis/akari_TSD/2006_04/FIS_{instrument}_20060415*hdf5')
#print(len(data_files))
npix = 12 * nside ** 2
tot_map = np.zeros(npix)
tot_hitmap = np.zeros(npix)

#pixel_flags = [0, 5] # bad pixel, reset_pixflag
pixel_flags = [0, 5, 6] # bad_pixel, reset_pixflag, reset anomaly
status_flags = [16, 17, 18] # three lamps on

def process_file(data_file, return_last_timeline=False, npix=npix, pixel_flags=pixel_flags, status_flags=status_flags):

    hdf_file = h5py.File(data_file)
    packet_id = np.array(hdf_file['FIS_OBS']['PACKETID'])
    cds = {}
    coadd = {}

    ra = hdf_file['GADS']['RA']
    dec = hdf_file['GADS']['DEC']
    adu = np.array(hdf_file['FIS_OBS']['DET'])[:, detector-1]
    aftime = np.array(hdf_file['FIS_OBS']['AFTIME'])
    is_cds = np.where(packet_id == 65)
    is_coadd = np.where(packet_id == 66)
#    status_flags = np.array(hdf_file['FIS_OBS']['STATUS'])
#    bad_pixel = np.array(hdf_file['FIS_OBS']['PIX_FLAG'])[:, detector-1, 0]
#    reset_pixflag = np.array(hdf_file['FIS_OBS']['PIX_FLAG'])[:, detector-1, 5]
    postprocess_flags = np.ones(len(adu), dtype=bool)
    for pixflag in pixel_flags:
        postprocess_flags &= ~np.array(hdf_file['FIS_OBS']['PIX_FLAG'])[:, detector-1, pixflag]
    for statusflag in status_flags:
        postprocess_flags &= ~np.array(hdf_file['FIS_OBS']['STATUS'])[:, statusflag]
#    reset_pixflag[1:] = reset_pixflag[:-1] | reset_pixflag[1:]
#    lamps_on = status_flags[:, 16] | status_flags[:, 17] | status_flags[:, 18]
#    postprocess_flags &= ~reset_pixflag
#    postprocess_flags &= ~lamps_on
    base_data = {}
    base_data['postprocess_flags'] = postprocess_flags

    cds['adu'] = adu[is_cds]
    if len(cds['adu']) > 0:
        cds['postprocess_flags'] = postprocess_flags[is_cds]
        cds['ra'] = ra[is_cds]
        cds['dec'] = dec[is_cds]
        cds['coordinates'] = SkyCoord(ra=cds['ra']*u.deg, dec=cds['dec']*u.deg,
                                      frame='icrs')
        cds['aftime'] = aftime[is_cds]
        cds['lon'], cds['lat'] = cds['coordinates'].galactic.l.value, cds['coordinates'].galactic.b.value
        cds['pixs'] = healpy.ang2pix(nside, cds['lon'], cds['lat'], lonlat=True)
        cds['adu_postprocess'] = cds['adu'][cds['postprocess_flags']]
        cds['pixs_postprocess'] = cds['pixs'][cds['postprocess_flags']]

    coadd['adu'] = adu[is_coadd] / 6 # Is this right?
    if len(coadd['adu']) > 0:
        coadd['postprocess_flags'] = postprocess_flags[is_coadd][1:]
        coadd['deramped_data'] = (coadd['adu'][1:] - coadd['adu'][0:-1]) # This is going to be wrong when there is a cds gap
        coadd['ra'] = ra[is_coadd][1:]
        coadd['dec'] = dec[is_coadd][1:]
        coadd['coordinates'] = SkyCoord(ra=coadd['ra']*u.deg, dec=coadd['dec']*u.deg, frame='icrs')
        coadd['aftime'] = aftime[is_coadd][1:]
        coadd['lon'], coadd['lat'] = coadd['coordinates'].galactic.l.value, coadd['coordinates'].galactic.b.value
        coadd['pixs'] = healpy.ang2pix(nside, coadd['lon'], coadd['lat'], lonlat=True)
        coadd['deramped_data_postprocess'] = coadd['deramped_data'][coadd['postprocess_flags']]
        coadd['pixs_postprocess'] = coadd['pixs'][coadd['postprocess_flags']]

    return cds, coadd, base_data


if __name__ == '__main__':

    data_dir = '/mn/stornext/d23/cmbco/cg/AKARI/akari_TSD_hdf5/'
    data_file_list = glob.glob(f'{data_dir}/*/FIS_{instrument}*hdf5')

    total_datalen = 0
    used_datalen = 0
    with Pool(processes=8) as pool:
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
