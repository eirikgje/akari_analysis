import h5py
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
import healpy
from akari_defaults import nside, detector, pixel_flags, status_flags, reset_flag


def get_flagged_tods(data_file, pixel_flags=pixel_flags,
                     status_flags=status_flags, nside=nside,
                     detector=detector):

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
    reset_pixflag = ~np.array(hdf_file['FIS_OBS']['PIX_FLAG'])[:, detector-1, reset_flag]
    flags = np.ones(len(adu), dtype=bool)
    for pixflag in pixel_flags:
        flags &= ~np.array(hdf_file['FIS_OBS']['PIX_FLAG'])[:, detector-1, pixflag]
    for statusflag in status_flags:
        flags &= ~np.array(hdf_file['FIS_OBS']['STATUS'])[:, statusflag]
#    reset_pixflag[1:] = reset_pixflag[:-1] | reset_pixflag[1:]
#    lamps_on = status_flags[:, 16] | status_flags[:, 17] | status_flags[:, 18]
#    postprocess_flags &= ~reset_pixflag
#    postprocess_flags &= ~lamps_on
    base_data = {}
    base_data['flags'] = flags

    cds['adu'] = adu[is_cds]
    if len(cds['adu']) > 0:
        cds['flags'] = flags[is_cds]
        cds['ra'] = ra[is_cds]
        cds['dec'] = dec[is_cds]
        cds['coordinates'] = SkyCoord(ra=cds['ra']*u.deg, dec=cds['dec']*u.deg,
                                      frame='icrs')
        cds['aftime'] = aftime[is_cds]
        cds['lon'], cds['lat'] = cds['coordinates'].galactic.l.value, cds['coordinates'].galactic.b.value
        cds['pixs'] = healpy.ang2pix(nside, cds['lon'], cds['lat'], lonlat=True)
        cds['adu_flagged'] = cds['adu'][cds['flags']]
        cds['pixs_flagged'] = cds['pixs'][cds['flags']]

#    coadd['adu'] = adu[is_coadd] / 6 # Is this right?
    coadd['adu'] = adu[is_coadd]
    if len(coadd['adu']) > 0:
        coadd['flags'] = flags[is_coadd]
        coadd['ra'] = ra[is_coadd]
        coadd['dec'] = dec[is_coadd]
        coadd['coordinates'] = SkyCoord(ra=coadd['ra']*u.deg, dec=coadd['dec']*u.deg, frame='icrs')
        coadd['aftime'] = aftime[is_coadd]
        coadd['lon'], coadd['lat'] = coadd['coordinates'].galactic.l.value, coadd['coordinates'].galactic.b.value
        coadd['pixs'] = healpy.ang2pix(nside, coadd['lon'], coadd['lat'], lonlat=True)
        coadd['adu_flagged'] = coadd['adu'][coadd['flags']]
        coadd['pixs_flagged'] = coadd['pixs'][coadd['flags']]
        coadd['reset_flag'] = reset_pixflag[is_coadd]
        coadd['reset_flag_flagged'] = coadd['reset_flag'][coadd['flags']]

    return cds, coadd, base_data
