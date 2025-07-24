from akari_data_loading import get_flagged_tods
from akari_defaults import pixel_flags, status_flags, detector
import healpy

def process_file(data_file, return_last_timeline=False,
                 pixel_flags=pixel_flags, status_flags=status_flags):

    cds, coadd, base_data = get_flagged_tods(data_file, pixel_flags, status_flags)
    if len(coadd['adu']) > 0:
        coadd['adu'] = coadd['adu'] / 6.0
        coadd['deramped_data'] = (coadd['adu'][1:] - coadd['adu'][0:-1]) # This is going to be wrong when there is a cds gap
#        coadd['postpocess_flags'] = coadd['postprocess_flags'][1:]
#        coadd['ra'] = coadd['ra'][1:]
#        coadd['dec'] = coadd['dec'][1:]
#        coadd['coordinates'] = coadd['coordinates'][1:]
#        coadd['aftime'] = coadd['aftime'][1:]
#        coadd['lon'], coadd['lat'] = coadd['lon'][1:], coadd['lat'][1:]
#        coadd['pixs'] = coadd['pixs'][1:]
        coadd['pixs_flagged'] = coadd['pixs'][1:][coadd['flags'][1:]]
        coadd['deramped_data_flagged'] = coadd['deramped_data'][coadd['flags'][1:]]
        if len(coadd['deramped_data_flagged']) == len(coadd['reset_flag_flagged']):
            coadd['deramped_data_reset_flagged'] = coadd['deramped_data_flagged'][coadd['reset_flag_flagged']]
        else:
            coadd['deramped_data_reset_flagged'] = coadd['deramped_data_flagged'][coadd['reset_flag_flagged'][1:]]

    return cds, coadd, base_data


def flux_slope_corr(data_file, pixel_flags=pixel_flags,
                    status_flags=status_flags, return_data=False,
                    detector=detector):

    cds, coadd, base_data = get_flagged_tods(data_file, pixel_flags,
                                             status_flags, nside=512,
                                             detector=detector)

    dirbe_data = healpy.read_map('/home/eirik/data/akari_analysis/CG_DIRBE_08_I_n0512_DR2.fits')
    coadd['pixs_flagged'] = coadd['pixs'][coadd['flags']]
    diffs = coadd['adu'][1:] - coadd['adu'][:-1]
    diffs_flagged = diffs[coadd['flags'][1:]]
    print(len(diffs_flagged))
    coadd['diffs_flagged'] = diffs_flagged
    coadd['adu_flagged'] = coadd['adu'][coadd['flags']]
    coadd['adu_reset_flagged'] = coadd['adu_flagged'][coadd['reset_flag_flagged']]
    coadd['pixs_reset_flagged'] = coadd['pixs_flagged'][coadd['reset_flag_flagged']]
    if len(diffs_flagged) < len(coadd['reset_flag_flagged']):
        coadd['diffs_reset_flagged'] = diffs_flagged[coadd['reset_flag_flagged'][1:]]
    else:
        coadd['diffs_reset_flagged'] = diffs_flagged[coadd['reset_flag_flagged']]

    flux_vals = dirbe_data[coadd['pixs_reset_flagged']]
    curr_voltages = []
    responses = []
    for i in range(1, len(coadd['adu_flagged'])):
        curr_voltages.append(coadd['adu_reset_flagged'][i-1])
        responses.append(coadd['diffs_reset_flagged'][i-1] / flux_vals[i])
#        fluxes.append(flux_vals[i])
#        slopes.append(coadd['adu'][i] - coadd['adu'][i-1])

    if return_data:
        return curr_voltages, responses, coadd, flux_vals
    else:
        return curr_voltages, responses
