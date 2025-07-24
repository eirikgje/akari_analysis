from akari_data_processing import flux_slope_corr
from akari_defaults import num_processes
import numpy as np
import glob


data_files = glob.glob(f'/home/eirik/data/akari_analysis/akari_TSD/2006_04/FIS_SW_200604132*')
print(data_files)

nbins = 100
detector = 95


def mask_and_bin(voltages, responses, bins):
    mask = (responses < 0) | (np.isinf(responses))
    mask = ~mask
    masked_voltages = voltages[mask]
    masked_responses = responses[mask]
    bin_assign = np.searchsorted(bins, masked_voltages)
    print(bin_assign)
    nbins = len(bins) - 1
    binsum = np.bincount(bin_assign, weights=masked_responses, minlength=nbins)
    print(bin_assign)
    hits = np.bincount(bin_assign, minlength=nbins)
    binmean = np.zeros_like(binsum)
    binmean[np.where(hits != 0)] = binsum[np.where(hits !=0)] / hits[np.where(hits != 0)]
    mean_subtracted = masked_responses - binmean[bin_assign]
    binvar = np.bincount(bin_assign, weights=mean_subtracted**2, minlength=nbins)
    binvar[np.where(hits != 0)] = binvar[np.where(hits !=0)] / hits[np.where(hits != 0)]
    binstd = np.sqrt(binvar)
    return binmean, binstd


def main(data_files=data_files, nbins=nbins, detector=detector):

    #    data_file = (f'/home/eirik/data/akari_analysis/akari_TSD/2006_04/FIS_SW_20060422190000_gb.hdf5')
    tot_voltages = []
    tot_responses = []
    tot_means = []
    tot_stds = []

    for data_file in data_files:
        voltages, responses, coadd, flux_vals = flux_slope_corr(
                data_file, return_data=True, detector=detector)
        voltages = np.array(voltages)
        responses = np.array(responses)
        tot_voltages.append(voltages)
        tot_responses.append(responses)
    allvoltages = []
    allresponses = []
    for i in range(len(tot_voltages)):
        allvoltages.extend(tot_voltages[i])
        allresponses.extend(tot_responses[i])
    allvoltages = np.array(allvoltages)
    allresponses = np.array(allresponses)
    binmin = np.min(allvoltages)
    binmax = np.max(allvoltages)
    bins = np.linspace(binmin, binmax, nbins+1)
    meanofall, stdofall = mask_and_bin(allvoltages, allresponses, bins)
    for voltages, responses in zip(tot_voltages, tot_responses):
        currmean, currstd = mask_and_bin(voltages, responses, bins)
        tot_means.append(currmean)
        tot_stds.append(currstd)
    return meanofall, stdofall, tot_means, tot_stds, tot_voltages, tot_responses, bins
#        tot_voltages.extend(list(voltages))
#        tot_responses.extend(list(responses))
