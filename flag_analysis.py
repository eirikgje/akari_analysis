import h5py
import numpy as np
import healpy
import glob

detector = 95
data_files = glob.glob('/home/eirik/data/akari_analysis/akari_TSD/2006_04/FIS_SW_*.hdf5')


def fraction_analysis(data_files=data_files, detector=detector):
    sums = np.zeros(32)
    numsamps = 0

    for data_file in data_files:
        hdf_file = h5py.File(data_file)
        pixflags = np.array(hdf_file['FIS_OBS']['PIX_FLAG'])[:, detector-1, :]
        numsamps += len(pixflags)
        for i in range(32):
            sums[i] += np.sum(pixflags[:, i])

    for i in range(32):
        print(f'Detector {i+1}, Fraction of excluded data: {sums[i]/numsamps}')


def coverage_analysis(data_file=data_files[200], detector=detector):
    coverage_mat = np.zeros((32, 32), dtype=bool)
    hdf_file = h5py.File(data_file)
    pixflags = np.array(hdf_file['FIS_OBS']['PIX_FLAG'])[:, detector-1, :]
    for i in range(32):
        for j in range(i, 32):
            coverage_mat[i, j] = not np.any(pixflags[:, j][~pixflags[:, i]])
            if j != i:
                coverage_mat[j, i] = not np.any(pixflags[:, i][~pixflags[:, j]])
    return coverage_mat


def is_bad_pixel_always_set_of_others(data_file=data_files[200], detector=detector):
    hdf_file = h5py.File(data_file)
    pixflags = np.array(hdf_file['FIS_OBS']['PIX_FLAG'])[:, detector-1, :]
#    question_answer = pixflags[:, 0]
    mask = pixflags[:, 0]
    any_flags = np.array([False]*int(np.sum(pixflags[:, 0])))
    for i in range(1, 32):
        if i == 26:
            continue
        any_flags = any_flags | pixflags[mask, i]
    return any_flags

if __name__ == '__main__':
#    fraction_analysis()
#    coverage_mat = coverage_analysis()
    any_flags = is_bad_pixel_always_set_of_others()

