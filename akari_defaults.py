nside = 256
num_processes = 8

instrument = 'SW'
detector = 95
#data_files = [f'/home/eirik/data/akari_analysis/akari_TSD/2006_04/FIS_{instrument}_20060415090000_gb.hdf5']
#data_files = glob.glob(f'/home/eirik/data/akari_analysis/akari_TSD/2006_04/FIS_{instrument}_20060415*hdf5')
#print(len(data_files))
npix = 12 * nside ** 2

reset_flag = 6

#pixel_flags = [0]
#pixel_flags = [0, 5] # bad pixel, reset_pixflag
#pixel_flags = [0, 5, 6] # bad_pixel, reset_pixflag, reset anomaly
pixel_flags = [el for el in range(1, 32) if el not in (0, 6)]
#pixel_flags = [el for el in range(1, 32)]# if el not in (0, 6)]
status_flags = [16, 17, 18] # three lamps on
