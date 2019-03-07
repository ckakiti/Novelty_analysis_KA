import numpy as np
import scipy.io as sio

save_directory='/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry_DLC/Miami/temp/Miami_190125_nidaq'
filename = '/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry_DLC/Miami/temp/Miami_190125_nidaq.dat'
nch=3
dtype='<f8'
with open(filename, "rb") as file_read:
    dat = np.fromfile(file_read, dtype)

nidaq_dict = {}

for i in range(nch-1):
    nidaq_dict['ch{:02d}'.format(i)] = dat[i::nch]

nidaq_dict['tstep'] = dat[nch-1::nch]

sio.savemat(save_directory, nidaq_dict)