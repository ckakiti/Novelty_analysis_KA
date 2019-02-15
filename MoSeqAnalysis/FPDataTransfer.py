import numpy as np
import scipy.io as sio

save_directory='/media/alex/DataDrive1/MoSeqData/SynTest/session_20190213104251/MoSeqFP_Miami_190125.mat'
filename = '/media/alex/DataDrive1/MoSeqData/SynTest/session_20190213104251/nidaq.dat'
nch=4
dtype='<f8'
with open(filename, "rb") as file_read:
    dat = np.fromfile(file_read, dtype)

nidaq_dict = {}

for i in range(nch-1):
    nidaq_dict['ch{:02d}'.format(i)] = dat[i::nch]

nidaq_dict['tstep'] = dat[nch-1::nch]

sio.savemat(save_directory, nidaq_dict)