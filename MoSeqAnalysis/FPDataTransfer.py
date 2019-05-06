import numpy as np
import scipy.io as sio

#save_directory='/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry_DLC/Francisco/temp/Francisco_190220_nidaq'
#filename = '/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry_DLC/Francisco/temp/Francisco_190220_nidaq.dat'
save_directory = '/media/alex/DataDrive1/MoSeqData/Iku_photometry2/Iku_photometry2_MoSeq/Vegas/190427/session_20190427122558/MoSeqFP.mat'
filename = '/media/alex/DataDrive1/MoSeqData/Iku_photometry2/Iku_photometry2_MoSeq/Vegas/190427/session_20190427122558/nidaq.dat'
nch=4
dtype='<f8'
with open(filename, "rb") as file_read:
    dat = np.fromfile(file_read, dtype)

nidaq_dict = {}

for i in range(nch-1):
    nidaq_dict['ch{:02d}'.format(i)] = dat[i::nch]

nidaq_dict['tstep'] = dat[nch-1::nch]

sio.savemat(save_directory, nidaq_dict)