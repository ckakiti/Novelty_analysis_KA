# Preprocessing
Follow these instructions if you've acquired data from the Xbox Kinect using the Datta lab's custom video capture software.
Processes raw data so it can be analysed with DLC or MoSeq

## Raw data structure
```
.
+-- Group_name
|   +-- animal1_name
|       +-- date1
|           +-- session_yyyymmddhhmmss
|               +-- depth.dat
|               +-- depth_ts.txt
|               +-- metadata.json
|               +-- rgb.mp4
|               +-- rgb_ts.txt
|       +-- date2
|       +-- date3
|   +-- animal2_name
|   +-- animal3_name
```

## Separating DLC- and MoSeq-specific files
#### 1. create blank folder on computer that can run DLC and MoSeq
 - this document assumes you're using "alex" computer and copying data first to ```/media/alex/DataDrive1/MoSeqData/```
 - name of folder will be referred to as ```groupname``` in this document
#### 2. within ```groupname``` create another blank folder labeled ```groupname_MoSeq```
 - this folder structure is important for later scripts
#### 3. copy raw data to ```groupname_MoSeq```
 - tip: if you just want to transfer certain files (e.g. only .mp4 and timestamp files), use rsync:
```
cd /location/of/raw/data
rsync -a --include '*/' --include '*.mp4' --include 'rgb_ts*' --exclude '*' . /media/alex/DataDrive1/MoSeqData/groupname/groupname_MoSeq/
```
#### 4. in MATLAB, run [MoSeqMoveRGB.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MoSeqAnalysis/MoSeqMoveRGB.m)
 - edit lines 3 and 5 to match your data
 - this will create another folder within ```groupname``` called ```groupname_DLC``` that contains all rgb.mp4 and rgb_ts.txt files (transferred out of ```groupname_MoSeq```, folder structure preserved)
#### 5. copy ```groupname_DLC``` to folder that you want to do DLC analysis in
 - default is ```/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/```
#### 6. in MATLAB, run [Convert_video.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Convert_Video.m) to convert .mp4 files to .avi 
 - for some reason if you run DLC on raw .mp4 files, it reads the wrong frame rate (not the case if you first convert to .avi)
