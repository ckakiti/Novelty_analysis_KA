# Preprocessing
Follow these instructions if you've acquired data from the Xbox Kinect using the Datta lab's custom video capture software.
Processes raw data so it can be analysed with DLC or MoSeq

## Raw data structure
```
.
+-- group_name
|   +-- animal1_name
|       +-- date1
|           +-- session_yyyymmddhhmmss
|               -- depth.dat
|               -- depth_ts.txt
|               -- metadata.json
|               -- rgb.mp4
|               -- rgb_ts.txt
|       +-- date2
|       +-- date3
|   +-- animal2_name
|   +-- animal3_name
```

## Separating DLC- and MoSeq-specific files
#### 1. create blank folder on computer that can run DLC and MoSeq
 - this top folder will be referred to as ```groupname``` for the rest of instructions
 - these instructions assume you're using the computer "alex" (Linux) and copying raw data into ```/media/alex/DataDrive1/MoSeqData/```
#### 2. within ```groupname``` create another blank folder labeled ```groupname_MoSeq```
 - this folder structure is important for later scripts
#### 3. copy raw data to ```groupname_MoSeq```
 - tip: if you just want to transfer certain files (e.g. only .mp4 and timestamp files), use rsync:
```
cd /location/of/raw/data
rsync -a --include '*/' --include '*.mp4' --include 'rgb_ts*' --exclude '*' . /media/alex/DataDrive1/MoSeqData/groupname/groupname_MoSeq/
```
 - rsync for specific file paths and pruning empty directories:
```
rsync -avr --include '*/' --include '*/190413/*/proc/results_00.yaml' --include '*/190413/*/proc/results_00.h5' --exclude '*' --prune-empty-dirs . /media/alex/DataDrive1/MoSeqData/combine3L/Capoeira_MoSeq/
```

#### 4. in MATLAB, run [MoSeqMoveRGB.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MoSeqAnalysis/MoSeqMoveRGB.m)
 - edit lines 3 and 5 to match your data
 - this will create another folder within ```groupname``` called ```groupname_DLC``` that contains all rgb.mp4 and rgb_ts.txt files (transferred out of ```groupname_MoSeq```, folder structure preserved)

## DLC-specific preprocessing
#### 5. copy ```groupname_DLC``` to folder that you want to do DLC analysis in
 - default is ```/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/```
#### 6. in MATLAB, run [Convert_video.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Convert_Video.m)
 - for some reason if you run DLC on raw .mp4 files, it reads the wrong frame rate (not the case if you first convert to .avi)
 - edit line 7 to match your data
#### 7. follow instructions to either run DLC on .avi files with [existing network](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Using_DLC_in_UchidaLab_Korleki.md) or [create a new network](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Training_a_new_network.md)
#### 8. follow workflow to [analyze DLC-generated files](https://github.com/ckakiti/Novelty_analysis_KA#dlc-workflow)
## MoSeq-specific preprocessing
#### 9. create copy of [Shell_Script_Template.sh](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MoSeqAnalysis/Shell_Script_Template.sh) and move it to ```groupname```
#### 10. edit ```Shell_Script.sh``` (created in step 9)
 - this is a document to streamline extraction, modeling, and visualization of MoSeq data (without having to execute instructions line-by-line)
 - for extraction, you'll need to create a new line for each .dat file in your dataset (see instructions in template file to quickly get path for each file)
 - extraction and PCA steps depend on whether your data contains optogenetic fibers in them or not
 - PCA step depends on whether you're running a small or large dataset (try instructions for small dataset first, and if your computer runs out of memory use instructions for large dataset)
#### 11. activate moseq environment
 - run ```source activate moseq_1.2``` (for instructions to create this environment, see original MoSeq documentation)
#### 12. run ```Shell_Script.sh``` in terminal (within moseq environment)
 - will execute all commands you specified when editing the file in step 10
#### 13. follow workflow to [analyze MoSeq-generated files](https://github.com/ckakiti/Novelty_analysis_KA#moseq-workflow)
