# Instructions for using DeepLabCut and MoSeq

[Preprocessing raw video data](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Preprocessing.md)

[Training a New Network for DeepLabCut](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Training_a_new_network.md)

[Running an existing network (Korleki)](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Using_DLC_in_UchidaLab_Korleki.md)

[MoSeq Commands](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/MoSeq_Example_Command.md)

[Running MoSeq on the cluster](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/MoSeq_on_cluster.md)

# Novelty analysis code

Code for novelty behavior analysis of processed data (see [Preprocessing.md](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Preprocessing.md) for instructions on processing raw video data)

Following sections:
- [DLC workflow](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#dlc-workflow)
- [MoSeq workflow](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#moseq-workflow)
- [Fiber photometry workflow](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#fiber-photometry-workflow)

## DLC Workflow
0. [Config_NovAna.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#0-config_novanam)
1. [MarkObjPos.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#1-markobjposm)
2. [Analysis.m](https://github.com/ckakiti/Novelty_analysis_KA#2-analysism)
3. [TimeStatistic.m](https://github.com/ckakiti/Novelty_analysis_KA#3-timestatisticm)
4. [MiceIndex.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#4-miceindexm)
5. [Plot_compare.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#5-plot_comparem)
6. [bout_analysis.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#4-bout_analysism)
7. analy_novelty_multi_sessions_2012

## DLC Script Details
### 0. [Config_NovAna.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Config_NovAna.m)
- basic information for current dataset including name of network used to run DeepLabCut, video file format (.mp4 or .avi), vid length and width, fps, radius_cm, and angle radius
- IMPORTANT: info is different for each dataset, double check before running other analysis scripts
- Input: manual editing of file
- Output: parameters saved in .m file to be loaded in later scripts

### 1. [MarkObjPos.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MarkObjPos.m)
- Manually select arena boundaries and object position for all sessions
- Has option to select LED position and extrapolate object positions for other 3 corners of arena (4obj)
- Input: Config_NovAna, .csv output files from DLC
- Output: Arena_Obj_Pos.m and/or Arena_Obj_Pos_4obj.m

### 2. [Analysis.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Analysis.m)
- Extracting DLC-based information - amount of time spent near object, angle relative to object, velocity, etc.
- Input: Config_NovAna, Arena_Obj_Pos.mat
- Output: Trajectories and Heatmap from first 10 min of each session, .mat file containing Labels and radius_cm

### 3. [TimeStatistic.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/FurtherAnalysis/TimeStatistic.m)
- create summary array of time spent near object, time spent in periphery of arena, or total distance run
- similar to summary array created by bout_analysis.m and area_analysis.m
- Input: Config_NovAna, .mat files from Analysis.m
- Output: TimeStatistic.csv, TimeStatistic_body_periph, TimeStatistic_nose_totalDistCut

### 4. MiceIndex.m
- contains information about mouse name and novelty condition for each mouse (optional: MSid and date for MoSeq analysis, manual transfer from moseq2-index.yaml file generated by MoSeq analysis)
- unique for each dataset
- Input: manual editing of file
- Output: parameters ('Mice' structure) saved in m. file to be loaded in later scripts

### 5. [Plot_compare.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/FurtherAnalysis/Plot_compare.m)
- plotting different scalar statistics across sessions and mice (e.g. time spent around object, orientation, bout num, bout len, total  distance run, area covered)
- Input: Config_NovAna, MiceIndex, summary file (e.g. TimeStatistic, TimeStatistic_nose_totalDistCut, boutAnalysis_nose)
- Output: tif files (timeNearObj_10min.tif, orientToObj_10min.tif, boutNum_10min_nose.tif, boutLen_10min_nose.tif, totalDistanceRun_10min_nose.tif, etc)

### 6. [bout_analysis.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/FurtherAnalysis/bout_analysis.m)
- similar to TimeStatistic, create summary array of number of bouts and average bout length across sessions and days
- can also create structure containing frame numbers for each poke and approach (needed for MoSeqEventAlignedAnalysis.m)
- Input: Config_NovAna, MiceIndex, .mat files from Analysis.m
- Output: boutAnalysis_nose.csv, PokesApproaches.mat, DatasetName_poke_labels_N1.csv

### 7. video_novelty_multi_labels_2011.m (Mitsuko code)
- creates DLC_label.mat (containing Labels and session_start)
- Input: akiti_miceID_210318.xlsx (animal name, group, condition), ArenaObjPos.m, DLC csv file
- Output: DLC_label.mat

### 8. analy_novelty2103 (Mitsuko code)
- creates bout_multi.mat (needed for plotting multi-day behavior imagesc)
- Input: akiti_miceID_210318.xlsx, ArenaObjPos.mat, DLC_label.mat
- Output: GROUPNAME_analy_novelty2103_ka.mat, bout_multi.mat, imagesc behavior plots (approach w tail behind/exposure, approach freq, etc)

### 9. analy_novelty_multi_sessions_2012 (Mitsuko code)
- creates distToObj, nose-tail, and bout plots for multiple mice and sessions
- Input: Config_NovAna, MiceIndex.mat, DLC_label.mat (generated by video_novelty_multi_labels_2011.m)
- Output: 

## MoSeq Workflow
0. [Shell_Script.sh](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#0-shell_scriptsh)
1. [ModelDataTransfer.py](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#1-modeldatatransferpy)
2. [MiceIndex_blank.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#2-miceindex_blankm)
3. [extract_uuid](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#3-extract_uuid)
4. [MoSeqEventAlignedAnalysis.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#4-moseqeventalignedanalysism)
5. [MoSeqEventBasedAnalysis.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#5-moseqeventbasedanalysism)

## MoSeq Script Details
### 0. [Shell_Script.sh](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MoSeqAnalysis/Shell_Script_Template.sh)
- contains command lines to run the following: moseq2-extract, moseq2-pca, moseq2-model, and moseq2-viz

### 1. [ModelDataTransfer.py](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MoSeqAnalysis/ModelDataTransfer.py)
- Transfers the data generated by MoSeq to MATLAB readable format
- Input: moseq2-index.yaml, my_model.p
- Output: MoSeqDataFrame.mat

### 2. MiceIndex_blank.m
- by hand, edit [MiceIndex template](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MoSeqAnalysis/Mice_Index_Template.m) to include list of mice names and which novelty condition they were assigned

### 3. [extract_uuid](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MoSeqAnalysis/extract_uuid.m)
- automated way to extract all uuids from moseq2-index.yaml file
- follow directions in header before running code (need to use visual studio)
- Inputs: MiceIndex_blank.m, moseq2-index_extract
- Output: MiceIndex.mat

### 4. [MoSeqEventAlignedAnalysis.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MoSeqAnalysis/MoSeqEventAlignedAnalysis.m)
- plot syllable usage at/around time of poke (automatically identified pokes)
- Input: MoSeqDataFrame.mat, MiceIndex.mat, DATASETNAME_poke_labels_N1.csv (generated by [bout_analysis.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/FurtherAnalysis/bout_analysis.m))
- Output: plot of syllable usage aligned to poke

### 5. [MoSeqEventBasedAnalysis.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MoSeqAnalysis/MoSeqEventBasedAnalysis.m)
- bar plot of total syllable usage (s) for each mouse, divided by novelty condition
- plot syllable usage across time bins within 1 day, divide by novelty condition
- plot average syllable usage per novelty condition across novelty days
- Inputs: MoSeqDataFrame.mat, MiceIndex.mat, a syllable of interest to plot
- Outputs: bar plot of total syllable usage, time course of usage within day, time course of usage across days

### 6. [MoSeqGeneralSyllableAnalysis.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MoSeqAnalysis/MoSeqGeneralSyllableAnalysis.m)
- plot fractional syllable usage (choice to select subset of mice and/or days, first 10 min or all frames)
- calculates how many syllables explain 90% of behavior (needed for statistical analysis)
- Inputs: MiceIndex.mat, MoSeqDataFrame.mat
- Outputs: plot of fractional syllable usage (sorted by usage, all days), plot of cumulative fractional usage

## Fiber Photometry Workflow
0. Preprocessing requirements
Mitsuko's code:
1. photometry_water_tone_Korleki_2105.m
2.
3.
4. photometry_novelty_Korleki_2105.m
Korleki's code:
a. [KW_1_LED.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#1-kw_1_ledm)
b. [FPdat_import.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#2-fpdat_importm)
c. [bout_analysis.m / by hand](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#3-bout_analysism)
d. [MoSeqFP_raster.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#4-moseqfp_rasterm)

### 0. Preprocessing requirements
- for each animal, need:
  - Arena_Obj_Pos.mat, including marking LED location (generated by [MarkObjPos.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MarkObjPos.m))
- for each photometry session, need:
  - FP_data and tone_data (generated by custom LabView code on day of recording)
  - DLC csv file

### a. [KW_1_LED.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/FurtherAnalysis/KW_1_LED.m)
- extract frames when LED was on (in rgb video)
- NOTE: TAKES A LONG TIME TO RUN
- Input:
  - video file (.mp4/.avi)
  - ArenaObjPos.mat (from MarkObjPos.m, containing LED position)
- Output: 
  -  save_###/save_###_max (pixel intensity for each frame of video - takes a while to run)
  -  LED_on (onset of LED flash across video)

### b. [FPdat_import.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/FPdat_import.m)
- read and save raw photometry signal (generated by custom program)
- Input:
  - LED_on (from KW_1_LED.m)
  - FP analog data file (generated by custom program, either .dat or no extension)
- Output:
  - mouse_date_FP.mat (containing GCaMP, tdTom, TTL, TTL_on, pos_within_TTL)
  - mouse_date_rgb_ts (containing rgb indices and corresponding FP indices)

### c. [bout_analysis.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/FurtherAnalysis/bout_analysis.m) OR hand-annotate
- annotate rgb/depth video with significant events (e.g. pokes, bout start, reward delivery)
- Input: Config_NovAna, MiceIndex, .mat files from Analysis.m
- Output: PokesApproaches.mat, DatasetName_poke_labels_N1.csv

### d. [MoSeqFP_raster.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/MoSeqAnalysis/MoSeqFP_raster.m)
- raster plot of photometry signal aligned to significant events
- Input:
  - Config_NovAna
  - depthts.txt, rgbts.txt (generated by computer during recording)
  - mouse_date_FP.mat (generated by FPdat_import.m)
  - csv of significant events (generated by bout_analysis.m or by hand; convert into FP frames within code)
- Output:
  - FP_mouse_session.tif
  - FP_raster_GCaMP_mouse_session.tif
  - FP_raster_tdTom_mouse_session.tif
