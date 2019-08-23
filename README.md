# Novelty_analysis

Matlab code for novelty behavior analysis

## Workflow
0. MoveFromDir.m (usually do this manually)
1. [MarkObjPos.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/README.md#1-markobjposm)
2. [Analysis.m](https://github.com/ckakiti/Novelty_analysis_KA#2-analysism)
3. [TimeStatistic](https://github.com/ckakiti/Novelty_analysis_KA#3-timestatisticm)
4. [VideoLabeling.m](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/VideoLabeling.m)

## Script details
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


# Sample commands for using DeepLabCut and MoSeq

[Running an existing network (Korleki)](https://github.com/Rxie9596/Novelty_analysis/blob/master/Docs/Using_DLC_in_UchidaLab_Korleki.md)


[Training a New Network for DeepLabCut](https://github.com/Rxie9596/Novelty_analysis/blob/master/Docs/Training_a_new_network.md)

[MoSeq Sample Commands](https://github.com/Rxie9596/Novelty_analysis/blob/master/Docs/MoSeq_Example_Command.md)
