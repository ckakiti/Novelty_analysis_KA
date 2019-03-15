% photometry analysis

clear
clc
close all

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry_DLC/Francisco

file_curr = 'Francisco_190220_rgbDeepCut_resnet50_MoSeqNovelty_RetrainSep17shuffle1_1030000.csv';

Labels = csvread(file_curr, 3, 0);