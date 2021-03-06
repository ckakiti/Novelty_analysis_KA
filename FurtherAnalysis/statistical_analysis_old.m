clear
close all
clc

Config_NovAna

%cd('/media/alex/TOSHIBA EXT/DLC_Analyzed_Videos/NewHope-ROTJ/')
cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/7day_preexposure_combine')
timeStat = readtable('TimeStatistic.csv');
timeStat2 = timeStat{:,3:end};

cond2name = 'cont';
cond1name = 'stim';
cond2 = [1:6] + 12; % add 12 for orientation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cond1 = [7:12] + 12;

cond2_H2 = timeStat2(cond2, 2);
cond2_N1 = timeStat2(cond2, 3);
cond2_N2 = timeStat2(cond2, 4);
cond1_H2 = timeStat2(cond1, 2);
cond1_N1 = timeStat2(cond1, 3);
cond1_N2 = timeStat2(cond1, 4);
%
%
% paired-sample t-test
%   5% significance level
%  95% confidence interval 
[cond1_H2N1_h, cond1_H2N1_p, cond1_H2N1_ci, cond1_H2N1_stats] = ttest(cond1_H2, cond1_N1);
[cond1_N1N2_h, cond1_N1N2_p, cond1_N1N2_ci, cond1_N1N2_stats] = ttest(cond1_N1, cond1_N2);

[cond2_H2N1_h, cond2_H2N1_p, cond2_H2N1_ci, cond2_H2N1_stats] = ttest(cond2_H2, cond2_N1);
[cond2_N1N2_h, cond2_N1N2_p, cond2_N1N2_ci, cond2_N1N2_stats] = ttest(cond2_N1, cond2_N2);

[cond1_cond2_N1_h, cond1_cond2_N1_p, cond1_cond2_N1_ci, cond1_cond2_N1_stats] ...
    = ttest(cond1_N1, cond2_N1);

% save('statistics_NewHope-ROTJ_timeNearObj')
% save('statistics_7day_preexposure_combine_timeNearObj')

% NewHope-ROTJ time around obj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stim stats
%   mean_H2 0.0293
%   mean_N1 0.0412
%   h       0
%   p       0.2677
%   ci      [-0.0362 0.0126]
%   stats   tstat -1.2468
%           df     5
%           sd     0.0232

% cont stats
%   mean_H2 0.0220
%   mean_N1 0.1462
%   h       1
%   p       0.0024
%   ci      [-0.1803 -0.0680]
%   stats   tstat -5.6822
%           df     5
%           sd     0.0535

% effect size: 0.1050
%  mean(cond2_N1)-mean(cond1_N1)

% SEM (cont): 0.0218
%  cond2_stats.sd/sqrt(length(cond2))

% NewHope-ROTJ orientation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stim stats
%   mean_H2 0.0483
%   mean_N1 0.1150
%   h       1
%   p       3.7547e-05
%   ci      [-0.0792 -0.0541]
%   stats   tstat -13.6719
%           df     5
%           sd     0.0119

% cont stats
%   mean_H2 0.0422
%   mean_N1 0.1064
%   h       1
%   p       8.1179e-04
%   ci      [-0.0872 -0.0413]
%   stats   tstat -7.1870
%           df     5
%           sd     0.0219


% 7-day preexposure time around obj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stim stats
%   mean_H2 0.0272
%   mean_N1 0.0366
%   h       0
%   p       0.1101
%   ci      [-0.0218 0.0030]
%   stats   tstat -1.9400
%           df     5
%           sd     0.0118

% cont stats
%   mean_H2 0.0259
%   mean_N1 0.1020
%   h       1
%   p       0.0456
%   ci      [-0.1501 -0.0022]
%   stats   tstat -2.6469
%           df     5
%           sd     0.0705

% 7-day preexposure orientation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stim stats
%   mean_H2 0.0451
%   mean_N1 0.1027
%   h       1
%   p       0.0107
%   ci      [-0.0675 -0.0144]
%   stats   tstat -3.9655
%           df     5
%           sd     0.0253

% cont stats
%   mean_H2 0.0452
%   mean_N1 0.0862
%   h       1
%   p       0.0044
%   ci      [-0.0878 -0.0275]
%   stats   tstat -4.9119
%           df     5
%           sd     0.0287
