clear
close all
clc

Config_NovAna_trim

path_to_your_file = ['/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/'...
    'Behavior/Standard setup/CombineAnalysis/Dataset_20190723'];
%    'Behavior/Standard setup/CombineAnalysis/for_statistics'];
cd(path_to_your_file)
% timeStat = readtable('TimeStatistic_combine3_tail.csv');
% timeStat = readtable('TimeStatistic_SAP_combine3_8cm_norm_incl.csv');
timeStat = readtable('TimeStatistic_combine3_tail.csv');

% timeStat_half = timeStat(1:height(timeStat)/2,:); disp('top half')
timeStat_half = timeStat(height(timeStat)/2+1:height(timeStat),:); disp('bottom half')
timeStat_half = timeStat_half{:,3:end};


path_to_MiceIndex = ['/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/' ...
    'Behavior/Standard setup/CombineAnalysis/Dataset_20190723/MiceIndex_combine3'];
%     'Behavior/Standard setup/CombineAnalysis/Dataset_20190723/'...
%     'MiceIndex_combine3.m'];
run(path_to_MiceIndex)

% name of novelty conditions
cond2name = 'cont'; %'lesions';
cond1name = 'stim'; %'controls';

% get indices of mice in each novelty condition
detectCond = cat(1, Mice.novelty);
cond2 = find(detectCond=='C');
cond1 = find(detectCond=='S');
disp(['cond1: ' num2str(cond1')])
disp(['cond2: ' num2str(cond2')])

cond1_H1 = timeStat_half(cond1, 1);
cond1_H2 = timeStat_half(cond1, 2);
cond1_H  = (cond1_H1+cond1_H2)/2;
cond1_N1 = timeStat_half(cond1, 3);
cond1_N4 = timeStat_half(cond1, 4);

cond2_H1 = timeStat_half(cond2, 1);
cond2_H2 = timeStat_half(cond2, 2);
cond2_H  = (cond2_H1+cond2_H2)/2;
cond2_N1 = timeStat_half(cond2, 3);
cond2_N4 = timeStat_half(cond2, 4);

cond2_norm = cond2_N1-cond2_H2;%./cond2_H;
cond1_norm = cond1_N1-cond1_H2;%./cond1_H;

%% t-tests
stat_summary(1).data = [cond1_H2 cond1_N1]; % paired
stat_summary(2).data = [cond2_H2 cond2_N1]; % paired
stat_summary(3).data = [cond1_N1 cond2_N1]; % unpaired
% [cond1_N1, cond1_N4];
% [cond2_N1, cond2_N4];
% test_names = {'paired', 'paired', 'unpaired'}; % parametric (ttest)
test_names = {'signrank', 'signrank', 'ranksum'}; % non-parametric
[stat_summary.name] = deal(test_names{:});

for compIter = 1:size(stat_summary,2)
    clear curr_h curr_p curr_ci curr_stats
    
    if(strcmp(stat_summary(compIter).name, 'paired'))
        % paired-sample t-test
        %   5% significance level
        %  95% confidence interval
        disp('paired t-test')
        
        [curr_h, curr_p, curr_ci, curr_stats] = ...
            ttest(stat_summary(compIter).data(:,1), ...
                  stat_summary(compIter).data(:,2));
        
    elseif(strcmp(stat_summary(compIter).name, 'unpaired'))
        % unpaired-sample t-test
        disp('unpaired t-test')
        
        [curr_h, curr_p, curr_ci, curr_stats] = ...
            ttest2(stat_summary(compIter).data(:,1), ...
                   stat_summary(compIter).data(:,2));
        
    elseif(strcmp(stat_summary(compIter).name, 'signrank'))
        % Wilcoxon signed-rank test (nonparametric version of paired ttest)
        disp('wilcoxon signed-rank test')
        
        [curr_p, curr_h, curr_stats] = ...
            signrank(stat_summary(compIter).data(:,1), ...
                     stat_summary(compIter).data(:,2));
                 
    elseif(strcmp(stat_summary(compIter).name, 'ranksum'))
        % Wilcoxon rank sum test (nonparametric version of unpaired ttest)
        disp('wilcoxon rank sum test')
        
        [curr_p, curr_h, curr_stats] = ...
            ranksum(stat_summary(compIter).data(:,1), ...
                    stat_summary(compIter).data(:,2));
               
    else
        error('test type not recognized')
    end
    
    stat_summary(compIter).raw_h     = curr_h;
    stat_summary(compIter).raw_p     = curr_p;
    stat_summary(compIter).raw_stats = curr_stats;
%         stat_summary(compIter).raw_ci    = curr_ci;

end

% multiple comparisons correction (bonferroni)
[correct_p, correct_h] = bonf_holm(cat(1,stat_summary.raw_p));
cell_p = num2cell(correct_p);
cell_h = num2cell(correct_h);
[stat_summary.correct_p] = deal(cell_p{:});
[stat_summary.correct_h] = deal(cell_h{:});

disp(correct_p)
disp(correct_h)

if(0)
    save('Dataset_20190723_statistics_timeNear')
    save('Dataset_20190723_statistics_orient')
%     save('statistics_NewHope-ROTJ_timeNearObj')
%     save('statistics_7day_preexposure_combine_timeNearObj')
end

% paired ttest -> Wilcoxon signed-rank test / paired samples Wilcoxon test
% unpaired ttest -> Wilcoxon rank sum test / unpaired two-samples Wilcoxon
%      test / Mann-Whitney test

%% power test to figure out number of lesion mice needed for future (sampsizepwr)
mean1 = mean(cond2_norm);
mean2 = mean(cond1_norm);
std1  = std(cond2_norm);
%mean1 = mean(cond2_N1); %cond2
%mean2 = mean(cond1_N1); %cond1
%std1  = std(cond2_N1);  %cond2

pwr   = 0.8;

n_needed = sampsizepwr('t2',[mean1 std1],mean2,pwr,[]);
disp(['n_needed: ' num2str(n_needed) '; pwr = ' num2str(pwr)])

% with pwr=0.9
% using std for controls: n_needed = 4
% using std for lesions: n_needed = 39

% with pwr=0.8
% using std for controls: n_needed = 4
% using std for lesions: n_needed = 30


pwr_calc = sampsizepwr('t2',[mean1 std1], mean2, [], 5);
disp(['calculated power: ' num2str(pwr_calc)])

% calculated pwr = 0.18

%% NewHope-ROTJ time around obj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
