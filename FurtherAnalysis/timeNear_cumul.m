%% plot cumulative histograms
clear
close all
clc

Config_NovAna

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/StandardSetup_combine/

run('MiceIndex_combine4')
detectCond = cat(1, Mice.novelty);
cond2 = find(detectCond=='C');
cond1 = find(detectCond=='S');
reorder = [cond2; cond1];

files    = dir;
whichDir = [files.isdir];
nameDir  = files(whichDir);
nameDir  = {nameDir.name};
nameDir(ismember(nameDir,{'.','..','temp'})) = [];

timeNearN1 = [];
global_max = 0;

durTotal = 10; % duration of analysis (min)
curr_day = 3;  % N1
disp(['Duration of analysis: ' num2str(durTotal) 'min'])
disp(['analyzing day: ' num2str(curr_day)])

Dis_ts_frame=500;
Dis_te_frame=durTotal.*60.*fps+Dis_ts_frame;
startframe  = Dis_ts_frame;
endframe    = Dis_te_frame;

%%
for fileiter = 1:length(reorder) %1:length(nameDir)
    cd(nameDir{reorder(fileiter)})
    disp(nameDir{reorder(fileiter)})
    cd Analyzed_Data_1obj_12cm_tail
    
    conv = dir('*Converted.mat');
    load(conv(3).name)
%     isNear = Labels(1:26000,21);
%     timeNearN1 = [timeNearN1 isNear];
    endframe = min(endframe, size(Labels,1));
    [counts, centers] = hist(Labels(startframe:endframe,17),0:70);
    count_cumul = cumsum(counts);
    timeNearN1 = [timeNearN1; count_cumul];
    
    if(max(Labels(:,17) > global_max))
        global_max = max(Labels(:,17));
    end
    
    cd ../..
end

timeNearN1_cumul = timeNearN1/15;
% timeNearN1_cumul = cumsum(timeNearN1)/15;
% xaxis_s          = (1:26000)/15;
% timeNearN1_x     = repmat(xaxis_s',1,length(nameDir));

disp('end')

%%
load('SAP_plus_dist_all.mat')

bin_size = 100;
cumul_sap_dist = [];

for miceiter = 1:length(dist_of_sap_all)
    for dayiter = curr_day %1:length(dist_of_sap_all(miceiter).sapDist)
        curr_saps = dist_of_sap_all(miceiter).sapDist{dayiter};
        
        [counts, centers] = hist(curr_saps,0:bin_size);
        count_cumul = cumsum(counts);
        cumul_sap_dist = [cumul_sap_dist; count_cumul]; 
        
        if(max(curr_saps) > global_max)
            global_max = max(curr_saps);
        end
    end
end

close all
clc

fig1 = figure(1);
hold on
cont = plot(centers, cumul_sap_dist(cond2,:), 'r-');
stim = plot(centers, cumul_sap_dist(cond1,:), 'k-');
line([radius_cm radius_cm], [0 max(cumul_sap_dist(:))], 'LineStyle', '--')

legend([cont(1) stim(1)], {'cont', 'stim'}, 'location','east') %nameDir(reorder)
title(['Cumulative SAP distance from object (' num2str(durTotal) 'min): N1'])
xlim([0 global_max])
xlabel('Distance from obj (cm)')
ylabel('Cumulative time at distance from obj (s)')
set(gca, 'FontSize', 16)
set(fig1, 'Position', [200 200 700 600])

disp('plotted')

if(0)
    saveas(gca, ['SAPnearN1_' num2str(durTotal) 'min_cumulative.tif'])
end

%%
close all
clc

fig1 = figure(1);
hold on
l = plot(centers, timeNearN1_cumul(1:length(cond2), :), 'r-');
c = plot(centers, timeNearN1_cumul(length(cond2)+1:end,:), 'k-');
line([radius_cm radius_cm], [0 max(timeNearN1_cumul(:))], 'LineStyle', '--')

legend([l(1) c(1)], {'6OHDA', 'saline'}, 'location','east') %nameDir(reorder)
title(['Cumulative distance from object (' num2str(durTotal) 'min): N1'])
xlim([0 global_max])
xlabel('Distance from obj (cm)')
ylabel('Cumulative time at distance from obj (s)')
set(gca, 'FontSize', 16)
set(fig1, 'Position', [200 200 700 600])

disp('plotted')

if(0)
    saveas(gca, ['timeNearN1_' num2str(durTotal) 'min_cumulative.tif'])
end