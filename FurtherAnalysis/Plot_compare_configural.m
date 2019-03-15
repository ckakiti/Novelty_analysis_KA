clear
close all
clc

Config_NovAna
radius_cm = 10; %8; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

group_name = {'Redwall'};%, 'Quarky_conf_DLC', 'Ice&Fire_conf_DLC'};
mouse_name = {'Luke', 'Martin', 'Matthias', 'Mattimeo'};%, 'Strange', 'Up', 'Arya', 'Jon'};
%    {'Arya', 'Brienne', 'Daenerys', 'Jon'};
%    {'Bottom', 'Charm', 'Strange', 'Up'};
timeStat_all = zeros(12, 6, length(mouse_name)); %%%%%%%%%%%%%%%%%%%%%%%%%%

for group = 1:length(group_name)
    for mouse = 1:length(mouse_name)
        cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' group_name{group}])
        
        if(isfolder(mouse_name{mouse}))    
            disp(mouse_name{mouse})
            cd(mouse_name{mouse})
            timeStat = readtable('TimeStatistic.csv');
            
            if(strcmp(mouse_name{mouse}, 'Aldehyde'))
                timeStat_all(2:12,:,mouse) = timeStat{:,2:end};
            else
                timeStat_all(:,:,mouse) = timeStat{:,2:end};
            end
        end
    end
end

XTick = {'1' '2' '3' '4' '5' '6' '7' '8' 'N1' 'N2' 'N3' 'N4'};
%XTick = {'H1' 'H2' '1' '2' '3' '4'};
X     = 1:size(timeStat_all,1);

%% time spent around certain positions in arena
currMouse = 4;

disPosFig=figure(1);
set(disPosFig, 'Position', [600 600 1200 450])

hold on
plot(X, timeStat_all(:,1,currMouse)', ...
    'color', [1 0 1], 'Marker', '*', 'LineStyle', '-')
plot(X, timeStat_all(:,2,currMouse)', ...
    'color', [1 0 0], 'Marker', '*', 'LineStyle', '-')
plot(X, timeStat_all(:,3,currMouse)', ...
    'color', [0 0 1], 'Marker', '*', 'LineStyle', '-')
plot(X, timeStat_all(:,4,currMouse)', ...
    'color', [0 0.5 1], 'Marker', '*', 'LineStyle', '-')
title([mouse_name{currMouse}, '; Time spent around objs (' ...
    num2str(round(Dis_te_frame/fps/60))  'min); Rad = ' num2str(radius_cm) 'cm'])
set(gca, 'FontSize', 14)
xlabel('Session')
ylabel('Fraction of time around position')
legend('Pos 1', 'Pos 2', 'Pos 3', 'Pos 4')
xlim([0 X(end)+1])
xticks(X);
xticklabels(XTick);

%cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' group_name{1}])% '/' ...
%    mouse_name{currMouse}])
%saveas(disPosFig,'timeNearPos')
%saveas(disPosFig,['timeNearPos_' mouse_name{currMouse} '.tif'])

%% time spent around objects (not position)
currMouse = 4;
novelObjs = 1:2;

timeStat_switch = timeStat_all;
timeStat_switch(9:12,novelObjs(1),currMouse) = timeStat_all(9:12,novelObjs(2),currMouse);
timeStat_switch(9:12,novelObjs(2),currMouse) = timeStat_all(9:12,novelObjs(1),currMouse);

disObjFig=figure(1);
set(disObjFig, 'Position', [600 600 1200 450]) %[600 600 1200 450]) %[600 600 1200 950])

%subplot(3,1,1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
plot(X, timeStat_switch(:,1,currMouse)', ...
    'color', [1 0 1], 'Marker', '*', 'LineStyle', '-')
plot(X, timeStat_switch(:,2,currMouse)', ...
    'color', [1 0 0], 'Marker', '*', 'LineStyle', '-')
plot(X, timeStat_switch(:,3,currMouse)', ...
    'color', [0 0 1], 'Marker', '*', 'LineStyle', '-')
plot(X, timeStat_switch(:,4,currMouse)', ...
    'color', [0 0.5 1], 'Marker', '*', 'LineStyle', '-')
title([mouse_name{currMouse}, '; Time spent around objs (' ...
    num2str(round(Dis_te_frame/fps/60))  'min); Rad = ' num2str(radius_cm) 'cm'])
set(gca, 'FontSize', 14)
xlabel('Session')
ylabel('Fraction of time around obj')
legend('Obj 1', 'Obj 2', 'Obj 3', 'Obj 4')
xlim([0 X(end)+1])
xticks(X);
xticklabels(XTick);

if(0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,2)
plot(X, timeStat_all(:,5,currMouse)', ...
    'color', [0 0 0], 'Marker', '*', 'LineStyle', '-')
set(gca, 'FontSize', 14)
xlabel('Session')
ylabel('Total distance run (cm)')
ylim([0 inf])
xlim([0 X(end)+1])
xticks(X);
xticklabels(XTick);

subplot(3,1,3)
plot(X, timeStat_all(:,6,currMouse)', ...
    'color', [0 0 0], 'Marker', '*', 'LineStyle', '-')
set(gca, 'FontSize', 14)
xlabel('Session')
ylabel('Fraction of arena area covered')
ylim([0 1])
xlim([0 X(end)+1])
xticks(X);
xticklabels(XTick);
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' group_name{1}])% '/' ...
%    mouse_name{currMouse}])
saveas(disObjFig,['timeNearObj_' mouse_name{currMouse} '.tif'])
%saveas(disObjFig,['timeDistArea_' mouse_name{currMouse} '.tif']) %timeNearObj

%% total distance run
statColumn = 6; % 5=totalDist, 6=runArea

statAll = squeeze(timeStat_all(:,statColumn,:));

currStat=figure(1);
set(currStat, 'Position', [600 600 1200 450])

plot(repmat(X,4,1)', statAll, ...
    'color', [0 0 0], 'Marker', '*', 'LineStyle', '-')
set(gca, 'FontSize', 14)
xlabel('Session')
%ylabel('Total distance run (cm)')
ylabel('Fraction of area covered')
%ylim([0 inf])
ylim([0 1])
xlim([0 X(end)+1])
xticks(X);
xticklabels(XTick);

title([group_name{1} ': Total distance run across all mice (10 min)'], 'interpreter', 'none')
%saveas(currStat, 'runArea_all.tif')

%% summary plot #1 (scatter plot with 3 groups along x-axis)
currMouse = 1;

baseline  = 8; %5:8
novelty   = 9; %9:12
novelObjs = 1:2;
familObjs = 3:4;

avgStat_all = zeros(2,4,length(mouse_name));
avgStat_avg = zeros(4,length(mouse_name));
for mouse = 1:length(mouse_name)
    avgStat_all(:,:,mouse) = [mean(timeStat_all(baseline,1:4,mouse),1); ...
                              mean(timeStat_all(novelty,1:4,mouse),1)];
    avgStat_avg(1,mouse) = mean(avgStat_all(1,familObjs,mouse),2);
    avgStat_avg(2,mouse) = mean(avgStat_all(1,novelObjs,mouse),2);
    avgStat_avg(3,mouse) = mean(avgStat_all(2,familObjs,mouse),2);
    avgStat_avg(4,mouse) = mean(avgStat_all(2,novelObjs,mouse),2);
end

avgStatOne = figure(2);

hold on
errorbar(1, mean(avgStat_all(1,:,currMouse),2), ...
            std(avgStat_all(1,:,currMouse),0,2)/size(avgStat_all,2), ...
    'LineWidth', 2, 'Marker', 'x')
errorbar(3, mean(avgStat_all(2,novelObjs,currMouse),2), ...
            std(avgStat_all(2,novelObjs,currMouse),0,2)/size(avgStat_all,2), ...
    'LineWidth', 2, 'Marker', 'x')
errorbar(2, mean(avgStat_all(2,familObjs,currMouse),2), ...
            std(avgStat_all(2,familObjs,currMouse),0,2)/size(avgStat_all,2), ...
    'LineWidth', 2, 'Marker', 'x')
willBeNovel = scatter(repmat(1.1, length(avgStat_all(1,novelObjs,currMouse)),1)', ...
    avgStat_all(1,novelObjs,currMouse), 'k', 'filled');
willBeFamil = scatter(repmat(1.1, length(avgStat_all(1,familObjs,currMouse)),1)', ...
    avgStat_all(1,familObjs,currMouse), 'k');
scatter(repmat(3.1, length(avgStat_all(2,novelObjs,currMouse)),1)', ...
    avgStat_all(2,novelObjs,currMouse), 'k', 'filled')
scatter(repmat(2.1, length(avgStat_all(2,familObjs)),1)', ...
    avgStat_all(2,familObjs,currMouse), 'k')
xlim([0 4])
ylim([0 max(max(timeStat_all(baseline(1):novelty(end),1:4,currMouse)))])
ylabel('Time around objects (s)')
title(['Time around objects (', mouse_name{currMouse}, '; r=', num2str(radius_cm), 'cm)'], ...
    'Interpreter', 'none')
legend([willBeNovel willBeFamil], {'Will be novel', 'Will be familiar'}, 'Location', 'Northwest')
set(gca, 'FontSize', 14, 'XTick', 1:3, ...
    'XTickLabel', {['sessions ' num2str(baseline(1)) ':' num2str(baseline(end))], ...
    'familiar', 'novel'})

%cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Ice&Fire_conf_DLC/' ...
%    mouse_name{currMouse} '/181223/Analyzed_Data'])
%saveas(avgStatOne,['avgStat_' mouse_name{currMouse} '_8-9.tif'])
%disp('saved')

%% summary plot #2 (scatter/line plot, will be novel vs will be familiar   and   novel vs familiar)
cond_6OHDA = 1:2; % [2 3]; %%%%%%%%%%%%%%%%%%%%%%%%%
cond_saline = 3:4; % [1 4]; %%%%%%%%%%%%%%%%%%%%%%%%

avgStatAll = figure(3);
set(avgStatAll, 'Position', [480 450 1250 550])

hold on
errorbar(1, mean(avgStat_avg(1,:),2), ...
            std(avgStat_avg(1,:),0,2)/size(avgStat_avg,2), ...
    'LineWidth', 2, 'Marker', 'x', 'Color', 'k')
errorbar(2, mean(avgStat_avg(2,:),2), ...
            std(avgStat_avg(2,:),0,2)/size(avgStat_avg,2), ...
    'LineWidth', 2, 'Marker', 'x', 'Color', 'k')
errorbar(3, mean(avgStat_avg(3,:),2), ...
            std(avgStat_avg(3,:),0,2)/size(avgStat_avg,2), ...
    'LineWidth', 2, 'Marker', 'x', 'Color', 'k')
errorbar(4, mean(avgStat_avg(4,:),2), ...
            std(avgStat_avg(4,:),0,2)/size(avgStat_avg,2), ...
    'LineWidth', 2, 'Marker', 'x', 'Color', 'k')

basePlot6OHDA = plot([repmat(1.1, size(avgStat_avg,2)/2, 1), ...
                      repmat(2.1, size(avgStat_avg,2)/2, 1)]', ...
                     [avgStat_avg(1,cond_6OHDA); avgStat_avg(2,cond_6OHDA)], 'ko-');%'ro-'); %%%%%
                 
basePlotSaline = plot([repmat(1.1, size(avgStat_avg,2)/2, 1), ...
                       repmat(2.1, size(avgStat_avg,2)/2, 1)]', ...
                      [avgStat_avg(1,cond_saline); avgStat_avg(2,cond_saline)], 'ko-');
                  
novelPlot6OHDA = plot([repmat(3.1, size(avgStat_avg,2)/2, 1), ...
                       repmat(4.1, size(avgStat_avg,2)/2, 1)]', ...
                      [avgStat_avg(3,cond_6OHDA); avgStat_avg(4,cond_6OHDA)], 'ko-');%'ro-'); %%%%%
                  
novelPlotSaline = plot([repmat(3.1, size(avgStat_avg,2)/2, 1), ...
                        repmat(4.1, size(avgStat_avg,2)/2, 1)]', ...
                       [avgStat_avg(3,cond_saline); avgStat_avg(4,cond_saline)], 'ko-');
                   
% scatter(repmat(1.1, size(avgStat_avg,2), 1)', ...
%     avgStat_avg(1,:), 'k')%, 'filled')
% scatter(repmat(2.1, size(avgStat_avg,2), 1)', ...
%     avgStat_avg(2,:), 'k', 'filled')
% scatter(repmat(3.1, size(avgStat_avg,2), 1)', ...
%     avgStat_avg(3,:), 'k')%, 'filled')
% scatter(repmat(4.1, size(avgStat_avg,2), 1)', ...
%     avgStat_avg(4,:), 'k', 'filled')

%set(basePlot, {'color'}, {[1 0 0]; [0 0.8 0]; [0 0 1]; [0 0.8 1]})
%set(novelPlot, {'color'}, {[1 0 0]; [0 0.8 0]; [0 0 1]; [0 0.8 1]})

xlim([0 5])
ylim([0 max(avgStat_avg(:))])
ylabel('Time around objects (s)')
title(['Time around objects (n=4 mice); r=', num2str(radius_cm), 'cm'], ...
    'Interpreter', 'none')
%legend([basePlot6OHDA; basePlotSaline], {'6OHDA', '', 'saline', ''})
set(gca, 'FontSize', 14, 'XTick', 1:size(avgStat_avg,2), ...
    'XTickLabel', {'will be familiar (8)', 'will be novel (8)' 'familiar (9)', 'novel (9)'})


%cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Configural_stimuli_DLC/')
%saveas(avgStatAll,'avgStat_allMice_8-9.tif')
%disp('saved')

%% summary plot #3 (scatter/line plot, will be familiar vs familiar   and   will be novel vs novel)
cond_6OHDA = [2 3]; %%%%%%%%%%%%%%%%%%%%%%%%%%
cond_saline = [1 4]; % 1:8; %%%%%%%%%%%%%%%%%%%%%%%%

avgStatAll2 = figure(4);
set(avgStatAll2, 'Position', [480 450 1250 550])

hold on
% errorbar(1, mean(avgStat_avg(1,:),2), ...
%             std(avgStat_avg(1,:),0,2)/size(avgStat_avg,2), ...
%     'LineWidth', 2, 'Marker', 'x')
% errorbar(2, mean(avgStat_avg(2,:),2), ...
%             std(avgStat_avg(2,:),0,2)/size(avgStat_avg,2), ...
%     'LineWidth', 2, 'Marker', 'x')
% 
% errorbar(3, mean(avgStat_avg(3,:),2), ...
%             std(avgStat_avg(3,:),0,2)/size(avgStat_avg,2), ...
%     'LineWidth', 2, 'Marker', 'x')
% errorbar(4, mean(avgStat_avg(4,:),2), ...
%             std(avgStat_avg(4,:),0,2)/size(avgStat_avg,2), ...
%     'LineWidth', 2, 'Marker', 'x')

fam6OHDA = plot([repmat(1.1, size(avgStat_avg,2)/2, 1), repmat(2.1, size(avgStat_avg,2)/2, 1)]', ...
    [avgStat_avg(1,cond_6OHDA); avgStat_avg(3,cond_6OHDA)], 'ko-');%'ro-');%, 'filled') %%%%%

famSaline = plot([repmat(1.1, size(avgStat_avg,2)/2, 1), repmat(2.1, size(avgStat_avg,2)/2, 1)]', ...
    [avgStat_avg(1,cond_saline); avgStat_avg(3,cond_saline)], 'ko-');%, 'filled')

nov6OHDA = plot([repmat(3.1, size(avgStat_avg,2)/2, 1), repmat(4.1, size(avgStat_avg,2)/2, 1)]', ...
    [avgStat_avg(2,cond_6OHDA); avgStat_avg(4,cond_6OHDA)], 'ko-', 'MarkerFaceColor', 'k'); % 'ro-'%%%

novSaline = plot([repmat(3.1, size(avgStat_avg,2)/2, 1), repmat(4.1, size(avgStat_avg,2)/2, 1)]', ...
    [avgStat_avg(2,cond_saline); avgStat_avg(4,cond_saline)], 'ko-', 'MarkerFaceColor', 'k');

xlim([0 5])
ylim([0 max(avgStat_avg(:))])
ylabel('Fractional time around objects')
title(['Fractional time around objects (n=' num2str(size(avgStat_avg,2)) ...
    ' mice); r=', num2str(radius_cm), 'cm'], 'Interpreter', 'none')
%legend([fam6OHDA; famSaline], {'6OHDA', '', 'saline', ''})
set(gca, 'FontSize', 14, 'XTick', 1:4, ...
    'XTickLabel', {'will be familiar (8)', 'familiar (9)' 'will be novel (8)', 'novel (9)'})

%cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Configural_stimuli_DLC/')
%saveas(avgStatAll2,'avgStat_allMice_2_8-9.tif')
%disp('saved')

% paired ttest results: (compare all mice, session 8 vs 9)
%  h_familiar = 0; p_familiar = 0.3486 
%  h_novel    = 0; p_novel    = 0.5739 
