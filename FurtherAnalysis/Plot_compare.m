clear
close all
clc

Config_NovAna_Ghana

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/StandardLesion_combine/
run('MiceIndex_combineL')

oneStat = 0;

if(oneStat) 
    % only one statistic to plot (e.g. time in periphery or totalDistRun
%     timeStat = readtable('TimeStatistic_periph.csv');
%      timeStat = readtable('TimeStatistic_nose_totalDistCut.csv');
    timeStat = readtable('areaAnalysis_nose_center10.csv');
%     timeStat  = readtable('TimeStatistic_norm.csv');
    timeStat2 = timeStat{:,3:end};

%     names = timeStat{1:height(timeStat),1}';
    names = cat(1, Mice.name);
%     Y_dis = timeStat2(1:height(timeStat), :);
    Y_dis = timeStat2(1:size(timeStat2,1), :);
    Y_ang = Y_dis;
else
    % two different statistics to plot (e.g. timeNearObj/orient or boutNum/boutLen
%     timeStat = readtable('TimeStatistic_12cm_tail.csv');
%     timeStat = readtable('boutAnalysis_nose.csv');   
%     timeStat = readtable('areaAnalysis_nose_quad3-4.csv');
    timeStat = readtable('TimeStatistic_12cm_tail.csv');
    
    timeStat2 = timeStat{:,3:end};
    names = timeStat{1:height(timeStat)/2,1}';
%     names = cat(1, Mice.name);
    Y_dis = timeStat2(1:height(timeStat)/2, :);
    Y_ang = timeStat2(height(timeStat)/2+1:height(timeStat),:);
end

XTick = {'H1' 'H2' 'N1' 'N2' 'N3' 'N4' 'N5' 'N6' 'N7' 'N8' 'N9' 'N10'};
% XTick = {'0-2','2-4','4-6','6-8','8-10'};
% XTick = {'0-6','7-12','13-18','19-24','25-30'};


cond2name = '6OHDA';%'6OHDA';%'cont';%'lego@5min';%'GroupB';%'cont';%'6OHDA';
cond1name = 'saline';%'saline';%'stim';%'lego@0min';%'GroupA';%'stim';%'saline';

detectCond = cat(1, Mice.novelty);
% cond2 = find(detectCond=='C');
% cond1 = find(detectCond=='S');
cond2 = find(detectCond=='l');
cond1 = find(detectCond=='s');
disp(['cond1: ' num2str(cond1')])
disp(['cond2: ' num2str(cond2')])
cond2Color = [0.5 0.0 0.5];%[0.5 0.0 0.5]; %[0.5 0.0 0.5]; [0.5 0.0 0.5]};
cond1Color = [1.0 0.5 0.0]; %[1.0 0.5 0.0]; [1.0 0.5 0.0]};
% cond2Color = {[0.3 0.0 0.3]; [0.5 0.0 0.5]; [0.8 0.0 0.5]};
% cond1Color = {[1.0 0.1 0.0]; [1.0 0.5 0.0]; [1.0 0.7 0.0]};

XTick = XTick(1:length(Y_dis(1,:)));
X     = 1:length(Y_dis(1,:));

condavg=mean(Y_dis(cond2,:));
studavg=mean(Y_dis(cond1,:));
condstd=std(Y_dis(cond2,:));
studstd=std(Y_dis(cond1,:));

conaavg=mean(Y_ang(cond2,:));
stuaavg=mean(Y_ang(cond1,:));
conastd=std(Y_ang(cond2,:));
stuastd=std(Y_ang(cond1,:));


fps=15; %%%%%%%%%%%%%%%%%%%%%%%
totalTime = 10;
Dis_te_frame = fps*60*totalTime+500;
radius_cm = 12; %12 %8
disp(['Total time: ' num2str(round(Dis_te_frame/fps/60)) 'min'])
disp(['Radius: ' num2str(radius_cm) 'cm'])

%% for plotting time around object and orientation to object
Y_dis_title = ['Time spent at obj (', num2str(round(Dis_te_frame/fps/60))  ...
    'min); Rad = ' num2str(radius_cm) 'cm'];
Y_ang_title = ['Orientation to obj (', num2str(round(Dis_te_frame/fps/60)) ...
    'min); Deg = +-' num2str(angle_radius) char(176)];
Y_dis_ylabel = 'Time near object (frac)';
Y_ang_ylabel = 'Orientation to object (frac)';
x_label = 'Time within session (min)';%'Training day';%

close all 

% Time near obj %%%
disfig=figure(1);
set(disfig, 'Position', [600 600 1200 450])

subplot(1, 2, 2)
hold on
% d1 = plot(repmat(X, length(cond2), 1)', Y_dis(cond2,:)', 'Marker', '*', 'LineWidth',2);
% d2 = plot(repmat(X, length(cond1), 1)', Y_dis(cond1,:)', 'Marker', '*', 'LineWidth',2);
d1 = plot(repmat(X, length(cond2), 1)', Y_dis(cond2,:)', 'Color', cond2Color, 'Marker', '*', 'LineWidth',2);
d2 = plot(repmat(X, length(cond1), 1)', Y_dis(cond1,:)', 'Color', cond1Color, 'Marker', '*', 'LineWidth',2);
% set(d1, {'color'}, cond2Color)
% set(d2, {'color'}, cond1Color)
title(Y_dis_title)
set(gca, 'FontSize', 14)
% set(gca, 'FontSize', 18)
xlabel(x_label)
ylabel(Y_dis_ylabel)
% legend(names(cond1))
% legend([names(cond2) names(cond1)])
legend([d1(1) d2(1)], [cond2name ' (n=' num2str(length(cond2)) ')'], ...
  [cond1name ' (n=' num2str(length(cond1)) ')']) %{'cont', 'stim'})
xlim([0 X(end)+1])
% ylim([0 0.6])
% ylim([0 1])
ylimDis = disfig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);
% set(gca,'YTick',[0 0.2 0.4])

subplot(1, 2, 1)
hold on
errorbar(X,condavg,condstd, 'Color', cond2Color, 'LineWidth',2)
errorbar(X,studavg,studstd, 'Color', cond1Color, 'LineWidth',2)
title(Y_dis_title)
set(gca, 'FontSize', 14)
% set(gca, 'FontSize', 18)
xlabel(x_label)
ylabel(Y_dis_ylabel)
legend([cond2name ' (n=' num2str(length(cond2)) ')'], ...
   [cond1name ' (n=' num2str(length(cond1)) ')'])
xlim([0 X(end)+1])
ylim(ylimDis)
xticks(X);
xticklabels(XTick);
% set(gca,'YTick',[0 0.2 0.4])


% Orientation %%%
angfig=figure(2);
set(angfig, 'Position', [400 50 1200 450])

subplot(1, 2, 2)
hold on
% a1 = plot(repmat(X, length(cond2), 1)', Y_ang(cond2,:)', 'Marker', '*', 'LineWidth',2);
% a2 = plot(repmat(X, length(cond1), 1)', Y_ang(cond1,:)', 'Marker', '*', 'LineWidth',2);
a1 = plot(repmat(X, length(cond2), 1)', Y_ang(cond2,:)', 'Color', cond2Color, 'Marker', '*', 'LineWidth',2);
a2 = plot(repmat(X, length(cond1), 1)', Y_ang(cond1,:)', 'Color', cond1Color, 'Marker', '*', 'LineWidth',2);
% set(a1, {'color'}, cond2Color)
% set(a2, {'color'}, cond1Color)
title(Y_ang_title)
set(gca, 'FontSize', 14)
% set(gca, 'FontSize', 18)
xlabel(x_label)
ylabel(Y_ang_ylabel)
% legend(names(cond1))
% legend([names(cond2) names(cond1)])
legend([a1(1) a2(1)], [cond2name ' (n=' num2str(length(cond2)) ')'], ...
  [cond1name ' (n=' num2str(length(cond1)) ')']) %{'cont', 'stim'})
xlim([0 X(end)+1])
% ylim([0 0.3])%25]) %4])
ylimAng = angfig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);
% set(gca,'YTick',[0 0.1 0.2])

subplot(1, 2, 1)
hold on
errorbar(X,conaavg,conastd, 'Color', cond2Color, 'LineWidth',2)
errorbar(X,stuaavg,stuastd, 'Color', cond1Color, 'LineWidth',2)
title(Y_ang_title)
set(gca, 'FontSize', 14)
% set(gca, 'FontSize', 18)
xlabel(x_label)
ylabel(Y_ang_ylabel)
legend([cond2name ' (n=' num2str(length(cond2)) ')'], ...
   [cond1name ' (n=' num2str(length(cond1)) ')'])
xlim([0 X(end)+1])
ylim(ylimAng)
xticks(X);
xticklabels(XTick);
% set(gca,'YTick',[0 0.1 0.2])

if(0)
    saveas(disfig,'timeNearObj_10min.tif')
    saveas(angfig,'orientToObj_10min.tif')
    %ranksum(Y_dis(cond1,3), Y_dis(cond2,3))
    %signrank(Y_dis(cond1,2), Y_dis(cond1,3))
end

%% for plotting number of bouts and bout length
%    (based on file from bout_analysis.m)
% Y_dis_title = ['Number of bouts (' num2str(round(Dis_te_frame/fps/60)) 'min)'];
% Y_dis_ylabel = 'Number of bouts';
% 
% Y_ang_title = ['Mean bout length (' num2str(round(Dis_te_frame/fps/60)) 'min)'];
% Y_ang_ylabel = 'Mean bout length (s)';
% Y_ang_title = ['Median bout length (' num2str(round(Dis_te_frame/fps/60)) 'min)];
% Y_ang_ylabel = 'Median bout length (s)';

% Y_dis_title = ['Time spent in quadrant 3 (' num2str(round(Dis_te_frame/fps/60)) 'min)'];
% Y_dis_ylabel = 'Fraction of time in quad 3';
% Y_ang_title = ['Time spent in quadrant 4 (' num2str(round(Dis_te_frame/fps/60)) 'min)'];
% Y_ang_ylabel = 'Fraction of time in quad 4';

Y_dis_title = ['Total number of SAPs (' num2str(round(Dis_te_frame/fps/60)) 'min)'];
Y_dis_ylabel = 'Number of SAPs';
Y_ang_title = ['SAPs near obj (norm; ' num2str(round(Dis_te_frame/fps/60)) 'min; '...
    'Rad = ' num2str(radius_cm) 'cm)'];
Y_ang_ylabel = 'Number of SAPs / time near';

close all
boutNumfig=figure(1);
set(boutNumfig, 'Position', [600 600 1200 450])

subplot(1, 2, 2)
hold on
% d1 = plot(repmat(X, length(cond2), 1)', Y_dis(cond2,:)', 'Marker', '*');
% d2 = plot(repmat(X, length(cond1), 1)', Y_dis(cond1,:)', 'Marker', '*');
d1 = plot(repmat(X, length(cond2), 1)', Y_dis(cond2,:)', 'Color', cond2Color, 'Marker', '*', 'LineWidth',2);
d2 = plot(repmat(X, length(cond1), 1)', Y_dis(cond1,:)', 'Color', cond1Color, 'Marker', '*', 'LineWidth',2);
title(Y_dis_title)
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel(Y_dis_ylabel)
% legend([names(cond2,:); names(cond1,:)], 'location', 'northwest')
legend([d1(1) d2(1)], {'cont', 'stim'})
xlim([0 X(end)+1])
ylim([0 90])
ylimDis = boutNumfig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);

subplot(1, 2, 1)
hold on
errorbar(X,condavg,condstd, 'Color', cond2Color, 'LineWidth',2)
errorbar(X,studavg,studstd, 'Color', cond1Color, 'LineWidth',2)
title(Y_dis_title)
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel(Y_dis_ylabel)
legend([cond2name ' (n=' num2str(length(cond2)) ')'], ...
    [cond1name ' (n=' num2str(length(cond1)) ')'])
xlim([0 X(end)+1])
ylim(ylimDis)
xticks(X);
xticklabels(XTick);



boutLenfig=figure(2);
set(boutLenfig, 'Position', [400 50 1200 450])

subplot(1, 2, 2)
hold on
% a1 = plot(repmat(X, length(cond2), 1)', Y_ang(cond2,:)', 'Marker', '*');%, ...
% a2 = plot(repmat(X, length(cond1), 1)', Y_ang(cond1,:)', 'Marker', '*');%, ...
a1 = plot(repmat(X, length(cond2), 1)', Y_ang(cond2,:)', 'Color', cond2Color, 'Marker', '*', 'LineWidth',2);
a2 = plot(repmat(X, length(cond1), 1)', Y_ang(cond1,:)', 'Color', cond1Color, 'Marker', '*', 'LineWidth',2);
title(Y_ang_title)
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel(Y_ang_ylabel)
% legend([names(cond2,:); names(cond1,:)], 'location', 'northwest')
legend([a1(1) a2(1)], {'cont', 'stim'})
xlim([0 X(end)+1])
% ylim([0 1])
ylimAng = boutLenfig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);

subplot(1, 2, 1)
hold on
errorbar(X,conaavg,conastd, 'Color', cond2Color, 'LineWidth',2)
errorbar(X,stuaavg,stuastd, 'Color', cond1Color, 'LineWidth',2)
title(Y_ang_title)
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel(Y_ang_ylabel)
legend([cond2name ' (n=' num2str(length(cond2)) ')'], ...
    [cond1name ' (n=' num2str(length(cond1)) ')'])
xlim([0 X(end)+1])
ylim(ylimAng)
xticks(X);
xticklabels(XTick);

if(0)
    saveas(boutNumfig,'boutNum_10min_nose.tif')
    saveas(boutLenfig,'boutLen_10min_nose.tif')
    
    saveas(boutNumfig,'SAPnum_10min.tif')
    saveas(boutLenfig,'SAPnear_10min_12cm.tif')
    saveas(boutLenfig,'SAPnear_10min_12cm_norm.tif')
end

%    saveas(boutLenfig,'boutLenMed_10min_nose.tif')
%saveas(boutNumfig,'boutNum_10min.tif')
%saveas(boutLenfig,'avgBoutLen_10min.tif')
%saveas(boutLenfig,'medBoutLen_10min.tif')

% saveas(boutNumfig,'quad1time_10min.tif')
% saveas(boutLenfig,'quad2time_10min.tif')
% saveas(boutNumfig,'quad3time_10min.tif')
% saveas(boutLenfig,'quad4time_10min.tif')

%% for plotting oneStats: totalDistRun, periph, centerTime
% Y_dis_title = ['Total distance run over ' num2str(round(Dis_te_frame/fps/60)) ' min'];
% Y_dis_ylabel = 'Total distance run (cm)';

% Y_dis_title = ['Time spent in periphery over ' num2str(round(Dis_te_frame/fps/60)) ' min'];
% Y_dis_ylabel = 'Fraction in periphery';

Y_dis_title = ['Time spent in center (radius=10cm; ' num2str(round(Dis_te_frame/fps/60)) 'min)'];
Y_dis_ylabel = 'Fraction of time spent in center';

oneStatFig=figure(1);
set(oneStatFig, 'Position', [600 600 1200 450])

subplot(1, 2, 2)
hold on
plot(repmat(X, length(cond2), 1)', Y_dis(cond2,:)', 'Marker', '*')%, ...
%     'Color', cond2Color, 'Marker', '*')
plot(repmat(X, length(cond1), 1)', Y_dis(cond1,:)', 'Marker', '*')%, ...
%     'Color', cond1Color, 'Marker', '*')
title(Y_dis_title)
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel(Y_dis_ylabel)
legend([names(cond2,:); names(cond1,:)], 'location', 'northwest')
xlim([0 X(end)+1])
%ylim([0 1])
ylimDis = oneStatFig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);

subplot(1, 2, 1)
hold on
errorbar(X,condavg,condstd, 'Color', cond2Color)
errorbar(X,studavg,studstd, 'Color', cond1Color)
title(Y_dis_title)
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel(Y_dis_ylabel)
legend([cond2name ' (n=' num2str(length(cond2)) ')'], ...
    [cond1name ' (n=' num2str(length(cond1)) ')'])
xlim([0 X(end)+1])
ylim(ylimDis)
xticks(X);
xticklabels(XTick);

if(0)
    saveas(oneStatFig,'totalDistanceRun_10min_nose.tif')
    saveas(oneStatFig,'timePeriph_10min_body.tif')
    saveas(oneStatFig,'centerTime_10min_nose.tif')
end


%% box plot of time near + orient for single day
currday = 3;

currYdis_cond1 = Y_dis(cond1,:);
currYdis_cond2 = Y_dis(cond2,:);

%%currYdis_cond2 = Y_dis(:,currday);
%%currYdis_cond1 = Y_dis(:,2);
%%currYdis_cond3 = Y_dis(:,8);
% currTitle      = ['Time near object: Day ' num2str(currday) ' (tail)'];
% currYlim       = [0 0.6];

%cond2 = [2 3 5 6 8 9];
%cond2 = [2 8 9];
%cond3 = [3 5 6];

%currYdis_cond1 = Y_dis(cond1,currday)-Y_dis(cond1,2);
%currYdis_cond2 = Y_dis(cond2,currday)-Y_dis(cond2,2);
%currYdis_cond3 = Y_dis(cond3,currday)-Y_dis(cond3,2);

% currTitle = ['Time near object: Day ' num2str(currday) ' (tail; norm)'];
currTitle = ['Time near object: cont (tail)'];
currYlim  = [-0.0 0.2];


close all
boxPlot_oneDay = figure(1);
% boxplot([currYdis_cond1; currYdis_cond2], ...; currYdis_cond3], ...
%     [zeros(1,length(currYdis_cond1)), ones(1,length(currYdis_cond2))])%, repmat(2, 1, length(currYdis_cond3))])
boxplot(currYdis_cond2)%, ...; currYdis_cond3], ...
%    [zeros(1,length(currYdis_cond1)), ones(1,length(currYdis_cond2))])%, repmat(2, 1, length(currYdis_cond3))])
%boxplot(currYdis_cond1)
set(gca,'FontSize',16)
%set(gca,'XTickLabel',{'saline','6OHDA'})%,'6OHDA2'})
%set(gca,'XTickLabel',{cond1name,cond2name})
set(gca, 'XTickLabel', {'H1','H2','1','2','3','4','5','6','7','8','9','10'})
title(currTitle)
%xlabel('Condition')
xlabel('Day')
ylabel('Time spent near object (frac)')
ylim(currYlim)

if(0)
    saveas(boxPlot_oneDay,'timeNearObj_10min_box_combine.tif')
    saveas(boxPlot_oneDay,'timeNearObj_10min_tail_box_allDays_stim.tif')
    saveas(boxPlot_oneDay, ['timeNearObj_10min_tail_box_N' num2str(currday-2) '_norm.tif'])
end

%% cumulative histogram of bouts over one session
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Capoeira_DLC
load('PokesApproaches.mat')
binsX = 0:15*60:10000;

close all
cBout_fig = figure(1);
set(cBout_fig, 'Position', [700 300 570 450])
hold on
for mousei=1:length(Mice)
    
    for filei=3
        
        currBouts = Mice(mousei).(['Pokes_Day' num2str(filei)]);
        h = histcounts(currBouts,binsX);
        cumsum_bouts = cumsum(h);
        
        if(strcmp(Mice(mousei).novelty,'S'))
            plot(1:length(cumsum_bouts),cumsum_bouts,'color',cond1Color)
        else
            plot(1:length(cumsum_bouts),cumsum_bouts,'color',cond2Color)
        end
    end
end
% saveas(gca,'Capoeira_cBouts_N1.tif')

%%

mouse_name = {'Bottom', 'Charm', 'Strange', 'Up'};
timeStat_all = zeros(12, 4, length(mouse_name));

for mouse = 1:length(mouse_name)
    cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Quarky_conf_DLC/' ...
        mouse_name{mouse}])
    timeStat = readtable('TimeStatistic.csv');
    
    if(strcmp(mouse_name{mouse}, 'Aldehyde'))
        timeStat_all(2:12,:,mouse) = timeStat{:,end};
    else
        timeStat_all(:,:,mouse) = timeStat{:,2:end};
    end
end

XTick = {'1' '2' '3' '4' '5' '6' '7' '8' 'N1' 'N2' 'N3' 'N4'};
X     = 1:size(timeStat_all,1);


currMouse = 4;

baseline  = 8; %5:8
novelty   = 9; %9:12
novelObjs = 1:2;
familObjs = 3:4;

avgStat_all = zeros(2,4,length(mouse_name));
avgStat_avg = zeros(4,length(mouse_name));
for mouse = 1:length(mouse_name)
    avgStat_all(:,:,mouse) = [mean(timeStat_all(baseline,:,mouse),1); ...
                              mean(timeStat_all(novelty,:,mouse),1)];
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
ylim([0 max(max(timeStat_all(baseline(1):novelty(end),:,currMouse)))])
ylabel('Time around objects (s)')
title(['Time around objects (', mouse_name{currMouse}, '; r=', num2str(radius_cm), 'cm)'], ...
    'Interpreter', 'none')
legend([willBeNovel willBeFamil], ...
    {'Will be novel', 'Will be familiar'}, 'Location', 'Northwest')
set(gca, 'FontSize', 14, 'XTick', 1:3, ...
    'XTickLabel', {['sessions ' num2str(baseline(1)) ':' num2str(baseline(end))], ...
    'familiar', 'novel'})

%cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Quarky_conf_DLC/'])% ...
%    mouse_name{currMouse}])
%saveas(avgStatOne,['avgStat_' mouse_name{currMouse} '_8-9.tif'])
%disp('saved')

%%
disfig=figure(1);
set(disfig, 'Position', [600 600 1200 450])

subplot(1, 2, 2)
hold on
plot(repmat(X, size(Y_dis, 1), 1)', Y_dis', ...
    'Color', 'k', 'Marker', '*')
title(['Time spent at obj (' num2str(round(Dis_te_frame/fps/60))  ...
    'min); Rad = ' num2str(radius_cm) 'cm'])
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel('Fraction of time near obj')
% legend(names)
xlim([0 X(end)+1])
ylimDis = disfig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);


angfig=figure(2);
set(angfig, 'Position', [400 50 1200 450])

subplot(1, 2, 2)
hold on
plot(repmat(X, size(Y_ang, 1), 1)', Y_ang', ...
    'Color', 'k', 'Marker', '*')
title(['Orientation to obj (' num2str(round(Dis_te_frame/fps/60)) ...
    'min); Deg = +-' num2str(angle_radius) char(176)])
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel('Fraction of time oriented to obj')
% legend(names)
xlim([0 X(end)+1])
ylimAng = angfig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);

saveas(disfig,'timeNearObj_10min_combine.tif')
saveas(angfig,'orientToObj_10min_combine.tif')