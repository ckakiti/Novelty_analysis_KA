clear
close all
clc

Config_NovAna

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry2_DLC
%cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/NewHope-ROTJ/')
%cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_6OHDA_DLC/')
%cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/CvsS_180831_DLC/')

periphStat = 0;

if(periphStat)
    timeStat = readtable('TimeStatistic_periph.csv');
    timeStat2 = timeStat{:,3:end};

    names = timeStat{1:height(timeStat),1}';
    Y_dis = timeStat2(1:height(timeStat), :);
    Y_ang = Y_dis;
else
    timeStat = readtable('TimeStatistic.csv');%'boutAnalysis_body.csv');
    timeStat2 = timeStat{:,3:end};

    names = timeStat{1:height(timeStat)/2,1}';
    Y_dis = timeStat2(1:height(timeStat)/2, :); %1:8
    Y_ang = timeStat2(height(timeStat)/2+1:height(timeStat),:); %9:16
end

XTick = {'H1' 'H2' 'N1' 'N2' 'N3' 'N4' 'N5' 'N6' 'N7' 'N8' 'N9' 'N10'};% 'N11' 'N12'};

% cont = 1:4; % 1-4 rows are contextual novelty (CvsS); 1-3 quarky
% stim = 5:8; % 5:8 rows are stimulus novelty (CvsS); 4-6 quarky
cond2name = 'cont';%'6OHDA';
cond1name = 'stim';%'saline';
cond2 = [2 4 6];%1:6;%[2 4 6];%[2 5 6];%[];%1:4; %3:4; %[1:3 5:7 9:11]; %[2 3]; %1:3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
cond1 = [1 3 5];%7:12;%[1 3 5];%[1 3 4];%1:7;%5:8;%1:2; %[4 8 12];       %[1 4]; %4:6; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cond2Color = [0.5 0.0 0.5]; 
cond1Color = [1.0 0.5 0.0]; 

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

%% for plotting time around object and orientation to object
disfig=figure(1);
set(disfig, 'Position', [600 600 1200 450])

subplot(1, 2, 2)
hold on
plot(repmat(X, length(cond2), 1)', Y_dis(cond2,:)', ...
    'Color', cond2Color, 'Marker', '*')
plot(repmat(X, length(cond1), 1)', Y_dis(cond1,:)', ...
    'Color', cond1Color, 'Marker', '*')
title(['Time spent at obj (' num2str(round(Dis_te_frame/fps/60))  'min); Rad = ' num2str(radius_cm) 'cm'])
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel('Fraction of time near obj')
% legend(names)
xlim([0 X(end)+1])
ylimDis = disfig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);

subplot(1, 2, 1)
hold on
errorbar(X,condavg,condstd, 'Color', cond2Color)
errorbar(X,studavg,studstd, 'Color', cond1Color)
title(['Time spent at obj (' num2str(round(Dis_te_frame/fps/60))  'min); Rad = ' num2str(radius_cm) 'cm'])
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel('Fraction of time near obj')
legend([cond2name ' (n=' num2str(length(cond2)) ')'], ...
    [cond1name ' (n=' num2str(length(cond1)) ')'])
xlim([0 X(end)+1])
ylim(ylimDis)
xticks(X);
xticklabels(XTick);


angfig=figure(2);
set(angfig, 'Position', [400 50 1200 450])

subplot(1, 2, 2)
hold on
plot(repmat(X, length(cond2), 1)', Y_ang(cond2,:)', ...
    'Color', cond2Color, 'Marker', '*')
plot(repmat(X, length(cond1), 1)', Y_ang(cond1,:)', ...
    'Color', cond1Color, 'Marker', '*')
title(['Orientation to obj (' num2str(round(Dis_te_frame/fps/60)) 'min); Deg = +-' num2str(angle_radius) char(176)])
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel('Fraction of time oriented to obj')
% legend(names)
xlim([0 X(end)+1])
ylimAng = angfig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);

subplot(1, 2, 1)
hold on
errorbar(X,conaavg,conastd, 'Color', cond2Color)
errorbar(X,stuaavg,stuastd, 'Color', cond1Color)
title(['Orientation to obj (' num2str(round(Dis_te_frame/fps/60)) 'min); Deg = +-' num2str(angle_radius) char(176)])
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel('Fraction of time oriented to obj')
legend([cond2name ' (n=' num2str(length(cond2)) ')'], ...
    [cond1name ' (n=' num2str(length(cond1)) ')'])
xlim([0 X(end)+1])
ylim(ylimAng)
xticks(X);
xticklabels(XTick);

%saveas(disfig,'timeNearObj_10min')
%saveas(angfig,'orientToObj_10min')
%saveas(disfig,'timeNearObj_10min.tif')
%saveas(angfig,'orientToObj_10min.tif')

%% for plotting number of bouts and bout length
boutNumfig=figure(1);
set(boutNumfig, 'Position', [600 600 1200 450])

subplot(1, 2, 2)
hold on
plot(repmat(X, length(cond2), 1)', Y_dis(cond2,:)', ...
    'Color', cond2Color, 'Marker', '*')
plot(repmat(X, length(cond1), 1)', Y_dis(cond1,:)', ...
    'Color', cond1Color, 'Marker', '*')
title(['Number of bouts in ' num2str(round(Dis_te_frame/fps/60))  'min; Rad = ' num2str(radius_cm) 'cm'])
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel('Number of bouts')
% legend(names)
xlim([0 X(end)+1])
ylimDis = boutNumfig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);

subplot(1, 2, 1)
hold on
errorbar(X,condavg,condstd, 'Color', cond2Color)
errorbar(X,studavg,studstd, 'Color', cond1Color)
title(['Number of bouts in ' num2str(round(Dis_te_frame/fps/60))  'min; Rad = ' num2str(radius_cm) 'cm'])
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel('Number of bouts')
legend([cond2name ' (n=' num2str(length(cond2)) ')'], ...
    [cond1name ' (n=' num2str(length(cond1)) ')'])
xlim([0 X(end)+1])
ylim(ylimDis)
xticks(X);
xticklabels(XTick);


Y_ang = Y_ang/fps;
conaavg=mean(Y_ang(cond2,:));
stuaavg=mean(Y_ang(cond1,:));
conastd=std(Y_ang(cond2,:));
stuastd=std(Y_ang(cond1,:));

boutLenfig=figure(2);
set(boutLenfig, 'Position', [400 50 1200 450])

subplot(1, 2, 2)
hold on
plot(repmat(X, length(cond2), 1)', Y_ang(cond2,:)', ...
    'Color', cond2Color, 'Marker', '*')
plot(repmat(X, length(cond2), 1)', Y_ang(cond1,:)', ...
    'Color', cond1Color, 'Marker', '*')
title(['Average bout length (' num2str(round(Dis_te_frame/fps/60)) 'min); Rad = ' num2str(radius_cm) 'cm'])
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel('Average bout length (s)')
% legend(names)
xlim([0 X(end)+1])
ylimAng = boutLenfig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);

subplot(1, 2, 1)
hold on
errorbar(X,conaavg,conastd, 'Color', cond2Color)
errorbar(X,stuaavg,stuastd, 'Color', cond1Color)
title(['Average bout length (' num2str(round(Dis_te_frame/fps/60)) 'min); Rad = ' num2str(radius_cm) 'cm'])
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel('Average bout length (s)')
legend([cond2name ' (n=' num2str(length(cond2)) ')'], ...
    [cond1name ' (n=' num2str(length(cond1)) ')'])
xlim([0 X(end)+1])
ylim(ylimAng)
xticks(X);
xticklabels(XTick);

%saveas(boutNumfig,'boutNum_10min_body.tif')
%saveas(boutLenfig,'boutLen_10min_body.tif')

%%
periphFig=figure(1);
set(periphFig, 'Position', [600 600 1200 450])

subplot(1, 2, 2)
hold on
plot(repmat(X, length(cond2), 1)', Y_dis(cond2,:)', ...
    'Color', cond2Color, 'Marker', '*')
plot(repmat(X, length(cond1), 1)', Y_dis(cond1,:)', ...
    'Color', cond1Color, 'Marker', '*')
title(['Time spent in periphery over ' num2str(round(Dis_te_frame/fps/60)) ' min'])
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel('Fraction in periphery')
% legend(names)
xlim([0 X(end)+1])
ylim([0 1])
ylimDis = periphFig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);

subplot(1, 2, 1)
hold on
errorbar(X,condavg,condstd, 'Color', cond2Color)
errorbar(X,studavg,studstd, 'Color', cond1Color)
title(['Time spent in periphery over ' num2str(round(Dis_te_frame/fps/60)) ' min'])
set(gca, 'FontSize', 14)
xlabel('Training day')
ylabel('Fraction in periphery')
legend([cond2name ' (n=' num2str(length(cond2)) ')'], ...
    [cond1name ' (n=' num2str(length(cond1)) ')'])
xlim([0 X(end)+1])
ylim(ylimDis)
xticks(X);
xticklabels(XTick);

%saveas(periphFig,'timePeriph_10min_body.tif')

%%

mouse_name = {'Bottom', 'Charm', 'Strange', 'Up'};
timeStat_all = zeros(12, 4, length(mouse_name));

for mouse = 1:length(mouse_name)
    cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Quarky_conf_DLC/' ...
        mouse_name{mouse}])
    timeStat = readtable('TimeStatistic.csv');
    
    if(strcmp(mouse_name{mouse}, 'Aldehyde'))
        timeStat_all(2:12,:,mouse) = timeStat{:,2:end};
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
legend([willBeNovel willBeFamil], {'Will be novel', 'Will be familiar'}, 'Location', 'Northwest')
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
title(['Time spent at obj (' num2str(round(Dis_te_frame/fps/60))  'min); Rad = ' num2str(radius_cm) 'cm'])
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
title(['Orientation to obj (' num2str(round(Dis_te_frame/fps/60)) 'min); Deg = +-' num2str(angle_radius) char(176)])
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