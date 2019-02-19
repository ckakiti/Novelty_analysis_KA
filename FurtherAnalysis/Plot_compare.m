clear
close all
clc

Config_NovAna

%cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/7day_preexposure_combine
cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/NewHope-ROTJ/')
timeStat = readtable('boutAnalysis_body.csv'); %'TimeStatistic_body.csv');
timeStat2 = timeStat{:,3:end};

names = timeStat{1:height(timeStat)/2,1}';
XTick = {'H1' 'H2' 'N1' 'N2' 'N3' 'N4' 'N5' 'N6' 'N7' 'N8' 'N9' 'N10'};% 'N11' 'N12'};

% cont = 1:4; % 1-4 rows are contextual novelty (CvsS); 1-3 quarky
% stim = 5:8; % 5:8 rows are stimulus novelty (CvsS); 4-6 quarky
cond2name = 'cont';
cond1name = 'stim';
cond2 = 1:6; %1:4; %[2 3]; %[2 4 6 8 10 12]; %1:3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
cond1 = 7:12;%5:8; %[1 4]; %[1 3 5 7 9 11]; %4:6; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cond2Color = [0.5 0.0 0.5];
cond1Color = [1.0 0.5 0.0];

Y_dis = timeStat2(1:height(timeStat)/2, :); %1:8
Y_ang = timeStat2(height(timeStat)/2+1:height(timeStat),:); %9:16

% Y_dis=[0.0226666667	0.0091111111	0.0299444444	0.0453333333	0.0821666667	0.0513333333	0.0631111111	0.0626666667	0.0570555556	0.036	0.0305	0.046	0.0451666667	0.0636111111;
% 0.0188333333	0.0123333333	0.1467222222	0.0843333333	0.0858333333	0.0928333333	0.0566666667	0.0256111111	0.0645	0.086	0.0986111111	0.0602777778	0.03	0.0829444444;
% 0.0173888889	0.0293333333	0.1832222222	0.0732777778	0.0986666667	0.1905	0.1973333333	0.1200555556	0.1185	0.1745555556	0.0642777778	0.1025555556	0.056	0.0692222222;
% 0.0347777778	0.0192777778	0.0198888889	0.0082222222	0.0106111111	0.0062777778	0.0154444444	0.0179444444	0.0321111111	0.0382777778	0.0345555556	0.0753333333	0.0595555556	0.0344444444;
% 0.0230555556	0.0540555556	0.0050555556	0.0026666667	0.0035555556	0.0023333333	0.012	0.0062222222	0.0059444444	0.0038333333	0.0127777778	0.0079444444	0.0064444444	0.0338333333;
% 0.0262222222	0.0075555556	0.0262222222	0.0363888889	0.0221111111	0.0277222222	0.044	0.0358333333	0.0332777778	0.0229444444	0.0201111111	0.0169444444	0.0443333333	0.0287222222;];
% 
% Y_ang=[0.0434444444	0.062	0.0801111111	0.0912222222	0.0608333333	0.0642222222	0.0678888889	0.0637222222	0.0866111111	0.0537777778	0.0646111111	0.0901111111	0.0873333333	0.0829444444;
% 0.0255	0.0382777778	0.0921666667	0.0612222222	0.082	0.0754444444	0.0310555556	0.0337222222	0.0947222222	0.0555555556	0.1042222222	0.0876666667	0.0349444444	0.0454444444;
% 0.0304444444	0.0453888889	0.0993333333	0.0498888889	0.0628888889	0.0956666667	0.1246111111	0.0761666667	0.0931111111	0.0952777778	0.0612222222	0.0818888889	0.0624444444	0.1145;
% 0.0441111111	0.0341111111	0.0807777778	0.0593333333	0.0346111111	0.037	0.0398333333	0.0391111111	0.0552222222	0.0575	0.05	0.055	0.0487222222	0.0451111111;
% 0.0668888889	0.0798888889	0.1758888889	0.0786666667	0.104	0.1104444444	0.0599444444	0.0604444444	0.0621111111	0.0788888889	0.0568333333	0.096	0.1189444444	0.067;
% 0.0312222222	0.0401666667	0.0814444444	0.0788333333	0.0354444444	0.0331666667	0.0443888889	0.0391111111	0.0775555556	0.0762777778	0.0471666667	0.0445	0.0675	0.0480555556;];

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
plot(repmat(X, length(cond2), 1)', Y_ang(cond1,:)', ...
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