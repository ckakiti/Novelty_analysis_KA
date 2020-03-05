%% Divide session into smaller time bins (updated version of plotting4.m)

clear
clc
close all
cd /home/alex/Programs/Novelty_analysis_KA
Config_NovAna;

currMouse = 'Stability';
disp(['Current mouse: ' currMouse])

cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/DRILLS_DLC/' currMouse])
pathname = cd;
PathRoot=[pathname '/'];
filelist=dir([PathRoot,'*' videoname_format(end-3:end)]);
flen = length(filelist);
cd Analyzed_Data_1obj_30min
load('Arena_Obj_Pos.mat');

durTotal = 30; %10 %30 % duration of analysis (min)
bins     = 5;  %5      % number of bins to divide into
disp(['Duration of analysis: ' num2str(durTotal) 'min; Bins: ' num2str(bins)])

Dis_ts_frame=500;
Dis_te_frame=durTotal.*60.*fps+Dis_ts_frame;

% if lego is put in arena after 5 min, don't analyze that first 5 min
shift_time = input('Exclude first 5 min? 0/1: ');

% %
if(~exist('divide_sessions','dir'))
    mkdir('divide_sessions')
end

timeNear_split = zeros(flen, bins);
orient_split = zeros(flen, bins);

for fiter = 1:flen
    load([filelist(fiter).name(1:end-4) '.mat'])
    vn = filelist(fiter).name(1:(length(currMouse)+7));
    disp(vn)
    
    % if lego is put in arena after 5 min, don't analyze that first 5 min
    if(shift_time==1)
        timeShift = 5000;
        disp(timeShift)
        
        startframe = Dis_ts_frame+timeShift;
        endframe   = min(Dis_te_frame+timeShift, size(Labels,1));
    else
        startframe = Dis_ts_frame;
        endframe   = min(Dis_te_frame, size(Labels,1));
    end
    
    % specify number of bins and frame number of bin edge
%     bins = 5;
    bin_frame = linspace(startframe, endframe, bins+1);
    
    for k = 1:bins
        % calculate time spent near object within each time bin
        timeNear_split(fiter,k) = sum(Labels(bin_frame(k):bin_frame(k+1),21))/(bin_frame(k+1)-bin_frame(k));
        orient_split(fiter,k)   = sum(Labels(bin_frame(k):bin_frame(k+1),23))/(bin_frame(k+1)-bin_frame(k));
        
%         % plot trajectory within each time bin
%         Trafigure=figure('visible','off');
%         scatter(Labels(bin_frame(k):bin_frame(k+1),14), ...
%                 Labels(bin_frame(k):bin_frame(k+1),15), 6, 'filled');
%         rectangle('Position',[arena(fiter,1),arena(fiter,2),...
%             arena(fiter,3)-arena(fiter,1),...
%             arena(fiter,4)-arena(fiter,2)],'EdgeColor','r','linewidth',8)
%         rectangle('Position',[obj(fiter,1),obj(fiter,2),...
%             obj(fiter,3)-obj(fiter,1),obj(fiter,4)-obj(fiter,2)],...
%             'EdgeColor','r','linewidth',2)
%         
%         hold on
%         th = 0:pi/50:2*pi;
%         x  = obj_center(fiter,1);
%         y  = obj_center(fiter,2);
%         xunit = radius * cos(th) + x;
%         yunit = radius * sin(th) + y;
%         plot(xunit, yunit,'r--','linewidth',3)
%         
%         set(gca,'ydir','reverse')
%         title(['Trajectory ' vn ': bin ' num2str(k) '/' num2str(bins) ...
%             ', radius=' num2str(radius_cm) 'cm'], 'Interpreter', 'none');
%         xlim([0 video_xlen+50]);
%         ylim([0 video_ywid]);
%         set(Trafigure, 'position', [0 0 1200 900]);
%         
%         cd divide_sessions
%         saveas(Trafigure,['Trajectory_' vn '_bin' num2str(k) '.tif'])
%         cd ..
%         
%         close all
    end
end

% if(0)
cd divide_sessions
save([currMouse '_timeNearOrient_split'], 'timeNear_split', 'orient_split', 'startframe', 'endframe');
disp('saved')
% end

disp('end')

%% gather timeNear_split for N1 across all mice
clear
clc
close all

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/DRILLS_DLC/
whichDay = 3; %N1

folderpath = cd;
folderd = dir(folderpath);
isub = [folderd(:).isdir];
foldernames = {folderd(isub).name}'; 
foldernames(ismember(foldernames,{'.','..','temp'})) = []; 
folderlen=length(foldernames);

timeNear_splitAll = [];
orient_splitAll   = [];

for fiter = 1:folderlen
    cd(foldernames{fiter})
    cd Analyzed_Data_1obj_30min
    cd divide_sessions
    
    matFile = dir('*timeNearOrient_split.mat');
    load(matFile.name)
    
    timeNear_splitAll = [timeNear_splitAll; timeNear_split(whichDay,:)];
    orient_splitAll   = [orient_splitAll; orient_split(whichDay,:)];
    
    cd ../../..
end

DisOrAng = [repmat({'Distance'},folderlen,1); repmat({'Angle'},folderlen,1)];
Table  = table(cat(1,foldernames,foldernames),DisOrAng,[timeNear_splitAll;orient_splitAll]);

if(0)
    writetable(Table,'TimeStatistic_N1_split_30min.csv');
end

%%
cd divide_sessions
currMat = dir('*timeNearOrient_split.mat');
load(currMat.name)

whichSessions = [2 3];
X = 1:size(timeNear_split,2);
XTick = {'0-2','2-4','4-6','6-8','8-10'};
cond2Color = [0.5 0.0 0.5];
cond1Color = [1.0 0.5 0.0];

Y_dis_title = ['Time spent at obj: divided (' num2str(round(Dis_te_frame/fps/60))  ...
    'min); Rad = ' num2str(radius_cm) 'cm'];
Y_ang_title = ['Orientation to obj: divided (' num2str(round(Dis_te_frame/fps/60)) ...
    'min); Deg = +-' num2str(angle_radius) char(176)];
Y_dis_ylabel = 'Time near object (frac)';
Y_ang_ylabel = 'Orientation to object (frac)';

close all 


% Time near obj %%%
disfig=figure(1);
set(disfig, 'Position', [600 600 600 450])

hold on
d1 = plot(X, timeNear_split(whichSessions(1),:), 'Color', [0 0 0], 'Marker', '*', 'LineWidth',2);
d2 = plot(X, timeNear_split(whichSessions(2),:), 'Color', cond1Color, 'Marker', '*', 'LineWidth',2);
title(Y_dis_title)
set(gca, 'FontSize', 14)
xlabel('Time bin')
ylabel(Y_dis_ylabel)
legend([d1(1) d2(1)], {'H2', 'N1'})
xlim([0 X(end)+1])
ylim([0 0.45])
ylimDis = disfig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);



% Orientation %%%
angfig=figure(2);
set(angfig, 'Position', [400 50 600 450])

hold on
a1 = plot(X, orient_split(whichSessions(1),:), 'Color', [0 0 0], 'Marker', '*', 'LineWidth',2);
a2 = plot(X, orient_split(whichSessions(2),:), 'Color', cond1Color, 'Marker', '*', 'LineWidth',2);
title(Y_ang_title)
set(gca, 'FontSize', 14)
xlabel('Time bin')
ylabel(Y_ang_ylabel)
legend([a1(1) a2(1)], {'H2', 'N1'})
xlim([0 X(end)+1])
ylim([0 0.225])
ylimAng = angfig.CurrentAxes.YLim;
xticks(X);
xticklabels(XTick);

if(0)
    saveas(disfig,[currMouse '_timeNear_split.tif'])
    saveas(angfig,[currMouse '_orient_split.tif'])
    disp('saved')
end