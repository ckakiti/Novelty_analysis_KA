% Analyze time spent in periphery, center, or quadrants of arena

clear
close all
clc

Config_NovAna
cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Rim_KO_DLC/')
run('MiceIndex')

% only analyze first 10 minutes of video
start_min=0.5;
end_min=start_min+10;
startframe=start_min.*60.*fps;
endframe=end_min.*60.*fps;
startframe=round(startframe);
endframe=round(endframe);

% for smoothing
Swindow=40;
Xtime=startframe:endframe;

files    = dir;
whichDir = [files.isdir];
nameDir  = files(whichDir);
nameDir  = {nameDir.name};
nameDir(ismember(nameDir,{'.','..','Retrain_Sep17'})) = [];
days     = 4;

% time spent in periphery, within center circle, or different quadrants
periphTime = zeros(length(nameDir),days);
centerTime = zeros(length(nameDir),days);
quad1Time  = zeros(length(nameDir),days);
quad2Time  = zeros(length(nameDir),days);
quad3Time  = zeros(length(nameDir),days);
quad4Time  = zeros(length(nameDir),days);

for mousei = 1:length(nameDir)
    cd([nameDir{1,mousei}, '/Analyzed_Data'])
    load('Arena_Obj_Pos.mat','arena')
    
    matFiles = dir('*rgb.mat');
    matFiles = cat(1, matFiles.name);
    
    arenaCenterPos = [(arena(1,1)+arena(1,3))/2 (arena(1,2)+arena(1,4))/2];
    radius_cm_center = 10;
    
    for filei = 1:size(matFiles,1)
        load(matFiles(filei,:))

        % smooth nose position
        Dis=Labels(Xtime,2:3); % pix
        SDis=smoothdata(Dis,'rloess',Swindow);
        
        % calculate time spent within center (circle at center)
        SDis_inCenter  = sqrt((arenaCenterPos(1)-SDis(:,1)).^2 ...
                             +(arenaCenterPos(2)-SDis(:,2)).^2)/ppc;
        centerTime(mousei,filei) = sum(SDis_inCenter<radius_cm_center)/length(Xtime);
        
        
        % calculate time spent in each quadrant 
        %   1=topleft; 2=topright; 3=bottomright; 4=bottomleft;
        SDis_in1 = length( intersect( find(SDis(:,1)<arenaCenterPos(1)), ...
                                      find(SDis(:,2)<arenaCenterPos(2))) )/length(Xtime);
                                  
        SDis_in2 = length( intersect( find(SDis(:,1)>=arenaCenterPos(1)), ...
                                      find(SDis(:,2)<arenaCenterPos(2))) )/length(Xtime);
                                  
        SDis_in3 = length( intersect( find(SDis(:,1)>=arenaCenterPos(1)), ...
                                      find(SDis(:,2)>=arenaCenterPos(2))) )/length(Xtime);
                                  
        SDis_in4 = length( intersect( find(SDis(:,1)<arenaCenterPos(1)), ...
                                      find(SDis(:,2)>=arenaCenterPos(2))) )/length(Xtime);
                                  
        quad1Time(mousei,filei) = SDis_in1;
        quad2Time(mousei,filei) = SDis_in2;
        quad3Time(mousei,filei) = SDis_in3;
        quad4Time(mousei,filei) = SDis_in4;
    end
    
    cd ../..
end
disp('end')

%% sanity plotting: center circle
inCenter = SDis_inCenter<radius_cm_center;

close all

fig1 = figure(1);
hold on
scatter(SDis(:,1), SDis(:,2), 8, 'filled');
scatter(SDis(inCenter,1), SDis(inCenter,2), 8, 'r', 'filled');
rectangle('Position',[arena(1,1),arena(1,2),arena(1,3)-arena(1,1),arena(1,4)-arena(1,2)],...
    'EdgeColor','r','linewidth',1)

th = 0:pi/50:2*pi;
x  = arenaCenterPos(1);
y  = arenaCenterPos(2);
xunit = radius_cm_center*ppc * cos(th) + x;
yunit = radius_cm_center*ppc * sin(th) + y;
plot(xunit, yunit)

set(gca,'ydir','reverse')
title(['TimeInCenter_radius=' num2str(radius_cm_center) ...
    'cm_Mouse' nameDir{1,mousei} '_Day' num2str(filei)], 'Interpreter', 'none')
xlim([0 video_xlen])
ylim([0 video_ywid])

if(0)
    saveas(fig1,['TimeInCenter_radius=' num2str(radius_cm_center) ...
        'cm_Mouse' nameDir{1,mousei} '_Day' num2str(filei) '.png'])
end

%% sanity plotting: quadrants
SDis_in1pts = intersect( find(SDis(:,1)<arenaCenterPos(1)), ...
                         find(SDis(:,2)<arenaCenterPos(2)));

SDis_in2pts = intersect( find(SDis(:,1)>=arenaCenterPos(1)), ...
                         find(SDis(:,2)<arenaCenterPos(2)));

SDis_in3pts = intersect( find(SDis(:,1)>=arenaCenterPos(1)), ...
                         find(SDis(:,2)>=arenaCenterPos(2)));

SDis_in4pts = intersect( find(SDis(:,1)<arenaCenterPos(1)), ...
                         find(SDis(:,2)>=arenaCenterPos(2)));

close all
fig2 = figure(2);
hold on
q1 = scatter(SDis(SDis_in1pts,1), SDis(SDis_in1pts,2), 8, 'r', 'filled');
q2 = scatter(SDis(SDis_in2pts,1), SDis(SDis_in2pts,2), 8, 'g', 'filled');
q3 = scatter(SDis(SDis_in3pts,1), SDis(SDis_in3pts,2), 8, 'b', 'filled');
q4 = scatter(SDis(SDis_in4pts,1), SDis(SDis_in4pts,2), 8, 'c', 'filled');
rectangle('Position',[arena(1,1),arena(1,2),arena(1,3)-arena(1,1),arena(1,4)-arena(1,2)],...
    'EdgeColor','r','linewidth',1)

set(gca,'ydir','reverse')
title(['TimeInQuadrants_Mouse' nameDir{1,mousei} ...
    '_Day' num2str(filei)], 'Interpreter', 'none')
xlim([0 video_xlen])
ylim([0 video_ywid])
legend([q1 q2 q3 q4], {'Quad1', 'Quad2', 'Quad3', 'Quad4'})

% saveas(fig2,['TimeInQuads_Mouse' nameDir{1,mousei} '_Day' num2str(filei) '.png'])

%%
numOrLen = cell(2*length(nameDir),1);
numOrLen(1:length(nameDir)) = {'centerTime'};
numOrLen(length(nameDir)+1:end) = {'centerTime'};
TableLabels = [cat(2,nameDir,nameDir)', numOrLen];
TableInfo   = [centerTime;centerTime];
Table       = [cell2table(TableLabels) array2table(TableInfo)];
disp('end')
%writetable(Table,'areaAnalysis_nose_center.csv');
%%
numOrLen = cell(2*length(nameDir),1);
numOrLen(1:length(nameDir)) = {'quad1Time'};
numOrLen(length(nameDir)+1:end) = {'quad2Time'};
TableLabels = [cat(2,nameDir,nameDir)', numOrLen];
TableInfo   = [quad1Time;quad2Time];
Table       = [cell2table(TableLabels) array2table(TableInfo)];

numOrLen2 = cell(2*length(nameDir),1);
numOrLen2(1:length(nameDir)) = {'quad3Time'};
numOrLen2(length(nameDir)+1:end) = {'quad4Time'};
TableLabels2 = [cat(2,nameDir,nameDir)', numOrLen2];
TableInfo2   = [quad3Time;quad4Time];
Table2       = [cell2table(TableLabels2) array2table(TableInfo2)];

% writetable(Table,'areaAnalysis_nose_quad1-2.csv');
% writetable(Table2,'areaAnalysis_nose_quad3-4.csv');






