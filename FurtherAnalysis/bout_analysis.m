% Analyze bouts of interaction with object

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

% number of bouts and bout length
boutNum = zeros(length(nameDir),days);
boutLen = zeros(length(nameDir),days);

for mousei = 1:length(nameDir)
    cd([nameDir{1,mousei}, '/Analyzed_Data'])
    load('Arena_Obj_Pos.mat')
    
    matFiles = dir('*rgb.mat');
    matFiles = cat(1, matFiles.name);
    
    for filei = 1:size(matFiles,1)
        load(matFiles(filei,:))
        
        DisX=Labels(startframe:endframe,2); % pix, nose
        DisY=Labels(startframe:endframe,3); % pix, nose
        Dis=sqrt((obj_center(filei,1)-DisX).^2+(obj_center(filei,2)-DisY).^2)/ppc; % cm

        SDis=smoothdata(Dis,'rloess',Swindow);
        SDis_inRad = double(SDis<radius_cm);
        
        % find frames where mouse crosses within radius specified by Config_NovAna
%         crossTmp  = crossing(Labels(startframe:endframe,21), [], 0.5);
        crossTmp = crossing(SDis_inRad, [], 0.5);
        if(mod(length(crossTmp),2))
            crossTmp = crossTmp(1:end-1);
        end
        cross_In  = crossTmp(1:2:end);
        cross_Out = crossTmp(2:2:end);
        
        % find first increase of velocity before each frame number in crossTmp
        %    diffXY  = diff(Labels(:,2:3));
        %    accelXY = diff(diffXY);
        
        if(~isempty(cross_In))
            boutNum(mousei,filei) = length(cross_In);
            boutLen(mousei,filei) = mean(cross_Out-cross_In)/fps; % seconds
        end
    end
    
    cd ../..
end

% % sanity check plotting
% close all
% set(gcf, 'Position', [46 351 1872 450])
% plot(SDis)
% hold on
% line([0 endframe], [radius_cm radius_cm])
% scatter(cross_In, repmat(radius_cm, 1, length(cross_In)))
% scatter(cross_Out, repmat(radius_cm, 1, length(cross_Out)))
% 
% figure(2)
% hist(cross_Out-cross_In,100)

%%
numOrLen = cell(2*length(nameDir),1);
numOrLen(1:length(nameDir)) = {'boutNum'};
numOrLen(length(nameDir)+1:end) = {'boutLen'};
TableLabels = [cat(2,nameDir,nameDir)', numOrLen];
TableInfo   = [boutNum;boutLen];
Table       = [cell2table(TableLabels) array2table(TableInfo)];
%writetable(Table,'boutAnalysis_nose.csv');
%clear all

