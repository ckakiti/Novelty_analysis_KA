% Analyze bouts of interaction with object

clear
close all
clc

Config_NovAna

% only analyze first 10 minutes of video
start_min=0.005;
end_min=start_min+10;
startframe=start_min.*60.*fps;
endframe=end_min.*60.*fps;
startframe=round(startframe);
endframe=round(endframe);

cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/NewHope-ROTJ/')
%cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/7day_preexposure_combine/')
files    = dir;
whichDir = [files.isdir];
nameDir  = files(whichDir);
nameDir  = {nameDir.name};
nameDir  = nameDir(3:end);
days     = 12;

% number of bouts and bout length
boutNum = zeros(length(nameDir),days);
boutLen = zeros(length(nameDir),days);

for mousei = 1:length(nameDir)
    cd([nameDir{1,mousei}, '/Analyzed_Data/body_r8'])
    
    matFiles = dir('*.mat');
    matFiles = cat(1, matFiles.name);
    
    for filei = 1:size(matFiles,1)
        load(matFiles(filei,:))
        
        % find frames where mouse crosses within radius specified by Config_NovAna
        crossTmp  = crossing(Labels(startframe:endframe,21), [], 0.5);
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
            boutLen(mousei,filei) = mean(cross_Out-cross_In);
        end
    end
    
    cd ../../..
end

numOrLen = cell(2*length(nameDir),1);
numOrLen(1:length(nameDir)) = {'boutNum'};
numOrLen(length(nameDir)+1:end) = {'boutLen'};
TableLabels = [cat(2,nameDir,nameDir)', numOrLen];
TableInfo   = [boutNum;boutLen];
Table       = [cell2table(TableLabels) array2table(TableInfo)];
%writetable(Table,'boutAnalysis_body.csv');
%clear all