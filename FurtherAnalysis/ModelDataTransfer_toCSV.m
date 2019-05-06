% Save MoSeq xPos/yPos info in csv file (requires MoSeqDataFrame.mat from ModelDataTransfer.py)
clear
close all
clc

Config_NovAna

currMouse = 'Nashville';
currDate  = '190425';

%% if loading from h5 file
cd(['/media/alex/DataDrive1/MoSeqData/Iku_photometry2/Iku_photometry2_MoSeq/' currMouse])

% load object and arena position info
cd(['./' currDate])
currSessionName = dir('session*');
cd(currSessionName.name)
cd('./proc/Analyzed_Data')
load('Arena_Obj_Pos.mat')
cd ..

offsetX = 80;
offsetY = 80;

xPos = h5read('results_00.h5', '/scalars/centroid_x_px')+offsetX;
yPos = h5read('results_00.h5', '/scalars/centroid_y_px')+offsetY;

xyPos = sqrt((obj_center(1,1)-xPos).^2+(obj_center(1,2)-yPos).^2);

close all
fig1 = figure(1);
hold on
plot(xPos(1:18000), yPos(1:18000), 'k.')
plot(obj_center(1), obj_center(2), 'r*')
plot(obj(1,[1 3]), obj(1,[2 4]), 'b*')
plot(obj(1,[1 3]), obj(1,[4 2]), 'b*')

th = 0:pi/50:2*pi;
x  = obj_center(1,1);
y  = obj_center(1,2);
xunit = radius * cos(th) + x;
yunit = radius * sin(th) + y;
plot(xunit, yunit)

axis square
set(gca, 'YDir', 'reverse')

disp('file output format:')
[xPos(1:10) yPos(1:10) xyPos(1:10)]
disp(['radius: ' num2str(radius)])

%% save
cd Analyzed_Data
csvwrite(['MoSeqPos_' currMouse '_' currDate], [xPos yPos xyPos])
saveas(fig1, ['MoSeqTrace_' currMouse '_' currDate 'r=' num2str(radius) '.tif'])
close all

%% if loading from MoSeqDataFrame
% load MoSeqDataFrame (contains raw xy position info)
load('MoSeqDataFrame.mat')

% load depth video (for sanity check)
cd(['./' currDate])
currSessionName = dir('session*');
cd(currSessionName.name)
cd('./proc')
vn = dir('*.mp4');
vn = vn.name;
video=VideoReader(vn);

% load object and arena position info
cd(['./Analyzed_Data'])
load('Arena_Obj_Pos.mat')

session_uuid_start = zeros(1,size(MoSeqDataFrame.session_uuid,1));
groupName = cell(1,size(MoSeqDataFrame.session_uuid,1));
for findSub = 1:size(MoSeqDataFrame.session_uuid,1)
    uuid_subset = ismember(MoSeqDataFrame.uuid, MoSeqDataFrame.session_uuid(findSub,:), 'rows');
    uuid_start  = find(uuid_subset,1,'first');
    session_uuid_start(1,findSub) = uuid_start;
    groupName{1,findSub} = MoSeqDataFrame.group(uuid_start,:);
end

whichGroup = find(strcmp(groupName, 'N1  ')); % 'Hab1'; 'Hab2'; 'N1  ';
sessionIndStart = session_uuid_start(whichGroup);
session_uuid_sort = sort([session_uuid_start length(MoSeqDataFrame.uuid)]);
sessionIndSort = find(sessionIndStart==session_uuid_sort);
sessionInd = sessionIndStart:session_uuid_sort(sessionIndSort+1)-1;

offsetX = 80;
offsetY = 80;

xPos = MoSeqDataFrame.centroid_x_px(sessionInd)+offsetX;
yPos = MoSeqDataFrame.centroid_y_px(sessionInd)+offsetY;

%% check to make sure xPos and yPos align with depth video and obj_center (from MarkObjPos.m)
close all
for i = 500:700
    disp(i)
    video.CurrentTime = (i-1)/video.FrameRate;
    frame = readFrame(video);
    
    imshow(frame)
    hold on
    plot(xPos(i), yPos(i), 'r*')
    pause(0.1)
end

%% navigating h5 file format
h5info('results_00.h5');
h5disp('results_00.h5');
h5read('results_00.h5', '/metadata/uuid');

xPos = h5read('results_00.h5', '/scalars/centroid_x_px');
yPos = h5read('results_00.h5', '/scalars/centroid_y_px');

