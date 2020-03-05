%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

whichMouse   = 'Francisco';
whichSession = 'N1'; %N1; %Hab2

if(strcmp(whichSession, 'Hab2'))
    cd /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Francisco/190219/session_20190219153503
    MSid='619acd8f-e08e-4ef3-9ab1-c7ff0ce8b4d9';
else
    cd /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Francisco/190220/session_20190220120254
    MSid='c689b8b6-70c8-4a3b-850e-20cc3e7f8f76';
end

filename = 'depth_ts.txt';

delimiter = ' ';

% Format for each line of text:
%   column2: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*q%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, ...
    'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% Create output variable
depthts = [dataArray{1:end-1}];
% same length as number of frames in depth video (30fps) - NOT rgb video (15fps)

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('MoSeqFP.mat'); % created by FPDataTransfer.py (need to have fiber photometry data)
% loads fiber photometry data: ch00 = GCaMP, ch01=tdTom, tstep=time stamp
% of each frame of ch00/cho01 which corresponds to the numbers in depthts

cd ../..

load('MoSeqDataFrame.mat'); % created by ModelDataTransfer.py (need to have run MoSeq successfully)
% metadata for depth videos concatenated across all days (e.g. length of 3*54100=162300 for 3 days)

% load object center position (identified manually through MarkObjPos.m)
%load('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry_DLC/Francisco/Analyzed_Data/Arena_Obj_Pos.mat', 'obj_center')

ObjPos=[-200 -105];%[-200 -125]; %-obj_center(3,:); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalcount=length(MoSeqDataFrame.model_label);

FSize=20;

%MSid='7c39030c-3ec5-4366-ac15-513773fa5607'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MSid='c689b8b6-70c8-4a3b-850e-20cc3e7f8f76'; % Francisco 190220 novelty day 1 
%MSid='619acd8f-e08e-4ef3-9ab1-c7ff0ce8b4d9'; % Francisco 190219 habituation day 2

MoSeqDataFrame.SODis=sqrt((MoSeqDataFrame.centroid_x_mm-ObjPos(1)).^2+(MoSeqDataFrame.centroid_y_mm-ObjPos(2)).^2);

FilteredXpos=[];  % only for this MSid
FilteredYpos=[];  % only for this MSid
FilteredSODis=[]; % only for this MSid
FilteredVel=[];   % only for this MSid

for frameiter=1:length(MoSeqDataFrame.SODis)
    if strcmp(MSid, MoSeqDataFrame.uuid(frameiter,:))
        FilteredXpos=[FilteredXpos double(MoSeqDataFrame.centroid_x_mm(frameiter))];
        FilteredYpos=[FilteredYpos double(MoSeqDataFrame.centroid_y_mm(frameiter))];
        FilteredSODis=[FilteredSODis MoSeqDataFrame.SODis(frameiter)];
        FilteredVel=[FilteredVel MoSeqDataFrame.velocity_3d_mm(frameiter)];
    end
end

len=length(tstep);
deplen=length(FilteredSODis);%length(depthts);
Dis=zeros(1,len); % distance from object for each frame in FP data (convert depth frame to FP frame)
Vel=zeros(1,len);
depth_index=1;
for iter=1:len
    if depthts(depth_index)<tstep(iter)
        depth_index=depth_index+1;
    end
    if depth_index<=deplen
        Dis(iter)=FilteredSODis(depth_index);
        Vel(iter)=FilteredVel(depth_index);
    else
        disp(iter)
        disp(['break: ', num2str(iter)])
        break
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%acttime=csvread('NoveltyResponse.csv');
%acttime=acttime+1; % labeled frame number start at 0 while depthts start at 1

GCaMP=ch00;
tdTom=ch01;

%remove 60Hz noise
d = designfilt('bandstopiir','FilterOrder',2, ...
 'HalfPowerFrequency1',58,'HalfPowerFrequency2',62, ...
 'DesignMethod','butter','SampleRate',1000);
filtered_GCaMP = filtfilt(d,GCaMP);
GCaMP = filtered_GCaMP;
filtered_tdTom = filtfilt(d,tdTom);
tdTom = filtered_tdTom;

GCaMP=GCaMP-mean(GCaMP);
tdTom=tdTom-mean(tdTom);
GCaMP=GCaMP-GCaMP(1);
tdTom=tdTom-tdTom(1);

%lenbef=3000;
%lenaft=4000;

% figure;
% for actiter=1:size(acttime,1)

%     acttime_index=find(tstep>depthts(acttime(actiter)));
%     acttime_index=acttime_index(1);

%     Signal1=GCaMP(acttime_index-lenbef:acttime_index+lenaft);
%     Signal2=tdTom(acttime_index-lenbef:acttime_index+lenaft);
%     Signal3=Dis(acttime_index-lenbef:acttime_index+lenaft);
%     Signal4=Vel(acttime_index-lenbef:acttime_index+lenaft);
%     X=(0-lenbef):lenaft;

%     subplot(ceil(size(acttime,1)/2),2,actiter)
%     plot(X,Signal1,'Color','g')
%     hold on
%     plot(X,Signal2,'Color','r')
%     hold on
%     plot(X,Signal4/50,'Color','m')
%     xlabel('time (ms)')
%     ylabel('signals dF/F')
%     % legend('GCaMP','tdTom','Velocity')
% end


%% plotting
% find which rows of MoSeqDataFrame correspond to current MSid
whereMSid = find(strcmp(MSid, cellstr(MoSeqDataFrame.uuid)));

tspan_start = 0;
tspan_end   = 100000;
tspan_start = min(whereMSid+tspan_start);
tspan_end   = min(min(whereMSid+tspan_end),max(whereMSid));

Signal1=GCaMP(tspan_start:tspan_end);
Signal2=tdTom(tspan_start:tspan_end);
Signal3=Dis(tspan_start:tspan_end);
Signal4=Vel(tspan_start:tspan_end);

X=tstep-tstep(1);
X=X(tspan_start:tspan_end);


%% plot GCaMP, tdTom, velocity
figure(1);
plot(X,Signal1,'Color','g')
hold on
plot(X,Signal2,'Color','r')
hold on
plot(X,Signal4/50+0.3,'Color','m')
xlabel('time (s)')
ylabel('signals dF/F | Speed/50+0.3 mm/s')
legend('GCaMP','tdTom','Velocity')

%% plot GCaMP, tdTom, distance
figure(2);
plot(X,Signal1,'Color','g')
hold on
plot(X,Signal2,'Color','r')
hold on
plot(X,Signal3/1000+0.3,'Color','b','LineWidth',1.5)
xlabel('time (s)')
ylabel('signals dF/F | Distance to Obj (m)')
legend('GCaMP','tdTom','Distance')



%% Plot trajectory
Trafigure=figure(3);
scatter(FilteredXpos,FilteredYpos,'k.');
hold on
scatter(ObjPos(1),ObjPos(2),'MarkerEdgeColor','r')
% rectangle('Position',[arena(fiter,1),arena(fiter,2),arena(fiter,3)-arena(fiter,1),arena(fiter,4)-arena(fiter,2)],...
%           'EdgeColor','r','linewidth',4)
% rectangle('Position',[obj(fiter,1),obj(fiter,2),obj(fiter,3)-obj(fiter,1),obj(fiter,4)-obj(fiter,2)],...
%           'EdgeColor','r','linewidth',4)
%  set(gca,'ydir','reverse')
% xlim([0 video_xlen]);
% ylim([0 video_ywid]);
% set(Trafigure, 'position', [0 0 1200 900]);


%% plot GCaMP and tdTom vs velocity (and assess correlation if any)
close all

coeffs_GCaMP  = polyfit(Signal4, Signal1, 1);
coeffs_tdTom  = polyfit(Signal4, Signal2, 1);
fittedX       = linspace(min(Signal4), max(Signal4), 200);
fittedY_GCaMP = polyval(coeffs_GCaMP, fittedX);
fittedY_tdTom = polyval(coeffs_tdTom, fittedX);

fig4 = figure(4);
set(fig4, 'Position', [275 432 1059 450])
subplot(1,2,1)
hold on
scatter(Signal4, Signal1, 'g.') % GCaMP/velocity
plot(fittedX, fittedY_GCaMP, 'k-', 'LineWidth', 3);
title(['GCaMP vs velocity (30 min, ', whichMouse, ' ', whichSession, ')'])
xlabel('Velocity')
ylabel('GCaMP (mean subtracted)')
hold off

subplot(1,2,2)
hold on
scatter(Signal4, Signal2, 'r.') % tdTom/velocity
plot(fittedX, fittedY_tdTom, 'k-', 'LineWidth', 3);
title(['tdTom vs velocity (30 min, ', whichMouse, ' ', whichSession, ')'])
xlabel('Velocity')
ylabel('tdTom (mean subtracted)')
hold off

saveas(fig4, ['FPvsVelocity_' whichMouse '_' whichSession '.tif'])

%% calculating dF/F
GCaMP_F0 = median(GCaMP);
GCaMP_dFF = (GCaMP - GCaMP_F0)./GCaMP_F0;

figure(1)
hold on
plot(GCaMP(tspan_start:tspan_end), 'k')
plot(GCaMP_dFF(tspan_start:tspan_end), 'g')





