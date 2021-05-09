%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importing (1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

%Config_NovAna
radius_cm = 8;
ppc = 42/6.3;
fps = 15;
disp(['radius_cm: ' num2str(radius_cm), ' (unused)'])
disp(['ppc: ' num2str(ppc)])
disp(['fps: ' num2str(fps)])

whichMouse = 'no_mouse_MoSeq';
whichDate  = '210506';
whichSession = 'N1'; %Hab2 %N1 %R1(arm+delivery) %R2(consumption)
whichFile = 1;

% cd(['/media/alex/DataDrive1/MoSeqData/Iku_photometry2/Iku_photometry2_MoSeq/' ...
%    whichMouse '/' whichDate])
cd(['~/Dropbox (Uchida Lab)/Korleki Akiti/Behavior/FP_test_210506/' ...
    whichMouse '/' whichDate])

filelist = dir('session*');
cd(filelist(whichFile).name)

filename = 'depth_ts.txt';

delimiter = ' ';

% Format for each line of text:
%   column2: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%[^\n\r]';
%formatSpec = '%*q%f%[^\n\r]';

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
depthts = [dataArray{1}]; %%%%%%%%%%%%%%%%%%%%%%%%%

% get rgb_ts array as well
filename = 'rgb_ts.txt';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, ...
    'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
rgbts = [dataArray{1}];

% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

disp('section 1')


%% alternative section 1+2 (for Mitsuko_photometry, LED detection)
clear
pause(0.01)

close all
clc

Config_NovAna
disp(['radius_cm: ' num2str(radius_cm), ' (unused)'])
disp(['ppc: ' num2str(ppc)])
disp(['fps: ' num2str(fps)])

whichMouse = 'Nashville';
whichDate  = '190425';
whichSession = 'N1'; %Hab2 %N1 %R1(arm+delivery) %R2(consumption)
whichFile = 1;

cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry2_DLC/' ...
    whichMouse '/' whichDate])

disp('section 1')


FP_dir = dir('*_FP.mat');
load(FP_dir.name)
ch00=cart_GCaMP;
ch01=cart_tdTom;

rgbts_file = dir('*rgb_ts');
rgbts = csvread(rgbts_file.name);

poke_file = dir('*poke');
approach_file = dir('*boutStart');
acttime = csvread(poke_file.name);
acttime2 = csvread(approach_file.name);

% convert rgb/depth frame number to FP frame number
actlen=length(acttime); %size(acttime,1);
acttime_index=zeros(actlen,1);
acttime2_index=zeros(actlen,1);
for actiter=1:actlen
    acttime_curr = find(rgbts(:,1)==acttime(actiter));
    acttime2_curr = find(rgbts(:,1)==acttime2(actiter));
    
    if(~isempty(acttime_curr) && ~isempty(acttime2_curr))
        acttime_index(actiter,1) = rgbts(acttime_curr,2);
        acttime2_index(actiter,1) = rgbts(acttime2_curr,2);
    end
end

% acttime_index = pos_within_TTL;
% acttime2_index = acttime_index;

disp('section 2')

%% (2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('MoSeqFP.mat');

if(strcmp(whichSession,'N1'))
%     boutFile = dir('*boutStart'); %nose
%     pokeFile = dir('*poke'); %nose
    boutFile = dir('poke*');
    pokeFile = dir('poke*');
elseif(strcmp(whichSession,'R1'))
    boutFile = dir('RewardArm*');
    pokeFile = dir('RewardDelivery*');
elseif(strcmp(whichSession,'R2'))
    boutFile = dir('RewardResponse*');
    pokeFile = boutFile;
end

acttime=csvread(boutFile.name);
acttime2=csvread(pokeFile.name);
%LEDtime=csvread('LED_on');

%acttime=acttime+1; % labeled frame number start at 0 while depthts start at 1
%acttime=acttime+900; % when extracting, frames were trimmed by 900

%% convert poke frame number (in depthts) to FP frame number
rgb_FP_file = dir('*_rgb_ts');
rgb_FP = load(rgb_FP_file(1).name);

poke_depth_rgb    = acttime;
poke_depth_rgb_FP = zeros(length(acttime),1);
for ts_iter = 1:length(acttime)
    poke_depth_rgb(ts_iter,2) = find(rgbts>depthts(acttime(ts_iter),1),1);
    rgb_FP_idx = find(poke_depth_rgb(ts_iter,2)==rgb_FP(:,1));
    poke_depth_rgb_FP(ts_iter,1) = rgb_FP(rgb_FP_idx,2);
end
poke_depth_rgb_FP = [poke_depth_rgb poke_depth_rgb_FP];

acttime_index  = poke_depth_rgb_FP(:,3);
acttime2_index = poke_depth_rgb_FP(:,3);

%% convert rgb/depth frame number to FP frame number
actlen=length(acttime); %size(acttime,1);
acttime_index=zeros(actlen,1);
acttime2_index=zeros(actlen,1);

orienttime=zeros(actlen,1); % first orientation after start of bout
orienttime_index=zeros(actlen,1);

for actiter=1:actlen
    acttime_index(actiter,1)=find(tstep>rgbts(acttime(actiter)),1); %%%%%%%%%%%%%%%
    acttime2_index(actiter,1)=find(tstep>rgbts(acttime2(actiter)),1); %%%%%%%%%%%%%%%
    
%     oriented = acttime(actiter) + find(DLCPos.Labels(acttime(actiter):end,23),1);
%     orienttime(actiter,1) = oriented-1;
%     orienttime_index(actiter,1)=find(tstep>rgbts(orienttime(actiter)),1); %%%%%%%%%%%%%%%
end


% allOrient = crossing(DLCPos.Labels(:,23)-0.5);
% allOrientStart = allOrient(1:2:end); % find every time mice orient to object
% allorient_index=zeros(length(allOrientStart),1);
% for orientiter = 1:length(allOrientStart)
%     allorient_index(orientiter,1)=find(tstep>rgbts(allOrientStart(orientiter)),1);
% end
% allorient_index = allorient_index(allorient_index<796200);

%acttime_index=allorient_index;
%acttime2_index=allorient_index;

% convert LED_on frame number to FP frame number
% LEDlen=length(LEDtime); %size(acttime,1);
% LED_index=zeros(LEDlen,1);
% for LEDiter=1:LEDlen
%     LED_index(LEDiter,1)=find(tstep>rgbts(LEDtime(LEDiter)),1,'first');
% end

disp('section 2')

%% calculate velocity and acceleration

%cd ./proc/Analyzed_Data
%posFile = dir('MoSeqPos*');
%MoSeqPos = load(posFile.name);
%vel = MoSeqPos(:,3);
%cd ../..
cd Analyzed_Data_1obj %%%%%%%
posFile = dir('*Converted.mat');% dir('*0000.mat'); % dir('*Converted.mat');
DLCPos = load(posFile.name);
load('Arena_Obj_Pos.mat');
%obj_center = [-206.9 -129.2]; % for MoSeqPos Nashville 190425R centroid_mm
%obj_center = [0 0];           % for calculating vel/acc relative to 0

cd ..

if(isempty(dir('*rgb_ts')))
    tic
    cutoff = 9000; %60 %9000 %18000
    veltime_ind = zeros(cutoff,1);
    for veliter = 1:length(veltime_ind)
        disp(veliter)
        veltime_ind(veliter,1) = find(tstep>rgbts(veliter),1,'first'); %%%%%%%%%%%%%%%%
        %[d, ix] = min(abs(tstep-depthts(veliter)));
        %veltime_ind(veliter-900,1) = ix;
    end
    toc
else
    rgbts_file = dir('*rgb_ts');
    veltime_file = load(rgbts_file.name);
    veltime_crop = veltime_file(:,1);
    veltime_ind = veltime_file(:,2);
    
    cutoff = length(veltime_crop);
end

% smooth x and y positions
xPos_body = DLCPos.Labels(veltime_crop,2);%14);
yPos_body = DLCPos.Labels(veltime_crop,3);%15);
xPos_smooth = smooth(xPos_body,30);
yPos_smooth = smooth(yPos_body,30);

% first average x and y position every 15 frames (1s)
%  then calculate velocity (difference in distance / dt)
posWin = fps;
disp(['fps: ' num2str(fps)])
modPos = mod(length(xPos_smooth),posWin);
xPos_smooth_cut = xPos_smooth(1:end-modPos);
yPos_smooth_cut = yPos_smooth(1:end-modPos);

reshape_xPos = reshape(xPos_smooth_cut, posWin, []);
reshape_yPos = reshape(yPos_smooth_cut, posWin, []);
avg_xPos = mean(reshape_xPos,1);
avg_yPos = mean(reshape_yPos,1);

% calculating velocity+acceleration from averaged xPos/yPos
diff_xPos = diff(avg_xPos);
diff_yPos = diff(avg_yPos);
dt = posWin/fps;
vel = sqrt( (diff_xPos .* diff_xPos) + (diff_yPos .* diff_yPos) )/(dt*ppc); %units = cm/s

vel_smooth = smooth(vel, 4);
acc = diff(vel_smooth); %units = cm/s^2
% acc_smooth = smooth(acc, 4);

disp('calculated vel/acc')

%% (3) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GCaMP=ch00;
tdTom=ch01;

GCaMP_unproc = GCaMP; %%%%

% filter 60Hz noise from photometry signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1000 %700:10:1000
%remove 60Hz noise
d = designfilt('bandstopiir',...
               'FilterOrder',2, ...
               'HalfPowerFrequency1',58,... %58
               'HalfPowerFrequency2',62, ... %62
               'DesignMethod','butter',...
               'SampleRate',i); %1000
GCaMP_filtered = filtfilt(d,GCaMP);

% close all
% x = 1:700000;
% figure(1)
% hold on
% plot(x, GCaMP_unproc(x), 'k-')
% plot(x, GCaMP_filtered(x)+0.00, 'r-')
% xlim([471200 471300])
% legend({'Unprocessed', '60Hz filtered'})
% 
% disp(i)
% pause

end

GCaMP = GCaMP_filtered;
filtered_tdTom = filtfilt(d,tdTom);
tdTom = filtered_tdTom;

disp('section 3')

%% remove sudden rise nose then smooth (4)

%remove sudden rise noise
diff_GCaMP = diff(GCaMP);
std_GCaMP = std(diff_GCaMP);
ind_noise = find(diff_GCaMP>5*std_GCaMP);
exclude_green = size(ind_noise);

for i = 1:length(ind_noise)
    GCaMP(ind_noise(i))=GCaMP(ind_noise(i)-1);
    for k=1:1000
        if GCaMP(ind_noise(i)+k)-GCaMP(ind_noise(i)+k-1)>5*std_GCaMP
            GCaMP(ind_noise(i)+k) = GCaMP(ind_noise(i)+k-1);
        end
    end
end

diff_tdTom = diff(tdTom);
std_tdTom = std(diff_tdTom);
ind_noise = find(diff_tdTom>5*std_tdTom);
exclude_red = size(ind_noise);

for i = 1:length(ind_noise)
    tdTom(ind_noise(i))=tdTom(ind_noise(i)-1);
    for k=1:1000
        if tdTom(ind_noise(i)+k)-tdTom(ind_noise(i)+k-1)>5*std_tdTom
            tdTom(ind_noise(i)+k) = tdTom(ind_noise(i)+k-1);
        end
    end
end

GCaMP_noRise = GCaMP; %%%%

%smoothing
normG = smooth(GCaMP,50);
GCaMP = normG;
normR = smooth(tdTom,50);
tdTom = normR;

GCaMP_smooth = GCaMP; %%%%


%GCaMP_unproc;
%GCaMP_filtered;
%GCaMP_noRise;
%GCaMP_smooth;

% x = 1:700000;
% figure(1)
% hold on
% plot(x, GCaMP_unproc(x), 'k-')
% plot(x, GCaMP_filtered(x)+0.05, 'r-')
% plot(x, GCaMP_noRise(x)+0.1, 'g-')
% plot(x, GCaMP_smooth(x)+0.15, 'b-')
% legend({'Unprocessed', '60Hz filtered', 'Correct sudden rise', 'Smoothed'})

disp('section 4')

%% (5) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trace average

plotWin = -5000:5000;

plotind = repmat(plotWin,length(acttime_index),1)+acttime_index;
%plotind = repmat(plotWin,length(allorient_index),1)+allorient_index;%orienttime_index;
plotind(plotind<=0)=1;
rawTrace = GCaMP(plotind);
R_rawTrace = tdTom(plotind);

plotind2 = repmat(plotWin,length(acttime2_index),1)+acttime2_index;
plotind2(plotind2<=0)=1;
rawTrace2 = GCaMP(plotind2);
R_rawTrace2 = tdTom(plotind2);


% deltaF/F = [(rawTrace - Fsubtr) / Fdiv] * 100
%   rawTrace: rows = trials, columns = FP signal centered around event
%   Fsubtr: baseline signal for each trial - avg 1s from beginning of trial (diff for each trial)
%   Fdiv: mean of first sec of all trials

baseShift = 0; %0=regular; 1000=shift 1sec closer to boutStart

Fdiv    = mean(mean(rawTrace(:,(1:1000)+baseShift)));
Fsubtr  = mean(rawTrace(:,(1001:2000)+baseShift),2); %%%%%%%%      % using this time window in plotwin for baseline
deltaF  = ((rawTrace-Fsubtr)/Fdiv)*100;
deltaF2 = ((rawTrace2-Fsubtr)/Fdiv)*100;

R_Fdiv    = mean(mean(R_rawTrace(:,(1:1000)+baseShift)));
R_Fsubtr  = mean(R_rawTrace(:,(1001:2000)+baseShift),2);    % using this time window in plotwin for baseline
R_deltaF  = ((R_rawTrace-R_Fsubtr)/R_Fdiv)*100;
R_deltaF2 = ((R_rawTrace2-R_Fsubtr)/R_Fdiv)*100;

% %% correct Gardner GCaMP
% figure(1);
% hold on
% for i = 1:length(R_rawTrace)
%     plot(plotWin, R_rawTrace(i,:))
%     disp(i)
%     pause
% end
% 
% rawTrace(27,:) = [];
% R_rawTrace(27,:) = [];

disp(['F sub/div shift: ' num2str(baseShift) 'ms'])
disp('section 5')

%% plot (6)
m_plot = mean(deltaF);
s_plot = std(deltaF)/sqrt(length(acttime_index));

R_m_plot = mean(R_deltaF);
R_s_plot = std(R_deltaF)/sqrt(length(acttime_index));

m_plot2 = mean(deltaF2);
s_plot2 = std(deltaF2)/sqrt(length(acttime2_index));

R_m_plot2 = mean(R_deltaF2);
R_s_plot2 = std(R_deltaF2)/sqrt(length(acttime2_index));

disp('section 6')
%%
close all

if(strcmp(whichSession,'N1') || strcmp(whichSession,'N2') || strcmp(whichSession,'H2'))
    title_dFF_1 = 'GCaMP signals during approach';
    xlabel_dFF_1 = 'time - point of bout start (ms)';
    
    title_dFF_2 = 'GCaMP signals during poke'; %retreat'
    xlabel_dFF_2 = 'time - point of poke (ms)'; %retreat (ms)'
    
    xlabel_raster_1 = 'time - approach (s)';
    xlabel_raster_2 = 'time - point of poke (s)'; %retreat (s)';
     
%     title_dFF_1 = 'GCaMP signals during orientation';
%     xlabel_dFF_1 = 'time - point of orientation (ms)';
%     
%     title_dFF_2 = 'GCaMP signals during orientation'; %retreat'
%     xlabel_dFF_2 = 'time - point of orientation (ms)'; %retreat (ms)'
%     
%     xlabel_raster_1 = 'time - point of orientation (s)';
%     xlabel_raster_2 = 'time - point of orientation (s)'; %retreat (s)';
    
elseif(strcmp(whichSession,'R1'))
    title_dFF_1 = 'GCaMP signals to arm reaching in';
    xlabel_dFF_1 = 'time - point of arm reaching in (ms)';
    
    title_dFF_2 = 'GCaMP signals to reward delivery';
    xlabel_dFF_2 = 'time - point of reward delivery (ms)';
    
    xlabel_raster_1 = 'time - point of arm reaching in (s)';
    xlabel_raster_2 = 'time - point of reward delivery (s)';
    
elseif(strcmp(whichSession,'R2'))
    title_dFF_1 = 'GCaMP signals to reward consumption';
    xlabel_dFF_1 = 'time - point of consumption (ms)';
    
    title_dFF_2 = 'GCaMP signals to reward consumption';
    xlabel_dFF_2 = 'time - point of consumption (ms)';
    
    xlabel_raster_1 = 'time - reward consumption (s)';
    xlabel_raster_2 = 'time - reward consumption (s)';
end

% average GCaMP activity aligned to event (poke, bout start, reward delivery, etc.)
fig1 = figure(1);
suptitle([whichMouse, '-' whichDate])
set(gcf, 'Position', [217 511 1089 518])

fig_1 = subplot(1,2,1);
errorbar_patch(plotWin,m_plot,s_plot,[0 1 0]);
errorbar_patch(plotWin,R_m_plot,R_s_plot,[1 0 0]);
title(title_dFF_1)
legend('GCaMP','tdTom')
xlabel(xlabel_dFF_1)
ylabel('dF/F (%)')
ylim([-4 4])

fig_2 = subplot(1,2,2);
errorbar_patch(plotWin,m_plot2,s_plot2,[0 1 0]);
errorbar_patch(plotWin,R_m_plot2,R_s_plot2,[1 0 0]);
title(title_dFF_2)
legend('GCaMP','tdTom')
xlabel(xlabel_dFF_2)
ylabel('dF/F (%)')
fig_2.YLim = fig_1.YLim;


% GCaMP raster plot (10ms window average)
rasterWin = 10;
a1 = reshape(deltaF(:,1:(length(plotWin)-1)), length(acttime_index), rasterWin, []);
b1 = squeeze(mean(a1,2));
a2 = reshape(deltaF2(:,1:(length(plotWin)-1)), length(acttime_index), rasterWin, []);
b2 = squeeze(mean(a2,2));

a3 = reshape(R_deltaF(:,1:(length(plotWin)-1)), length(acttime_index), rasterWin, []); %tdTom
b3 = squeeze(mean(a3,2));
a4 = reshape(R_deltaF2(:,1:(length(plotWin)-1)), length(acttime_index), rasterWin, []); %tdTom
b4 = squeeze(mean(a4,2));

c = plotWin(1:1000:end)/1000;


fig2 = figure(2);
%suptitle(whichMouse)
set(gcf, 'Position', [200 300 1400 500])

fig2_1 = subplot(1,2,1);
imagesc(b1, [-4 4])
title(['GCaMP raster (' num2str(rasterWin) 'ms bin): ' whichMouse])

fig2_1.XTick = 0:100:1000;
fig2_1.XTickLabel = c;
fig2_1.XLabel.String = xlabel_raster_1;
fig2_1.YLabel.String = 'trial';
colorbar
axis square

fig2_2 = subplot(1,2,2);
imagesc(b2, [-4 4])
title(['GCaMP raster (' num2str(rasterWin) 'ms bin): ' whichMouse])
fig2_2.XTick = 0:100:1000;
fig2_2.XTickLabel = c;
fig2_2.XLabel.String = xlabel_raster_2;
fig2_2.YLabel.String = 'trial';
colorbar
axis square


% tdTom raster plot (10ms window average)
fig3 = figure(3);
%suptitle(whichMouse)
set(gcf, 'Position', [200 100 1400 500])

fig3_1 = subplot(1,2,1);
imagesc(b3, [-4 4])
title(['tdTom raster (' num2str(rasterWin) 'ms bin): ' whichMouse])
fig3_1.XTick = 0:100:1000;
fig3_1.XTickLabel = c;
fig3_1.XLabel.String = xlabel_raster_1;
fig3_1.YLabel.String = 'trial';
colorbar
axis square

fig3_2 = subplot(1,2,2);
imagesc(b4, [-4 4])
title(['tdTom raster (' num2str(rasterWin) 'ms bin): ' whichMouse])
fig3_2.XTick = 0:100:1000;
fig3_2.XTickLabel = c;
fig3_2.XLabel.String = xlabel_raster_2;
fig3_2.YLabel.String = 'trial';
colorbar
axis square

% save
wantToSave = input('Save? 0/1: ');
if(wantToSave)
%cd ./proc/Analyzed_Data
cd ./Analyzed_Data_1obj
saveas(fig1, ['FP_' whichMouse '_' whichSession '_DLC(r' num2str(radius_cm) ')_baseShift:'...
    num2str(baseShift) '.tif'])
saveas(fig2, ['FP_raster_GCaMP_' whichMouse '_' whichSession '_poke_DLC_baseShift:'...
    num2str(baseShift) '.tif'])
saveas(fig3, ['FP_raster_tdTom_' whichMouse '_' whichSession '_poke_DLC_baseShift:'...
    num2str(baseShift) '.tif'])
close all
end

%% plot GCaMP activity as a function of position, velocity, and acceleration
GCaMP_curr = GCaMP(veltime_ind);
GCaMP_ind  = 1:posWin:cutoff;

% round velocities to fall into 1 of 6 (or more) discrete bins
bin = (max(GCaMP_curr)-min(GCaMP_curr))/50; 
roundTargets = min(GCaMP_curr):bin:max(GCaMP_curr);
dataPosRound = interp1(roundTargets, roundTargets, GCaMP_curr, 'nearest');
z0 = zeros(size(veltime_ind));
col = dataPosRound';

close all
intervalXfps = 1:cutoff/2;
fig3 = figure(3);
set(gcf, 'Position', [45 366 1873 505])
suptitle([whichMouse ': ' whichSession])
subplot(1,3,1)
surface([xPos_smooth(intervalXfps),xPos_smooth(intervalXfps)]/ppc', ...
        [yPos_smooth(intervalXfps),yPos_smooth(intervalXfps)]/ppc', ...
        [z0(intervalXfps),z0(intervalXfps)], [col(intervalXfps);col(intervalXfps)]',...
    'marker', '.','markersize',2,'markerfacecol',[0 0 0],'edgecol','interp');
axis square
colorbar
colormap cool
%caxis([-4 4])

subplot(1,3,2)
%scatter(vel_smooth(1:(cutoff/posWin)), GCaMP_curr(1:posWin:cutoff), 'b.')
GCaMP_curr_vel = GCaMP_curr(GCaMP_ind(1:end-2));
scatter(vel_smooth(GCaMP_curr_vel>0), GCaMP_curr_vel(GCaMP_curr_vel>0), 'b.')
xlabel('Velocity (cm/s)')
ylabel('GCaMP (raw)')
%ylim([0.9 1.3])

subplot(1,3,3)
GCaMP_curr_acc = GCaMP_curr(GCaMP_ind(1:end-3));
scatter(acc(GCaMP_curr_acc>0), GCaMP_curr_acc(GCaMP_curr_acc>0), 'r.')
xlabel('Acceleration (cm/s^2)')
ylabel('GCaMP (raw)')
%ylim([0.9 1.3])

if(0)
cd('./proc/Analyzed_Data')
saveas(fig3, ['FP_posVelAccGCaMP_' whichMouse '_' whichSession '_DLC.tif'])
cd ../..
end

%% plot x, y position, velocity, and acceleration over first 10 min
intervalXfps = 1:cutoff;
timeMin = (1:cutoff)/fps/60;

close all
fig5 = figure(5);
set(gcf, 'Position', [45 350 1873 505])
title([whichMouse ': ' whichSession])
hold on
plot(timeMin, xPos_smooth(intervalXfps)/ppc')
plot(timeMin, yPos_smooth(intervalXfps)/ppc')
xlabel('Time (min)')
ylabel('Position (cm); smoothed, relative to (0,0)')
legend({'x-position', 'y-position'})

fig6 = figure(6);
set(gcf, 'Position', [45 250 1873 505])
title([whichMouse ': ' whichSession])
hold on
plot((1:cutoff/fps)/60, vel_smooth(1:cutoff/fps))
xlabel('Time (min)')
ylabel('Velocity (cm/s); smoothed, relative to object')

fig7 = figure(7);
set(gcf, 'Position', [45 150 1873 505])
title([whichMouse ': ' whichSession])
hold on
plot((1:cutoff/fps)/60, acc(1:cutoff/fps))
xlabel('Time (min)')
ylabel('Acceleration (cm/s^2); smoothed, relative to object')


if(0)
cd('./proc/Analyzed_Data')
saveas(fig5, ['MoSeqFP_posXY_' whichMouse '_' whichSession '_DLC.tif'])
saveas(fig6, ['MoSeqFP_vel_' whichMouse '_' whichSession '_DLC.tif'])
saveas(fig7, ['MoSeqFP_acc_' whichMouse '_' whichSession '_DLC.tif'])
cd ../..
end

%% plot bout trajectories and color code nose position based on GCaMP activity
xPos_nose = xPos_smooth;%MoSeqPos(:,1); %DLCPos.Labels(:,2)/ppc;
yPos_nose = yPos_smooth;%MoSeqPos(:,2); %DLCPos.Labels(:,3)/ppc;
fig4 = figure(4);
hold on
for iter = 1:length(acttime)
    currInt = acttime(iter):(acttime2(iter)+100);
    
    surface([xPos_nose(currInt),xPos_nose(currInt)]', ...
        [yPos_nose(currInt),yPos_nose(currInt)]', ...
        [z0(currInt),z0(currInt)]', [col(currInt);col(currInt)],...
        'marker', '.','markersize',2,'markerfacecol',[0 0 0],'edgecol','interp');
end
%plot(obj_center(1)/ppc, obj_center(2)/ppc, 'r*')
plot(obj_center(1), obj_center(2), 'r*')
set(gca, 'YDir', 'reverse')
axis square
colorbar

%% save (7)
cd ./proc/Analyzed_Data
saveas(fig1, ['MoSeqFP_' whichMouse '_' whichSession '_DLC(r' num2str(radius_cm) ').tif'])
saveas(fig2, ['MoSeqFP_raster_GCaMP_' whichMouse '_poke_DLC.tif'])
saveas(fig3, ['MoSeqFP_raster_tdTom_' whichMouse '_poke_DLC.tif'])
%saveas(fig3, ['MoSeqFP_posVelAccGCaMP_' whichMouse '_' whichSession '_DLC.tif'])
%saveas(fig4, ['MoSeqFP_boutTraceNoseGCaMP_' whichMouse '_' whichSession '_DLC.tif'])


%saveas(fig1, ['MoSeqFP_' whichMouse '_' whichSession '_MoSeq.tif'])
%saveas(fig1, ['MoSeqFP_' whichMouse '_ArmReachRewDelivery_MoSeq.tif'])
%saveas(fig2, ['MoSeqFP_raster_' whichMouse '_' whichSession '_MoSeq.tif'])
%saveas(fig2, ['MoSeqFP_raster_' whichMouse '_ArmReachRewDelivery_MoSeq.tif'])%RewardDelivery_MoSeq.tif'])
%saveas(fig3, ['MoSeqFP_posVelAccGCaMP_' whichMouse '_' whichSession '_MoSeq.tif'])
close all

%% sanity check for pokes (not part of rest of analysis)
close all
% manually identified light on for IR_Rse
%on_rgb    = [197 258 417 419 457];
%on_depth  = [234 306 488 491 537];
% p = polyfit(on_rgb, on_depth, 1);

%on_depth = load('NoveltyResponse_Nashville_Day3_boutStart'); %poke %boutStart
on_depth = csvread('RewardResponse_Nashville_190425R.csv');

% on_rgb_FP = zeros(1, length(on_rgb));
on_depth_FP = zeros(1, length(on_depth));
for iter = 1:length(on_depth)
%     [d, ix] = min(abs(rgbts(on_rgb(iter))-tstep));
%     on_rgb_FP(1,iter) = ix;
    [d2, ix2] = min(abs(depthts(on_depth(iter))-tstep));
    on_depth_FP(1,iter) = ix2;
end

% figure(1)
% hold on
% plot(on_rgb, on_depth, 'ko-')
% line([0 500], [p(2) p(1)*500+p(2)], 'color', 'r')
% xlim([0 550])
% ylim([0 550])
% axis square

fig2 = figure(2);
hold on
ch00_plot = plot(ch00, 'g');
ch01_plot = plot(ch01, 'r');
ch02_plot = plot(ch02, 'k');
depthon_line = line([on_depth_FP; on_depth_FP], repmat([-5 10]', 1, length(on_depth_FP)), 'color', 'b');
%rgbon_line = line([on_rgb_FP; on_rgb_FP], repmat([-5 10]', 1, length(on_rgb_FP)), 'color', 'c');
legend([ch00_plot ch01_plot ch02_plot depthon_line(1)], ... rgbon_line(1)], ...
    {'GCaMP', 'tdTom', 'LED', 'depth on'})%, 'rgb on'})
%xlim([1555000 1578000])
ylim([-6 13])

% oneOn numbers
%    event  rgb_frame  depth_frame  FP_frame  rgb_ts             depth_ts           tstep
%    on     434        510          18620     63691196802.64810  63691196805.48410  63691196805.47710
%    off    630        742          26750     63691196809.31710  63691196813.58610  63691196813.60710

%% sample 60Hz filter for generated sine wave
close all

Fs = 1000;                   % samples per second
dt = 1/Fs;                   % seconds per sample
StopTime = 1;             % seconds
t = (0:dt:StopTime-dt)';     % seconds
%%Sine wave:
Fc = 60;                     % hertz
x1 = cos(2*pi*Fc*t);
x2 = cos(2*pi*Fc/10*t);
x3 = x1+x2;

d = designfilt('bandstopiir',...
               'FilterOrder',2, ...
               'HalfPowerFrequency1',58,... %58
               'HalfPowerFrequency2',62, ... %62
               'DesignMethod','butter',...
               'SampleRate',1000);
filtered_x3 = filtfilt(d,x3);

% Plot the signal versus time:
figure;
hold on
plot(t,x3,'k-');
plot(t,filtered_x3,'r-');
xlabel('time (in seconds)');
title('Signal versus Time');
legend({'Unprocessed', '60Hz filtered'})
%zoom xon;

%%
if(0) 
% Response average trace for conparison
response = deltaF(:,4000:5000); %4000 is trigger
response = mean(response,2);
response = response';
ste_response = std(response)/sqrt(length(response));

% randon trace for conparison
randtime_index=round((rand(actlen,1).*0.8+0.1).*length(tstep),0);
randind = repmat(plotWin,length(randtime_index),1)+randtime_index;
rand_rawTrace = GCaMP(randind);
rand_F = mean(rand_rawTrace(:,1:1000),2);        %using this time window in plotwin for baseline
rand_deltaF = rand_rawTrace-rand_F;

rand_response = rand_deltaF(:,4000:5000); %4000 is trigger
rand_response = mean(rand_response,2);

rand_response = rand_response';
ste_rand_response = std(rand_response)/sqrt(length(rand_response));


X=1:2;
XTick={'Retreat Response' 'Random'};

figure
Yres=[mean(response) mean(rand_response)];
Yerr=[ste_response ste_rand_response];
bar(X,Yres)
hold on
errorbar(X,Yres,Yerr,'b.','LineWidth',1.5)
hold on
scatter(ones(1,actlen),response,'filled')
hold on
scatter(2.*ones(1,actlen),rand_response,'filled')
legend('GCaMP response','Standard Error');
title('GCaMP response')
ylabel('signals dF/F')
xticks(X);
xticklabels(XTick);


% raster plot green
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/1.5 scrsz(3)/1.5 scrsz(4)/1.5])
%1) bin the data
trialNum = size(deltaF,1); 
binSize = 100;
length_x = plotWin(end)-plotWin(1);
[sorted_response,sr_index]=sort(response','descend');
sorted_deltaF=deltaF(sr_index,:);
binedF = squeeze(mean(reshape(sorted_deltaF(:,1:length_x),trialNum, binSize,[]),2));
% imagesc(binedF,[-1 1]);
% imagesc(binedF,[-0.2 0.2]);
% imagesc(binedF,[-0.05 0.05]);
imagesc(binedF,[-0.15 0.15]);
colormap parula
colorbar
h=gca;
h.XTick = [0:10:(length_x/binSize)];
h.XTickLabel = {(plotWin(1)/1000):(plotWin(end)/1000)};
title('GCaMP signals during approach-retreat')
xlabel('time - point of retreat (s)');
hold on

% 2) plot the triggers
Xtri = [-plotWin(1)/binSize -plotWin(1)/binSize];
plot(Xtri,[0 trialNum+0.5],'r','LineWidth',2) 
end


function h = errorbar_patch(x,y,er,c)

    %ERRORBAR_PATCH    - errorbar by patch
    %
    % errorbar_patch(x,y,er,c)
    %
    %   input
    %     - x:    
    %     - y:    mean
    %     - er:    
    %     - color
    %     - N:          start y value for plot
    %     
    %   output
    %
    if nargin < 4
        c = [0 0 1];
    end
    
    if size(x,1)~= 1
        x = x';
    end
    
    if size(y,1)~= 1
        y = y';
    end
    
    if size(x,1)~= 1
        er = er';
    end
    
    X = [x fliplr(x)];
    Y = [y+er fliplr(y-er)];
    h1 = patch(X,Y,c,'edgecolor','none','FaceAlpha',0.2); hold on
    h2 = plot(x,y,'color',c,'LineWidth',1.5);
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    if nargout>0
        h = [h1 h2]; 
    end
    
end    