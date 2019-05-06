%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importing (1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

whichMouse = 'Vegas';
whichSession = 'R'; %Hab2 %N1 %R
whichFile = 1;

cd /media/alex/DataDrive1/MoSeqData/Iku_photometry2/Iku_photometry2_MoSeq/Vegas/190428R
%cd /media/alex/DataDrive1/MoSeqData/MSFP_Test/DataSet180922/session_20180922154525
%cd /media/alex/DataDrive1/MoSeqData/Iku_photometry/Iku_photometry_MoSeq/Francisco/190219/session_20190219153503

filelist = dir('session*');
cd(filelist(whichFile).name)

filename = 'depth_ts.txt';

delimiter = ' ';

% Format for each line of text:
%   column2: double (%f)
% For more information, see the TEXTSCAN documentation.
%formatSpec = '%f%[^\n\r]';
formatSpec = '%*q%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);

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
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
rgbts = [dataArray{1}];


% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

disp('section 1')

%% (2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('MoSeqFP.mat');

%boutFile = dir('RewardResponse*');
%pokeFile = boutFile;
%boutFile = dir('*boutStart');
%pokeFile = dir('*poke');
boutFile = dir('RewardArm*');
pokeFile = dir('RewardDelivery*');

acttime=csvread(boutFile.name);
acttime2=csvread(pokeFile.name);

%acttime=load(['NoveltyResponse_' whichMouse '_Day3_boutStart']); %poke %boutStart
%acttime=acttime+1; % labeled frame number start at 0 while depthts start at 1
%acttime=acttime+900; % when extracting, frames were trimmed by 900
%rescaling = length(depthts)/length(rgbts);
%acttime = round(acttime*rescaling);

% convert depth frame number to FP frame number
actlen=length(acttime); %size(acttime,1);
acttime_index=zeros(actlen,1);
acttime2_index=zeros(actlen,1);
for actiter=1:actlen
    acttime_index(actiter,1)=find(tstep>depthts(acttime(actiter)),1,'first');
    acttime2_index(actiter,1)=find(tstep>depthts(acttime2(actiter)),1,'first');
end

% tic
% veltime_ind = zeros(length(depthts)-900,1);
% for veliter = 901:length(depthts)
%     disp(veliter)
%     veltime_ind(veliter-900,1) = find(tstep>depthts(veliter),1,'first');
%     %[d, ix] = min(abs(tstep-depthts(veliter)));
%     %veltime_ind(veliter-900,1) = ix;
% end
% toc

cd ./proc/Analyzed_Data
posFile = dir('MoSeqPos*');
MoSeqPos = load(posFile.name);
vel = MoSeqPos(:,3);
cd ../..

disp('section 2')

%% (3) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

Fdiv    = mean(mean(rawTrace(:,1:1000)));
Fsubtr  = mean(rawTrace(:,1001:2000),2);        % using this time window in plotwin for baseline
deltaF  = ((rawTrace-Fsubtr)/Fdiv)*100;
deltaF2 = ((rawTrace2-Fsubtr)/Fdiv)*100;

R_Fdiv    = mean(mean(R_rawTrace(:,1:1000)));
R_Fsubtr  = mean(R_rawTrace(:,1001:2000),2);    % using this time window in plotwin for baseline
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


close all
fig1 = figure(1);
suptitle(whichMouse)
set(gcf, 'Position', [217 511 1089 518])

fig_1 = subplot(1,2,1);
errorbar_patch(plotWin,m_plot,s_plot,[0 1 0]);
errorbar_patch(plotWin,R_m_plot,R_s_plot,[1 0 0]);
title('GCaMP signals during arm reaching in') %approach') %to reward')
legend('GCaMP','tdTom')
xlabel('time - point of arm reaching in') %approach (ms)') %approach %retreat %consumption
ylabel('dF/F (%)')

fig_2 = subplot(1,2,2);
errorbar_patch(plotWin,m_plot2,s_plot2,[0 1 0]);
errorbar_patch(plotWin,R_m_plot2,R_s_plot2,[1 0 0]);
title('GCaMP signals during reward delivery') %retreat')
legend('GCaMP','tdTom')
xlabel('time - point of reward delivery') %retreat (ms)') %approach %retreat %consumption
ylabel('dF/F (%)')
fig_2.YLim = fig_1.YLim;

% raster plot (10ms window average)
rasterWin = 10;
a = reshape(deltaF2(:,1:(length(plotWin)-1)), length(acttime_index), rasterWin, []);
%a = reshape(R_deltaF(:,1:(length(plotWin)-1)), length(acttime_index), rasterWin, []); %tdTom
b = squeeze(mean(a,2));

fig2 = figure(2);
imagesc(b)
title(['Raster (' num2str(rasterWin) 'ms bin): ' whichMouse])
c = plotWin(1:1000:end)/1000;
fig2.CurrentAxes.XTickLabel = c(2:end);
fig2.CurrentAxes.XLabel.String = 'time - point of reward delivery'; %retreat (s)';
%fig2.CurrentAxes.XLabel.String = 'time - point of arm reaching in';
fig2.CurrentAxes.YLabel.String = 'trial';

%cutoff = 25000;
%fig3 = figure(3);
%scatter(vel(1:cutoff), GCaMP(veltime_ind(1:25000)), 'k.')
%xlabel('Velocity')
%ylabel('GCaMP')

%% save (7)
cd ./proc/Analyzed_Data
%saveas(fig1, ['MoSeqFP_' whichMouse '_' whichSession '_MoSeq.tif'])
saveas(fig1, ['MoSeqFP_' whichMouse '_ArmReachRewDelivery_MoSeq.tif'])
%saveas(fig2, ['MoSeqFP_raster_' whichMouse '_' whichSession '_MoSeq.tif'])
saveas(fig2, ['MoSeqFP_raster_' whichMouse '_ArmReach_MoSeq.tif'])%RewardDelivery_MoSeq.tif'])
saveas(fig3, ['MoSeqFP_velGCaMP_' whichMouse '_' whichSession '_MoSeq.tif'])

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

if(0)
LEDtest(1).eventname='on1';
LEDtest(2).eventname='on2';
LEDtest(3).eventname='on3';
LEDtest(4).eventname='on4';
LEDtest(1).rgbonframe = 66;
LEDtest(2).rgbonframe = 170;
LEDtest(3).rgbonframe = 22660;
LEDtest(4).rgbonframe = 22764;
LEDtest(1).depthonframe = 134;
LEDtest(2).depthonframe = 342;
LEDtest(3).depthonframe = 45158;
LEDtest(4).depthonframe = 45369;
LEDtest(1).rgbonTstepIdx = 6399;
LEDtest(2).rgbonTstepIdx = 13599;
LEDtest(3).rgbonTstepIdx = 1557188;
LEDtest(4).rgbonTstepIdx = 1564115;
LEDtest(1).depthonTstepIdx = 6399;
LEDtest(2).depthonTstepIdx = 13531;
LEDtest(3).depthonTstepIdx = 1556721;
LEDtest(4).depthonTstepIdx = 1563889;

% mouse in = 220; 

% LED on
cross_where = crossing(zscore(ch02));
LED_on      = cross_where(1:2:end); % FP index
LED_off     = cross_where(2:2:end); % FP index
end

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