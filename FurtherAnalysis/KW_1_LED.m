% Read the video file (.mp4) %%%%%%%%%%%%%%%%%%%
clear
close all
clc
%cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry2_DLC/Nashville/190423
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/FP_Test/FP_Test_190422/

% import video information
%tic
for i = 1:1
    
    folder = pwd;
    rgbfile = dir('*rgb.mp4');
    videoObject = VideoReader(rgbfile(1).name);
    numberOfFrames = round(videoObject.Duration*videoObject.FrameRate);
    vidHeight = videoObject.Height;
    vidWidth = videoObject.Width;
    
end

% import LED position (from MarkObjPos.m)
cd Analyzed_Data
load('Arena_Obj_Pos.mat')
cd ..

%% Find when the light was on %%%%%%%%%%%%%%%%%%%%%%%%%%
save_time = [5000:5000:numberOfFrames-1 numberOfFrames];
sync_light = zeros(1,numberOfFrames);
for frame = 1:numberOfFrames %3000
    disp(frame)
    
    % Extract one frame and make it grayscale
    
    %thisFrame = read(videoObject, frame);
    videoObject.currentTime = (frame-1)/videoObject.FrameRate;
    thisFrame = readFrame(videoObject);
    %thisFrame1 = rgb2gray(thisFrame);
    
    
    % Cut out the area outside of the LED
    
    % thisFrame1(1:420,:)= [];
    % thisFrame1(:,70:end)= [];
    % thisFrame1(:,1:45)= [];
    % thisFrame1(20:end,:)= [];
    thisFrame1 = thisFrame(LEDpos(2):LEDpos(4), LEDpos(1):LEDpos(3),:);
    thisFrame1 = rgb2gray(thisFrame1);
    
    
    % Find how bright the LED is on this frame
    
    sync_light(1,frame) = round(mean(mean(thisFrame1)));
    
    if (any(frame==save_time))
        disp(['save: ' num2str(frame)])
        csvwrite(['save_' num2str(frame)], sync_light)
    end
end

%clear frame thisFrame thisFrame1

%% Plot when the light was on %%%%%%%%%%%%%%%%%%%%%%%%%%
crosses = crossing(zscore(sync_light(1:5001)));

close all
figure(1)
hold on
plot(zscore(sync_light(1:5001)));
line([crosses; crosses], repmat([-1.5 2.5], length(crosses), 1)')

LED_on = crosses(1:2:end);
diff_LED_on = diff(LED_on);
missing_on = find(diff_LED_on>1.5*mode(diff_LED_on));
if(~isempty(missing_on))
    disp(['missing LED: ' num2str(length(missing_on))])
    line([LED_on(missing_on); LED_on(missing_on)], ...
        repmat([-2.5 3.5], length(missing_on), 1)', 'color', 'r')
end

for addIn = 1:length(missing_on)
    missingCurr = ( LED_on(missing_on(addIn)) + LED_on(missing_on(addIn)+1) ) / 2;
    LED_on      = [LED_on(1:missing_on(addIn)) missingCurr LED_on(missing_on(addIn)+1:end)];
    missing_on  = missing_on+1;
end
figure(2)
plot(LED_on, 'k.-')

%% save crossing (LED on/off) variable
csvwrite('LED_on', LED_on)

