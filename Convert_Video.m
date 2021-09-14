% Note: Have not tested yet (yuxie)
% Update: tested and modified by CKA 190501
clear
close all
clc

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Aviary_DLC/
%cd('/home/alex/Dropbox (Uchida Lab)/Korleki Akiti/ForMelissa/')

files    = dir;
whichDir = [files.isdir];
nameDir  = files(whichDir);
nameDir  = {nameDir.name};
nameDir(ismember(nameDir,{'.','..','temp'})) = [];
mouseList = nameDir;

for mouseiter=1:length(mouseList)
    cd(mouseList{mouseiter})
    disp(['Mouse: ' mouseList{mouseiter}])
    
    files    = dir;
    whichDir = [files.isdir];
    nameDir  = files(whichDir);
    nameDir  = {nameDir.name};
    nameDir(ismember(nameDir,{'.','..','25fps','temp'})) = [];
    dateList = nameDir;
    
    for dateiter = 8%1:length(dateList) %%%%%%%%%%%%%%%%%%%%%%%
        cd(dateList{dateiter})
        
        filelist = dir('*.mp4');
        for fileiter = 1:length(filelist)
            fileCurr = filelist(fileiter).name;
            
            vn = fileCurr;
            disp(['Analyzing: ' vn]);
            
            raw_video=VideoReader(vn);
            final_video = VideoWriter([vn(1:end-4) '_Converted.avi']);
            final_video.FrameRate = raw_video.FrameRate;
            open(final_video);
            videolength=round(raw_video.Duration.*raw_video.FrameRate);
            
            tic
            framenum = 1;
            h = waitbar(0,[num2str(round(100*framenum/videolength)) '%' '    |    ' num2str(framenum) '/' num2str(videolength)]);
            while hasFrame(raw_video)
                rawframe=readFrame(raw_video);
                writeVideo(final_video,rawframe);%finalframe);
                framenum = framenum + 1;
                waitbar(framenum/videolength,h,[num2str(round(100*framenum/videolength)) '%' '    |    ' num2str(framenum) '/' num2str(videolength)]);
            end
            close(h);
            toc
            close(final_video);
        end
        
        cd ..
    end
    
    cd ..
end

%% create video clips containing synchronized, concatenated rgb+depth frames
clear
close all
clc
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Capoeira_DLC/temp/Au_190413


% load timestamp files for synchronization
depthts_file = dir('*depth_ts.txt');
rgbts_file = dir('*rgb_ts.txt');
filename1 = rgbts_file.name;
filename2 = depthts_file.name;
delimiter = ' ';
formatSpec = '%f%[^\n\r]'; %'%*q%f%[^\n\r]';
fileID1 = fopen(filename1,'r');
fileID2 = fopen(filename2,'r');
dataArray1 = textscan(fileID1, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, ...
    'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
dataArray2 = textscan(fileID2, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, ...
    'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID1);
fclose(fileID2);
rgbts = [dataArray1{1}];
depthts = [dataArray2{1}];


% load video files
filelistRGB = dir('*Converted*.mp4');%rgb_25fps.mp4');  % original rgb output file (raw)
filelistDEP = dir('*results_00_Labeled_cropped.m4v'); % depth file from .dat file, extracted by MoSeq
fileCurrRGB = filelistRGB.name;
fileCurrDEP = filelistDEP.name;
vnRGB = fileCurrRGB;
vnDEP = fileCurrDEP;
raw_video_rgb = VideoReader(vnRGB);
raw_video_dep = VideoReader(vnDEP);


% cropping for DLC video
cropLeft_rgb   = 99;
cropRight_rgb  = 455;
cropTop_rgb    = 34-25;
cropBottom_rgb = 385+25;

% cropping for depth video
cropLeft_dep   = 105;
cropRight_dep  = 461;
cropTop_dep    = 69;
cropBottom_dep = 528;%470;

% cropping for depth aligned
cropLeft_align   = 0;
cropRight_align  = 80;
cropTop_align    = 0;
cropBottom_align = 80;

% iterate through rgb frames, then match depth to rgb
frames = (1:300) + 1150; % desired interval within rgb vid
min_frame = find(rgbts>depthts(901),1); % determined by number of depth frames cropped during extraction
frame_offset = [repmat(875,135,1); repmat(868,165,1)];

final_video = VideoWriter([vnRGB(1:end-4) '_temp.avi']); % _merged % _cropped %%%%%%%%%%%%%%%%
final_video.FrameRate = raw_video_rgb.FrameRate;
open(final_video);

for frameiter = 1:length(frames) 
    disp(frameiter)
    frame_rgb = frames(frameiter) + min_frame;
    
    curr_rgb_ts = rgbts(frame_rgb);
    frame_dep   = frameiter;
%     frame_dep   = find(depthts>curr_rgb_ts,1)-frame_offset(frameiter); %900;
%     [temp, frame_dep] = min(abs(depthts-curr_rgb_ts));
    
    raw_video_rgb.CurrentTime = (frame_rgb-1)/raw_video_rgb.FrameRate;
    raw_video_dep.CurrentTime = (frame_dep-1)/raw_video_dep.FrameRate; % -875
    frame_rgb = readFrame(raw_video_rgb);
    frame_dep = readFrame(raw_video_dep);
    frameCrop_rgb = frame_rgb;
    frameCrop_dep = frame_dep;
    
    frameCrop_rgb(:,1:cropLeft_rgb,:) = [];
    frameCrop_rgb(1:cropTop_rgb,:,:)  = [];
    frameCrop_rgb(:,cropRight_rgb-cropLeft_rgb:end,:) = [];
    frameCrop_rgb(cropBottom_rgb-cropTop_rgb:end,:,:) = [];
    
    frameCrop_dep(:,1:cropLeft_dep,:) = [];
    frameCrop_dep(1:cropTop_dep,:,:)  = [];
    frameCrop_dep(:,cropRight_dep-cropLeft_dep:end,:) = [];
    frameCrop_dep(cropBottom_dep-cropTop_dep:end,:,:) = [];
%     frameCrop_dep(:,1:cropLeft_align,:) = [];
%     frameCrop_dep(1:cropTop_align,:,:)  = [];
%     frameCrop_dep(:,cropRight_align-cropLeft_align:end,:) = [];
%     frameCrop_dep(cropBottom_align-cropTop_align:end,:,:) = [];

%     newFrame = [frameCrop_rgb frameCrop_dep];
    newFrame=frameCrop_dep;
    writeVideo(final_video,newFrame);
%     imshow(newFrame)
end
close(final_video);
disp('end')

%% iterate through depth frames, match rgb to depth (lesser quality in rgb part of video)
frames = 1:1000 + 0; % desired interval within depth vid
min_frame = 900;     % determined by number of frames cropped during extraction

for frameiter = 1:length(frames) 
    frame_dep = frames(frameiter);
    
    curr_dep_frame = frames(frameiter)+min_frame;
    curr_dep_ts    = depthts(curr_dep_frame);
    
    frame_rgb = find(rgbts>curr_dep_ts,1);
    
    raw_video_rgb.CurrentTime = (frame_rgb-1)/raw_video_rgb.FrameRate;
    raw_video_dep.CurrentTime = (frame_dep-1)/raw_video_dep.FrameRate;
    frame_rgb = readFrame(raw_video_rgb);
    frame_dep = readFrame(raw_video_dep);
    frameCrop_rgb = frame_rgb;
    frameCrop_dep = frame_dep;
    
    frameCrop_rgb(:,1:cropLeft_rgb,:) = [];
    frameCrop_rgb(1:cropTop_rgb,:,:)  = [];
    frameCrop_rgb(:,cropRight_rgb-cropLeft_rgb:end,:) = [];
    frameCrop_rgb(cropBottom_rgb-cropTop_rgb:end,:,:) = [];
    
    frameCrop_dep(:,1:cropLeft_dep,:) = [];
    frameCrop_dep(1:cropTop_dep,:,:)  = [];
    frameCrop_dep(:,cropRight_dep-cropLeft_dep:end,:) = [];
    frameCrop_dep(cropBottom_dep-cropTop_dep:end,:,:) = [];
    
    newFrame = [frameCrop_rgb frameCrop_dep];
    imshow(newFrame)
end
