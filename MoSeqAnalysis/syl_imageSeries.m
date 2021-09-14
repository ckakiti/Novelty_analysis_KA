%% create image series for particular syllables
clear
close all
clc

% load MiceIndex with MoSeq labels synced to DLC video
cd('/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/Behavior/others/MoSeq_combine3L/');
load('MiceIndex_wLabels_update.mat') % created in add_MoSeqLabels.m
curr_MoSeqLabels = Mice(1).ExpDay.moseq_align;

% load DLC avi video
cd('/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/Behavior/others/Standard setup/Capoeira_newLegos/Videos/Au/20190413');
vn='Au_190413_01_rgb_Converted.avi';
video=VideoReader(vn);
vidFrames = video.Duration*video.FrameRate;

% load DLC labels
cd('/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/Behavior/stimulus/Au/novel1')
fn=dir([vn(1:end-4) '*.csv']);
fn=fn.name;
Labels = csvread(fn,3,0);

% check to make sure length of MoSeq labels, 
% DLC video, and DLC labels are the same
sameLength = (length(curr_MoSeqLabels)==vidFrames) ...
    & (vidFrames==length(Labels));
if(~sameLength)
    disp('!! labels and video frames not the same length !!')
end

%% find all frames where syllable of interest is expressed
currSyl = 14; %%%%%%%%% syllable of interest
currSyl_loc = find(curr_MoSeqLabels==currSyl);
disp(['total frames for currSyl: ' num2str(length(currSyl_loc))])

% find breaks between syllable bouts
sylBout_break = find(diff(currSyl_loc)>1);
curr_break = 4; % which syllable bout to focus on for image series 
                %syl14/br2 %syl79/br5 - %syl14/br8 %syl79/br3 - %syl79/br5 %syl14/br4
curr_sylBout = currSyl_loc(sylBout_break(curr_break)+1:sylBout_break(curr_break+1));

% get center of mass (nose/tail avg) for all syllable frames
com_sylFrames = [(Labels(currSyl_loc,2)+Labels(currSyl_loc,11))/2 ...
                 (Labels(currSyl_loc,3)+Labels(currSyl_loc,12))/2];

frame_cropSeries = [];
%for frameiter = 1:length(currSyl_loc)
for frameiter = sylBout_break(curr_break)+1:sylBout_break(curr_break+1)
    %disp([num2str(frameiter) '/' num2str(length(currSyl_loc))])
    video.CurrentTime = (currSyl_loc(frameiter)-1)/video.FrameRate;
    frame=readFrame(video);
    
    cropSize = 40; % frames to display around mouse c.o.m.
    curr_com = com_sylFrames(frameiter,:);
    frame_crop = frame((curr_com(2)-cropSize):(curr_com(2)+cropSize),...
                       (curr_com(1)-cropSize):(curr_com(1)+cropSize),...
                       :);
    if(any(currSyl_loc(frameiter)==curr_sylBout))
        disp(frameiter)
        frame_crop_ex = frame_crop;
        frame_cropSeries = [frame_cropSeries frame_crop_ex];
    end
    
    %imshow(frame)
    %imshow(frame_crop)
end
disp('end')

close all
imshow(frame_cropSeries)