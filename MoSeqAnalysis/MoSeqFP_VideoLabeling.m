clear
clc
close all

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry_DLC/Miami/temp
load('Miami_190125_nidaq.mat')

vn='Miami_190125_depth.mp4';
raw_video=VideoReader(vn);

final_video = VideoWriter([vn(1:end-4) '_FP_Labeled.avi']);
final_video.FrameRate = raw_video.FrameRate;
open(final_video);
videolength=round(raw_video.Duration.*raw_video.FrameRate,0);

imageslen=3987;%55474;
imagesfile = '/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry_DLC/Miami/temp/jpgs_lessProcessed/';
%imagesfile='/Users/yuxie/Downloads/Temp/';

framenum = 900;  % >1 if video was cropped for MoSeq analysis
mywaitbar = waitbar(0,[num2str(round(100*framenum/videolength)) '%' '    |    ' num2str(framenum) '/' num2str(videolength)]);

for frameiter=1:imageslen-framenum

    if ~hasFrame(raw_video)
        disp(['no frame: ', num2str(frameiter)])
        break
    end

    rawframe=readFrame(raw_video);

    %sidebar1=uint8(255.*ones(500,65,3));
    %finalframe=cat(2,sidebar1,rawframe);

    %sidebar2=uint8(255.*ones(500,66,3));
    %finalframe=cat(2,finalframe,sidebar2); 23

    finalframe=rawframe;
    sidebar1=uint8(255.*ones(156,23,3));
    FPbar=imread([imagesfile 'fig_' num2str(framenum) '.jpg']);
    FPbar=cat(2,FPbar,sidebar1);
    
    finalframe=cat(1,finalframe,FPbar);

    writeVideo(final_video,finalframe);

    framenum = framenum + 1;
    waitbar(framenum/videolength,mywaitbar,[num2str(round(100*framenum/videolength)) '%' '    |    ' num2str(framenum) '/' num2str(videolength)]);
end
close(mywaitbar);
close(final_video);