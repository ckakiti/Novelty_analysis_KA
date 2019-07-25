% Note: Have not tested yet
% Update: tested and modified by CKA 190501
clear
close all
clc

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Mitsuko_photometry_190617/

mouseList = dir;
for mouseiter=3:length(mouseList)
    cd(mouseList(mouseiter).name)
    disp(['Mouse: ' mouseList(mouseiter).name])
    
    dateList = dir;
    for dateiter = 4:5 %3:length(dateList)
        cd(dateList(dateiter).name)
        
        filelist = dir('*.mp4');
        fileCurr = filelist.name;
        
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
        
        cd ..
    end
    
    cd ..
end