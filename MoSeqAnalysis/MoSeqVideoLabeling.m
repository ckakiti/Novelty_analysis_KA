% Notice: some of MoSeq videos and label length do not match

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

Labeling_mice=[2]
Labeling_days=[2]
BarHeight=30;
Mice_Index_path='/media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq/MiceIndex_isIn_head_depth.mat';
% Mice_Index_path='/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Configural_set/MiceIndex.mat';

%run(Mice_Index_path); %for .m file
load(Mice_Index_path); %for .mat file


tic
cd '/media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq/'
% cd '/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Configural_set/'
load('MoSeqDataFrame.mat');
% load('ActAlignedPercentage_poke_PW200.mat')
% cd('Videos');
% cd extra

%  for miceiter=1:length(Mice)
for miceiter=Labeling_mice
    cd Capoeira_MoSeq
    cd(Mice(miceiter).name);

    % for dayiter=1:length(Mice(miceiter).ExpDay)
    for dayiter=Labeling_days %%%%%%%%%%%%%%%%%%%%%%%
        cd(Mice(miceiter).ExpDay(dayiter).date(:))
        sessionName = dir('session*');
        %cd(sessionName(Labeling_days).name) %%%%%%%%%%%%%%%%%%%%%%%
        cd(sessionName(1).name)
        cd proc
        
        % find MSid index
        MSidindex=1;
        for indexiter=1:size(MoSeqDataFrame.session_uuid,1)
            if strcmp(MoSeqDataFrame.session_uuid(indexiter,:),Mice(miceiter).ExpDay(dayiter).MSid)
                break
            end
            MSidindex=MSidindex+1;
            if MSidindex==size(MoSeqDataFrame.session_uuid,1)+1
                error('MSid not found');
            end
        end

        Labels=double(MoSeqDataFrame.labels{MSidindex});
%         Labels = Mice(Labeling_mice).labels';
        
        vn='results_00.mp4';
%         vn='Au_190413_01_rgb_Converted.avi';
        Changepoint=(abs(Labels(2:end)-Labels(1:end-1))>0);
        Changepoint=[Changepoint zeros(1,6)];
        cmap=jet(100);
%         cmap=jet(12);

        raw_video=VideoReader(vn);
        final_video = VideoWriter([vn(1:end-4) '_Labeled.avi']);
        final_video.FrameRate = raw_video.FrameRate;
        open(final_video);
        videolength=round(raw_video.Duration.*raw_video.FrameRate,0);

        labellen=length(Labels);

        % if videolength<labellen-5 || videolength>labellen+5
        %     error('video and label mot match');
        % end

        framenum = 1;
        mywaitbar = waitbar(0,[num2str(round(100*framenum/videolength)) '%' '    |    ' num2str(framenum) '/' num2str(videolength)]);
        % while hasFrame(raw_video) 
        for frameiter=1:labellen

            if ~hasFrame(raw_video)
                break
            end

            rawframe=readFrame(raw_video);

            finalframe=rawframe;
            %Adding text
            textpos = [0,100;    % Frame Number
                    0,120;
                    0,140;    % Syllable Number
                    0,160;];

            insertedtext = {'Frame Number:';
                            num2str(framenum);
                            'Syllable';
                            num2str(Labels(framenum))};

            finalframe = insertText(finalframe,textpos,insertedtext,'BoxOpacity',0,'TextColor','white');

            %Adding Changepoint Notation

            if Changepoint(framenum)~=0 ||  Changepoint(framenum+1)~=0 || Changepoint(framenum+2)~=0
                finalframe = insertShape(finalframe,'Filledcircle',[40 90 10],'Color','Red');
                finalframe = insertShape(finalframe,'Filledcircle',[40 raw_video.Height-10 10],'Color','Red');
            end

            %Adding Syllable Bar
            
            middle_x=round(raw_video.Width./2,0);
            middle_y=raw_video.Height+1;

            Syllableline=255.*ones(3,raw_video.Width);
            for lineiter=1:raw_video.Width
                rpos=lineiter-middle_x;
                if rpos>-1 && rpos<1
                    syllablecolor=[1 0 0];
                else
                    plotframenum=framenum+rpos;

                    if plotframenum<=0 || plotframenum>labellen
                        syllableindex=-5;
                    else
                        syllableindex=Labels(plotframenum);
                    end

                    if syllableindex==-5
                        syllablecolor=[0 0 0];
                    elseif syllableindex<100 && syllableindex>=0
                        syllablecolor=cmap(syllableindex+1,:);
                    else
                        error(['syllableindex out of bund, index=' num2str(syllableindex)]);
                    end
                end

                Syllableline(:,lineiter)=Syllableline(:,lineiter).*(syllablecolor');
            end

            videoline=zeros(1,raw_video.Width,3);
            videoline(1,:,1)=Syllableline(1,:);
            videoline(1,:,2)=Syllableline(2,:);
            videoline(1,:,3)=Syllableline(3,:);

            Syllablebar=videoline;
            for appenditer=1:BarHeight-1
                Syllablebar=cat(1,Syllablebar,videoline);
            end

            finalframe=cat(1,finalframe,Syllablebar);
            writeVideo(final_video,finalframe);

            framenum = framenum + 1;
            waitbar(framenum/videolength,mywaitbar,[num2str(round(100*framenum/videolength)) '%' '    |    ' num2str(framenum) '/' num2str(videolength)]);
        end

        close(mywaitbar);
        close(final_video);
        toc

        cd ..
    end

    cd ../..
    cd ..
end
 cd ..
 
%%
pokes = csvread('Hiking_poke_labels_N1.csv',0,0);
pokes_Curr = pokes(pokes(:,1)==1,:);

pokes_curr_depth = pokes_Curr(:,3);
pokes_too_early = find((pokes_curr_depth-899)<0,1,'last');
pokes_curr_syl_MoSeq = Labels(pokes_curr_depth(pokes_too_early+1:end)-899)';
pokes_curr_syl_manual = pokes_Curr(pokes_too_early+1:end,4);

mismatch = find(pokes_curr_syl_MoSeq~=pokes_curr_syl_manual);