% Now using a loop to generate image, by implementing matrix operation using index information in actalignedusage could improve speed
% when syllable index is -5 it is a 'none' type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

%cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/CvsS_180831_DLC/
%cd /media/alex/DataDrive1/MoSeqData/Iku_photometry2/Iku_photometry2_MoSeq/Nashville/
cd /media/alex/DataDrive1/MoSeqData/Dataset_20191007/Data
% cd('/media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq')

fps=30;
PlotWidth=200;%500;
BarHeight=5;
cmap=jet(100);
fsize=24;

load('MoSeqDataFrame.mat');
% Mice_Index_path='./MiceIndex.m';%'../MiceIndex/MiceIndex.m';
% run(Mice_Index_path);
load('MiceIndex.mat') % to get this, need to run extract_uuid.m
%cd ../MoSeq/DRILLS_MoSeq
%load('test_NearObj_ts.mat')

AnalysisDay=3;      % 3 = first novelty day
%G1_Mice=1;
detectCond = cat(1, Mice.novelty);
G1_Mice=find(detectCond=='C')';%[1 2 3 4];
G2_Mice=find(detectCond=='S')';%[5 6 7 8];
whichMouse = 3;           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['whichMouse: ' num2str(whichMouse)])

frameCutoff = 18000;

trim_frame_start=899;
% AllActLabels=csvread('CvsS_poke_labels_N1_byHand.csv',1,3);%2);
AllAct_file = dir('*poke_labels_N1.csv');
AllAct=csvread(AllAct_file.name,0,0);
AllActLabels=AllAct(:,3);
AllActLabels=AllActLabels-trim_frame_start;


for miceiter=1:length(Mice)
    Mice(miceiter).datanum = length(find(AllAct(:,1)==miceiter));
end

Mice(1).index=1;
for miceiter=2:length(Mice)
    Mice(miceiter).index=Mice(miceiter-1).index+Mice(miceiter-1).datanum;
    endindex=Mice(miceiter).index+Mice(miceiter).datanum;
end

if endindex ~= length(AllActLabels)+1
    disp('index calculation error');
%     error('index calculation error');
end

for miceiter=1:length(Mice)
    Mice(miceiter).ExpDay(AnalysisDay).act=AllActLabels(Mice(miceiter).index:Mice(miceiter).index+Mice(miceiter).datanum-1,:);
end

disp('section 1')

%% get labels for all days

for miceiter=1:length(Mice)
    disp(Mice(miceiter).name)
    
    for day_iter = 1:length(Mice(miceiter).ExpDay)
        disp(day_iter)
        
        % find MSid index
        MSidindex=1;
        for indexiter=1:size(MoSeqDataFrame.session_uuid,1)
            if strcmp(MoSeqDataFrame.session_uuid(indexiter,:),Mice(miceiter).ExpDay(day_iter).MSid)
                break
            end
            MSidindex=MSidindex+1;
            if MSidindex==size(MoSeqDataFrame.session_uuid,1)+1
                error('MSid not found');
            end
        end
        
        Labels=double(MoSeqDataFrame.labels{MSidindex});
        Mice(miceiter).ExpDay(day_iter).labels = Labels;
    end

end
if(0)
    save('Mice_wLabels', 'Mice')
end


allLabels = []; % cropped by 10 min
allSyls   = [];

for miceiter = 1:length(Mice)
    for dayiter = 3 %3:length(Mice(miceiter).ExpDay)   %%%%%%%%%%%%%%%%%%%%%%%
        currLabels = cat(1, Mice(miceiter).ExpDay(dayiter).labels);
        allSyls = [allSyls currLabels];
        
%         frameCutoff = length(currLabels);              %%%%%%%%%%%%%%%%%%%%%%%
        [counts, centers] = hist(currLabels(1:frameCutoff),[-1 1:99]);
        freq = counts./sum(counts);
        allLabels = [allLabels; freq];
    end
end
disp(['Frame cutoff: ' num2str(frameCutoff)])
allLabels_S = allLabels(G2_Mice,:);
allLabels_C = allLabels(G1_Mice,:);

if(0)
    csvwrite('Dataset_20191007_sylExpr_allN_stim_30min.csv', allLabels_S)
    
    csvwrite('Dataset_20191007_sylExpr_N1.csv', allLabels)
    csvwrite('Dataset_20191007_sylExpr_N1_stim_30min.csv', allLabels_S)
    csvwrite('Dataset_20191007_sylExpr_N1_stim_10min.csv', allLabels_S)
end

disp('end')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic

for miceiter=1:length(Mice)


    % find MSid index
    MSidindex=1;
    for indexiter=1:size(MoSeqDataFrame.session_uuid,1)
        if strcmp(MoSeqDataFrame.session_uuid(indexiter,:),Mice(miceiter).ExpDay(AnalysisDay).MSid)
            break
        end
        MSidindex=MSidindex+1;
        if MSidindex==size(MoSeqDataFrame.session_uuid,1)+1
            error('MSid not found');
        end
    end

    Labels=double(MoSeqDataFrame.labels{MSidindex});
    Mice(miceiter).ExpDay(AnalysisDay).labels = Labels;
    labellen=length(Labels);

    Mice(miceiter).ExpDay(AnalysisDay).rasterimage=uint8(255.*ones(20,PlotWidth,3));
    Mice(miceiter).ExpDay(AnalysisDay).rasterimage=insertText(Mice(miceiter).ExpDay(AnalysisDay).rasterimage,[round(PlotWidth./2-35,0),0],[Mice(miceiter).name '  Day: ' num2str(AnalysisDay)],'BoxOpacity',0,'TextColor','black');

    for actiter=1:Mice(miceiter).datanum

        framenum=Mice(miceiter).ExpDay(AnalysisDay).act(actiter);%,2);

        %Adding Syllable Bar
        middle_x=round(PlotWidth./2,0);
        middle_y=1;

        Syllableline=255.*ones(3,PlotWidth);

        for lineiter=1:PlotWidth
            rpos=lineiter-middle_x;


            plotframenum=framenum+rpos;

            if plotframenum<=0 || plotframenum>labellen
                syllableindex=-5;
            else
                syllableindex=Labels(plotframenum);
            end


            % Adding to a matrix actalignedusage that store all syllable usage around the labeled point
            % actalignedpoint is the labeled point
            Mice(miceiter).ExpDay(AnalysisDay).actalignedusage(actiter,lineiter)=syllableindex;
            Mice(miceiter).ExpDay(AnalysisDay).actalignedpoint(actiter)=middle_x;

            % Calculatiing syllabel color corrisponding to color map
            if syllableindex==-5
                syllablecolor=[0 0 0];
            elseif syllableindex<100 && syllableindex>=0
                syllablecolor=cmap(syllableindex+1,:);
            else
                error(['syllableindex out of bund, index=' num2str(syllableindex)]);
            end

            % Adding current position pointer in the raster plot
            if rpos>-1 && rpos<1
                syllablecolor=[1 0 0];
            end

            Syllableline(:,lineiter)=Syllableline(:,lineiter).*(syllablecolor');
        end

        barline=uint8(zeros(1,PlotWidth,3));
        barline(1,:,1)=Syllableline(1,:);
        barline(1,:,2)=Syllableline(2,:);
        barline(1,:,3)=Syllableline(3,:);

        Syllablebar=barline;
        for appenditer=1:BarHeight-1
            Syllablebar=cat(1,Syllablebar,barline);
        end
        Mice(miceiter).ExpDay(AnalysisDay).rasterimage=cat(1,Mice(miceiter).ExpDay(AnalysisDay).rasterimage,Syllablebar);
    end
end


% Calculating general usage only used when plot width is cosistent across mice
Generalactalignedusage=[];
for miceiter=1:length(Mice)
    Generalactalignedusage=cat(1,Generalactalignedusage,Mice(miceiter).ExpDay(AnalysisDay).actalignedusage);
end

% General act aligned syllable usage GAASU, 
% meaning GAASU(r , c) = count of syllable r-1 in all data at frame c
Syllablebinedge=[-6,-0.5:1:99.5];
GAASU=zeros(101,PlotWidth);
for lineiter =1:PlotWidth
    GAASU(:,lineiter)=(histcounts(Generalactalignedusage(:,lineiter),Syllablebinedge))';
end
GAASU(size(GAASU,1)+1,:)=GAASU(1,:);
GAASU(1,:)=[];

% Percentage usage
PGAASU=GAASU./size(Generalactalignedusage,1);

% % Accumulated act aligned syllable usage AAASU
% AAASU=GAASU;
% for rowiter=size(AAASU,1):-1:2
%     AAASU(rowiter,:)=sum(AAASU(1:rowiter,:));
% end

% Calculating Group1 Group2 usage only used when plot width is cosistent across mice
G1actalignedusage=[];
for miceiter=G1_Mice%(whichMouse)        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G1actalignedusage=cat(1,G1actalignedusage,Mice(miceiter).ExpDay(AnalysisDay).actalignedusage);
end

G2actalignedusage=[];
for miceiter=G2_Mice%(whichMouse)        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G2actalignedusage=cat(1,G2actalignedusage,Mice(miceiter).ExpDay(AnalysisDay).actalignedusage);
end

% Group1 act aligned syllable usage G1AASU,  
% meaning G1AASU(r , c) = count of syllable r-1 in all data at frame c
Syllablebinedge=[-6,-0.5:1:99.5];
G1AASU=zeros(101,PlotWidth);
for lineiter =1:PlotWidth
    G1AASU(:,lineiter)=(histcounts(G1actalignedusage(:,lineiter),Syllablebinedge))';
end
G1AASU(size(G1AASU,1)+1,:)=G1AASU(1,:);
G1AASU(1,:)=[];

% Percentage usage
PG1AASU=G1AASU./size(G1actalignedusage,1);


% Group2 act aligned syllable usage G2AASU,
% meaning G2AASU(r , c) = count of syllable r-1 in all data at frame c
G2AASU=zeros(101,PlotWidth);
for lineiter =1:PlotWidth
    G2AASU(:,lineiter)=(histcounts(G2actalignedusage(:,lineiter),Syllablebinedge))';
end
G2AASU(size(G2AASU,1)+1,:)=G2AASU(1,:);
G2AASU(1,:)=[];

% Percentage usage
PG2AASU=G2AASU./size(G2actalignedusage,1);

disp('section 2')

% save variables for future analysis
%save(['ActAlignedPercentage_Day' num2str(AnalysisDay)], 'middle_x', 'PG1AASU', 'PG2AASU')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
timeline=((1:PlotWidth)-round(PlotWidth./2,0))./fps;

% Plot_Gactalignedusage=figure;
% Gareahandle=area(timeline,PGAASU','LineWidth',0.05);
% % areahandle(82).FaceColor=cmap(82,:);
% title('Syllable Usage Aligned by Human Labeled interaction','FontSize',fsize)
% xlabel('Time (s)','FontSize',fsize)
% ylabel('Usage Percentage','FontSize',fsize)
% axis([timeline(1),timeline(end),0,1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_G1actalignedusage=figure;
set(Plot_G1actalignedusage,'Position',[46 307 1874 390])
G1areahandle=area(timeline,PG1AASU','LineWidth',0.05);
set(gca,'FontSize',16)
% title('Syllable Usage of Contextual Novelty Mice','FontSize',fsize)
title('Contextual Novelty','FontSize',fsize)
% title(['Syllable Usage of ' Mice(G1_Mice(whichMouse)).name],'FontSize',fsize)
xlabel('Time (s)','FontSize',fsize)
ylabel('Syllable Usage (%)','FontSize',fsize)
axis([timeline(1),timeline(end),0,1])
set(gca,'YTick',[0 0.5 1])

% saveas(Plot_G1actalignedusage,'Planets_syllableUsage_cont.tif')
% saveas(Plot_G1actalignedusage,['Capoeira_' Mice(G1_Mice(whichMouse)).name '_syllableUsage_cont.tif'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_G2actalignedusage=figure;
set(Plot_G2actalignedusage,'Position',[46 107 1874 390])
G2areahandle=area(timeline,PG2AASU','LineWidth',0.05);
set(gca,'FontSize',16)
% title('Syllable Usage of Stimulus Novelty Mice','FontSize',fsize)
title('Stimulus Novelty','FontSize',fsize)
% title(['Syllable Usage of ' Mice(G2_Mice(whichMouse)).name],'FontSize',fsize)
xlabel('Time (s)','FontSize',fsize)
ylabel('Syllable Usage (%)','FontSize',fsize)
axis([timeline(1),timeline(end),0,1])
set(gca,'YTick',[0 0.5 1])

% saveas(Plot_G2actalignedusage,'Planets_syllableUsage_stim.tif')
% saveas(Plot_G2actalignedusage,['Capoeira_' Mice(G2_Mice(whichMouse)).name '_syllableUsage_stim.tif'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is for matching the syllable color with the color in the raster plot and labeled videos
% for coloriter=1:100
%     areahandle(coloriter).FaceColor=cmap(coloriter,:);
% end
% areahandle(101).FaceColor=[0 0 0];

disp('plotted')

%% plots of individual syllables for statistical analysis
close all
currSyl = 60; %58 %71 %39 94 15
currWinLen = 80;
currWin = middle_x-currWinLen:middle_x+currWinLen;
PG1AASU_currSyl = PG1AASU(currSyl,currWin);
PG1AASU_currSyl = smooth(PG1AASU_currSyl,10);
PG2AASU_currSyl = PG2AASU(currSyl,currWin);
PG2AASU_currSyl = smooth(PG2AASU_currSyl,10);

fig1 = figure(1);
set(fig1, 'Position', [46 345 1875 390])
hold on
group1 = plot(currWin-middle_x, PG1AASU_currSyl,'r');
group2 = plot(currWin-middle_x, PG2AASU_currSyl,'b');
ylim([0 1])
title(['Syllable ', num2str(currSyl)])
legend([group1, group2], {'Contextual' 'Stimulus'})

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir('EventAlignedAnalysis')
cd('EventAlignedAnalysis')

% Saving raster plot
for miceiter=1:length(Mice)
    imwrite(Mice(miceiter).ExpDay(AnalysisDay).rasterimage,[Mice(miceiter).name '_interaction_raster_plot.png']);
end
cd ..
toc

