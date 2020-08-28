%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all
clc

cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/Analysis-tools/B-SOID/Dataset_200501')

fps=15;
PlotWidth=200;
BarHeight=5;
syllablesTotal = 12; %%%%%%%%%%%%%%% depends on bsoid output (# of clusters) %%% 
cmap=jet(syllablesTotal);
fsize=24;
AnalysisDay = 3; % Day 3 = N1

% load results from bsoid, output of bsoid_master_vlp2.m
% (after model training, then prediction on 15fps DLC files)
load('test_mdl_200507_allData_15fpsLabel.mat',...
    'csv_files','f_15fps_test','labels_test15fps');
bsoid_labels_reshape = reshape(labels_test15fps, 2, [])';

% load poke labels for H2 and N1 for Dataset_20190723
Mice_Capoeira = load('PokesApproaches_Capoeira.mat','Mice');
Mice_Chess = load('PokesApproaches_Chess.mat','Mice');
Mice_Hiking = load('PokesApproaches_Hiking.mat','Mice');
AllMice = [Mice_Capoeira.Mice Mice_Chess.Mice Mice_Hiking.Mice];

% create structures for H2 and N1 (will merge later)
AllMice_H2 = struct('name', {AllMice.name}, ...
    'novelty', {AllMice.novelty}, ...
    'day', {2}, ...
    'act', {AllMice.Pokes_Day2});
AllMice_N1 = struct('name', {AllMice.name}, ...
    'novelty', {AllMice.novelty}, ...
    'day', {3}, ...
    'act', {AllMice.Pokes_Day3});

% merge H2 and N1 array
AllMice_merge = [AllMice_H2 AllMice_N1];

% sort alphabetically so they match order from bsoid analysis
[sort_names, sort_idx] = sort({AllMice_merge.name});
AllMice_sorted = AllMice_merge(sort_idx);

% add bsoid labels to structure
[AllMice_sorted.labels] = labels_test15fps{:};
[AllMice_sorted.features] = f_15fps_test{:};

datanums = cellfun(@length, {AllMice_sorted.act}, 'UniformOutput', false);
[AllMice_sorted.datanum] = datanums{:};
Mice_OneDay = AllMice_sorted(cat(1,AllMice_sorted.day)==AnalysisDay);

% Mice = AllMice_sorted; %H2 and N1     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mice=Mice_OneDay; % only H2 or only N1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detectCond = cat(1,Mice.novelty);
G1_Mice=find(detectCond=='C')'; % relative to Mice_bsoid/bsoid_files
G2_Mice=find(detectCond=='S')';
whichMouse = 9;           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['whichMouse: ' num2str(whichMouse) ...
    ' (' Mice(G1_Mice(whichMouse)).name '/' ...
         Mice(G2_Mice(whichMouse)).name ')'])

disp('section 1')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic

for miceiter=1:length(Mice)
    % extract bsoid labels
    Labels = Mice(miceiter).labels;
    labellen=length(Labels);

    Mice(miceiter).rasterimage=uint8(255.*ones(20,PlotWidth,3));
    Mice(miceiter).rasterimage=insertText(Mice(miceiter).rasterimage,...
        [round(PlotWidth./2-35,0),0],...
        [Mice(miceiter).name '  Day: ' num2str(AnalysisDay)],...
        'BoxOpacity',0,'TextColor','black');

    for actiter=1:Mice(miceiter).datanum

        framenum=Mice(miceiter).act(actiter);

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
            Mice(miceiter).actalignedusage(actiter,lineiter)=syllableindex;
            Mice(miceiter).actalignedpoint(actiter)=middle_x;

            % Calculatiing syllable color corrisponding to color map
            if syllableindex==-5
                syllablecolor=[0 0 0];
            elseif syllableindex<syllablesTotal && syllableindex>=0
                syllablecolor=cmap(syllableindex+1,:);
            else
                error(['syllableindex out of bound, index=' num2str(syllableindex)]);
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
        Mice(miceiter).rasterimage=cat(1,Mice(miceiter).rasterimage,Syllablebar);
    end
end


% Calculating general usage only used when plot width is cosistent across mice
Generalactalignedusage=[];
for miceiter=1:length(Mice)
    Generalactalignedusage=cat(1,Generalactalignedusage,Mice(miceiter).actalignedusage);
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
for miceiter=G1_Mice(whichMouse)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G1actalignedusage=cat(1,G1actalignedusage,Mice(miceiter).actalignedusage);
end

G2actalignedusage=[];
for miceiter=G2_Mice(whichMouse)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G2actalignedusage=cat(1,G2actalignedusage,Mice(miceiter).actalignedusage);
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
if(0)
    save(['ActAlignedPercentage_poke_PW' num2str(PlotWidth)], ...
        'middle_x', 'PG1AASU', 'PG2AASU')
end

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
% title('Contextual Novelty','FontSize',fsize)
title(['Syllable Usage of ' Mice(G1_Mice(whichMouse)).name],'FontSize',fsize)
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
% title('Stimulus Novelty','FontSize',fsize)
title(['Syllable Usage of ' Mice(G2_Mice(whichMouse)).name],'FontSize',fsize)
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
% if(0)
    saveas(Plot_G1actalignedusage,['BSOID_Dataset_20200501_C' ...
        num2str(whichMouse) '_' Mice(G1_Mice(whichMouse)).name '_syllableUsage.tif'])
%         'cont_syllableUsage.tif'])
    saveas(Plot_G2actalignedusage,['BSOID_Dataset_20200501_S' ...
        num2str(whichMouse) '_' Mice(G2_Mice(whichMouse)).name '_syllableUsage.tif'])
%         'stim_syllableUsage.tif'])
% end