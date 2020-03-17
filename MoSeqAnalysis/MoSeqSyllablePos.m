clear
clc
close all


dataset_name = 'Dataset_20190723';
analysis_len = 10; % analyze first 10 min or all 30 min of each session
disp(dataset_name)
disp([num2str(analysis_len) 'min'])


cd(['/media/alex/DataDrive1/MoSeqData/' dataset_name '/MoSeq'])
% cd(['/media/alex/DataDrive1/MoSeqData/' dataset_name '/Data'])

load('MoSeqDataFrame.mat')
Syllablebinedge=[-6,-0.5:1:99.5];

Mice_Index_path=['/media/alex/DataDrive1/MoSeqData/' dataset_name '/MoSeq/MiceIndex.mat'];
% Mice_Index_path=['/media/alex/DataDrive1/MoSeqData/' dataset_name '/Data/MiceIndex.mat'];
load(Mice_Index_path)
% run(Mice_Index_path);


load(['GeneralAnalysis_' dataset_name '_' num2str(analysis_len) 'min.mat'])

MarkerSize=5;
Fsize=20;

disp('done')

%%
% IntSyllable=[45 20 9 33 14 17 64 62 60 95 85 94 36 50 55 11 0 10 16 2];%[42 94 39 38 70 13];
% IntSyllable=[36 37 74 18 24];
% IntSyllable = [15 41 72 1 4 23 63 75 31 52 44 33 49 47 8];
IntSyllable = 9;% [57 52];

close all
SyllableDis=figure;
title('Syllable Position Distribution','FontSize',Fsize)
xlabel('x position (mm)','FontSize',Fsize)
ylabel('y position (mm)','FontSize',Fsize)
set(gca,'ydir','reverse')
set(SyllableDis, 'position', [0 0 1000 850]);

for syliter=1:length(IntSyllable)
    XPos=MoSeqDataFrame.centroid_x_mm(MoSeqDataFrame.model_label==IntSyllable(syliter));
    YPos=MoSeqDataFrame.centroid_y_mm(MoSeqDataFrame.model_label==IntSyllable(syliter));
    scatter(XPos,YPos,MarkerSize,'filled')
    hold on
    
    disp(IntSyllable(syliter))
    pause
end
%
legend(strcat('Syllable  ',num2str(IntSyllable')))
title('Spatial Distribution of Syllable Expression','FontSize',Fsize)
xlabel('x position (mm)','FontSize',Fsize)
ylabel('y position (mm)','FontSize',Fsize)
xlim([-400 300])
ylim([-400 300])
set(gca,'ydir','reverse')
set(SyllableDis, 'position', [0 0 1000 850]);

%%
IntSyllable = [15 41 72 1 4 23 63 75 31 52 44 33 49 47 8];

close all
SyllableDis=figure;
% title('Syllable Position Distribution','FontSize',Fsize)
% xlabel('x position (mm)','FontSize',Fsize)
% ylabel('y position (mm)','FontSize',Fsize)
% set(gca,'ydir','reverse')
% set(SyllableDis, 'position', [0 0 1000 850]);

for syliter=1:length(IntSyllable)
    XPos=MoSeqDataFrame.centroid_x_mm(MoSeqDataFrame.model_label==IntSyllable(syliter));
    YPos=MoSeqDataFrame.centroid_y_mm(MoSeqDataFrame.model_label==IntSyllable(syliter));
    
    subplot(3,5,syliter)
    scatter(XPos,YPos,2,'k','filled')
    title(num2str(IntSyllable(syliter)))
    xlim([-400 300])
    ylim([-400 300])
    axis square
    set(gca,'ydir','reverse')
end

set(SyllableDis, 'position', [0 0 1850 850]);

%%
saveas(SyllableDis, 'Dataset_20191007_SyllablePosDis_syl57_syl52.tif')

%% divide plots into expression during novelty vs habituation

%  *********WARNING: TAKES FOREVER TO RUN**************

IntSyllable = 52; %[57 52]; 9;
syliter = IntSyllable(1);

xPos = [];
yPos = [];

allSyl = find(MoSeqDataFrame.model_label==syliter);

tic
for miceiter=1:length(Mice)
    disp(miceiter)
    
    for dayiter = 2 %G1_Days
        curr_uuid = Mice(miceiter).ExpDay(dayiter).MSid;
        find_pos  = ismember(MoSeqDataFrame.uuid,curr_uuid,'rows');
        find_pos  = find(find_pos);
        crop_pos  = find_pos(1:(analysis_len*fps*60));
        
        curr_model_label = MoSeqDataFrame.model_label(crop_pos);
        where_syliter    = find(curr_model_label == syliter);
        intersect_allSyl = intersect(allSyl, crop_pos(where_syliter));
        
        add_xPos = MoSeqDataFrame.centroid_x_mm(intersect_allSyl);
        add_yPos = MoSeqDataFrame.centroid_y_mm(intersect_allSyl);
        
        xPos = [xPos add_xPos];
        yPos = [yPos add_yPos];
    end
end
toc

if(0)
    csvwrite([dataset_name '_syl' num2str(IntSyllable) '_N1_xPos_' ...
        num2str(analysis_len) 'min.csv'],xPos)
    csvwrite([dataset_name '_syl' num2str(IntSyllable) '_N1_yPos_' ...
        num2str(analysis_len) 'min.csv'],yPos)
end

%% plotting interesting syllable on specific day
IntSyllable = 9; %[52 57]; %9
whichDay = 'N1';
xPos_1 = csvread([dataset_name '_syl' num2str(IntSyllable(1)) '_' whichDay ...
    '_xPos_' num2str(analysis_len) 'min.csv']);
yPos_1 = csvread([dataset_name '_syl' num2str(IntSyllable(1)) '_' whichDay ...
    '_yPos_' num2str(analysis_len) 'min.csv']);
% xPos_2 = csvread([dataset_name '_syl' num2str(IntSyllable(2)) '_' whichDay ...
%    '_xPos_' num2str(analysis_len) 'min.csv']);
% yPos_2 = csvread([dataset_name '_syl' num2str(IntSyllable(2)) '_' whichDay ...
%    '_yPos_' num2str(analysis_len) 'min.csv']);

close all
SyllableDis=figure;

hold on
s1 = scatter(xPos_1,yPos_1,MarkerSize,'filled');
% s2 = scatter(xPos_2,yPos_2,MarkerSize,'filled');

title(['Spatial Distribution of Syllable Expression (' num2str(analysis_len) 'min)'], ...
    'FontSize',Fsize)
xlabel('x position (mm)','FontSize',Fsize)
ylabel('y position (mm)','FontSize',Fsize)
legend(s1(1), ['Syllable  ', num2str(IntSyllable(1))])
% legend([s1(1) s2(1)], {['Syllable  ', num2str(IntSyllable(1))] ...
%                        ['Syllable ', num2str(IntSyllable(2))]})
xlim([-400 300])
ylim([-400 300])
set(gca,'ydir','reverse')
set(SyllableDis, 'position', [0 0 1000 850]);

disp('plotted')

if(0)
    saveas(SyllableDis, [dataset_name '_SyllablePosDis_syl' ...
        num2str(IntSyllable(1)) '_' whichDay '_' num2str(analysis_len) 'min.tif'])

    saveas(SyllableDis, [dataset_name '_SyllablePosDis_syl' ...
        num2str(IntSyllable(1)) '_syl' num2str(IntSyllable(2)) '_' ...
        whichDay '_' num2str(analysis_len) 'min.tif'])
end



