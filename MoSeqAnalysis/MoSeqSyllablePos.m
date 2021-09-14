clear
clc
close all


dataset_name = 'Dataset_20190723';
analysis_len = 30; % analyze first 10 min or all 30 min of each session
disp(dataset_name)
disp([num2str(analysis_len) 'min'])


MSDF_path = ['/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/Behavior/others/Standard setup/' ...
    'CombineAnalysis/Just-in-case files/' dataset_name '/MoSeqDataFrame.mat'];
% '/media/alex/DataDrive1/MoSeqData/' dataset_name '/MoSeq';
load(MSDF_path)
Syllablebinedge=[-6,-0.5:1:99.5];

% Mice_Index_path=['/media/alex/DataDrive1/MoSeqData/' dataset_name '/MoSeq/MiceIndex.mat'];
Mice_Index_path = ['/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/Behavior/others/Standard setup/' ...
   'CombineAnalysis/Just-in-case files/Dataset_20190723/MiceIndex.mat'];
load(Mice_Index_path)


%load(['GeneralAnalysis_' dataset_name '_' num2str(analysis_len) 'min.mat'])

MarkerSize=5;
Fsize=20;

disp('done')

%%
group1 = {'Au', 'Ginga', 'Negativa', ...
    'Bishop', 'Knight', 'Rook', ...
    'Appalachian', 'Continental', 'Long'};
group2 = {'Esquiva', 'MeiaLua', 'Queixada', ...
    'King', 'Pawn', 'Queen', ...
    'Arizona', 'JohnMuir', 'Pacific'};

mouse_names_all = cellstr(MoSeqDataFrame.SubjectName);
mouse_names     = unique(mouse_names_all);

% assign group names
group1_loc  = all(ismember(mouse_names_all, group1),2);
group2_loc  = all(ismember(mouse_names_all, group2),2);
MoSeqDataFrame.group(group1_loc,:) = repmat('group1 ',length(find(group1_loc)),1);
MoSeqDataFrame.group(group2_loc,:) = repmat('group2 ',length(find(group2_loc)),1);

% get x/y position just for N1
centroid_x = [];
centroid_y = [];
for mouseiter = 1:length(mouse_names)
    % extract session uuids for single mouse
    curr_mouse_loc  = all(ismember(mouse_names_all, mouse_names(mouseiter)),2);
    uuid_all        = MoSeqDataFrame.uuid(curr_mouse_loc,:);
    uuid_unique     = unique(cellstr(uuid_all));
    
    % get x/y position for current mouse
    justN1 = all(ismember(uuid_all, uuid_unique(3)),2);
    centroid_x_all = MoSeqDataFrame.centroid_x_mm(curr_mouse_loc);
    centroid_y_all = MoSeqDataFrame.centroid_y_mm(curr_mouse_loc);
    centroid_x_N1 = centroid_x_all(justN1);
    centroid_y_N1 = centroid_y_all(justN1);
    
    centroid_x = [centroid_x centroid_x_N1];
    centroid_y = [centroid_y centroid_y_N1];
end

%%
% IntSyllable=[45 20 9 33 14 17 64 62 60 95 85 94 36 50 55 11 0 10 16 2];%[42 94 39 38 70 13];
% IntSyllable=[36 37 74 18 24];
% IntSyllable = [15 41 72 1 4 23 63 75 31 52 44 33 49 47 8];
IntSyllable = 9;% 9; 1; %[57 52];

close all
SyllableDis=figure;
title('Syllable Position Distribution','FontSize',Fsize)
xlabel('x position (mm)','FontSize',Fsize)
ylabel('y position (mm)','FontSize',Fsize)
set(gca,'ydir','reverse')
set(SyllableDis, 'position', [0 0 1000 850]);

for syliter=1:length(IntSyllable)
    % *** plot position just for N1 ***
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

%% subplot of multiple interesting syllables
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

%% before plotting position during novelty vs habituation:
%    create cell arrays for centroid_x_mm/centroid_y_mm
%    and indices for each uuid within the larger MSDF
clear uuid_index curr_uuid
indexiter = 1;
uuid_index(indexiter).uuid = MoSeqDataFrame.uuid(1,:);
uuid_index(indexiter).idx_start = 1;
uuid_index(indexiter).idx_end = [];

for MSDFiter = 1:length(MoSeqDataFrame.centroid_x_mm)
    curr_uuid = MoSeqDataFrame.uuid(MSDFiter,:);
    
    if(~strcmp(curr_uuid, uuid_index(indexiter).uuid))
        uuid_index(indexiter).idx_end = MSDFiter-1;
        
        indexiter = indexiter + 1;
        uuid_index(indexiter).uuid = curr_uuid;
        uuid_index(indexiter).idx_start = MSDFiter;
    end
end
uuid_index(end).idx_end = MSDFiter;

framenum = num2cell(cat(1,uuid_index.idx_end)+1-cat(1,uuid_index.idx_start));
[uuid_index.framenum] = deal(framenum{:});

[uuid_sort_string, uuid_sort_index] = sort({uuid_index.uuid});
uuid_index_sorted = uuid_index(uuid_sort_index);

for centroid_iter = 1:length(uuid_index_sorted)
    curr_start = uuid_index_sorted(centroid_iter).idx_start;
    curr_end   = uuid_index_sorted(centroid_iter).idx_end;
    
    uuid_index_sorted(centroid_iter).centroid_x_mm = ...
        MoSeqDataFrame.centroid_x_mm(curr_start:curr_end);
    uuid_index_sorted(centroid_iter).centroid_y_mm = ...
        MoSeqDataFrame.centroid_y_mm(curr_start:curr_end);
    
    uuid_index_sorted(centroid_iter).model_label = ...
        MoSeqDataFrame.model_label(curr_start:curr_end);
end
  
%% divide plots into expression during novelty vs habituation
IntSyllable = 3; %[57 52]; 9;
syliter = IntSyllable(1);

xPos = [];
yPos = [];

for miceiter=1:length(Mice)
    disp(miceiter)
    daynum = length(Mice(miceiter).ExpDay);
    
    for dayiter = 2 %1:daynum %G1_Days
        curr_uuid  = Mice(miceiter).ExpDay(dayiter).MSid;
        where_uuid = find(ismember(cat(1,uuid_index_sorted.uuid),...
            curr_uuid,'rows'));
        
        curr_labels = uuid_index_sorted(where_uuid).model_label;
        crop_labels = curr_labels(1:(analysis_len*fps*60)); %%%%%%%%%%%%%%%%%%%%%%
        
        where_syl = find(crop_labels==syliter);
        curr_xPos = uuid_index_sorted(where_uuid).centroid_x_mm(where_syl);
        curr_yPos = uuid_index_sorted(where_uuid).centroid_y_mm(where_syl);

        xPos = [xPos curr_xPos];
        yPos = [yPos curr_yPos];
    end
end

%% plotting interesting syllable on specific day
% IntSyllable = 9; %[52 57]; %9
whichDay = 'H2'; %'N1'; 'H2'; 'allDays';
xPos_1 = xPos;
yPos_1 = yPos;
% xPos_1 = csvread([dataset_name '_syl' num2str(IntSyllable(1)) '_' whichDay ...
%     '_xPos_' num2str(analysis_len) 'min.csv']);
% yPos_1 = csvread([dataset_name '_syl' num2str(IntSyllable(1)) '_' whichDay ...
%     '_yPos_' num2str(analysis_len) 'min.csv']);
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

%% previous code to divide plots into novelty vs habituation:
%  *********WARNING: TAKES FOREVER TO RUN**************
IntSyllable = 9; %[57 52]; 9;
syliter = IntSyllable(1);

xPos = [];
yPos = [];

allSyl = find(MoSeqDataFrame.model_label==syliter);

tic
for miceiter=1:length(Mice)
    disp(miceiter)
    
    for dayiter = 3 %G1_Days
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
                              