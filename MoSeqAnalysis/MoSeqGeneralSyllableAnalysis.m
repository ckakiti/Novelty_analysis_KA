%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

% setName = 'Standard_CapoeiraHikingChess';
setName = 'Dataset_20190723';%'7day_preexposure';

% Mice_Index_path='/Users/yuxie/Dropbox/YuXie/CvsS_180831/CvsS_180831_MoSeq/Mice_Index.m';
% Mice_Index_path='/media/alex/DataDrive1/MoSeqData/7day_preexposure_MoSeq/Mice_Index.m';
% Mice_Index_path='/media/alex/DataDrive1/MoSeqData/Capoeira/Capoeira_MoSeq/Mice_Index.m';
Mice_Index_path = '/media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq/MiceIndex.mat';
%Mice_Index_path = ['/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/Behavior/Standard setup/' ...
%   'CombineAnalysis/Just-in-case files/Dataset_20190723/MiceIndex.mat'];
% run(Mice_Index_path);
load(Mice_Index_path)

detectCond = cat(1, Mice.novelty);
cond2 = find(detectCond=='C');
cond1 = find(detectCond=='S');

G1_Mice = cond2';%1:3;%[2 4 6];%[1 2 3 4]; %7:12;
G2_Mice = cond1';%4:6;%[1 3 5];%[5 6 7 8]; %1:6;
G1_Days = [3 4 5 6]; %3:7;
G2_Days = [3 4 5 6]; %3:7;

% G3 Base line (habituation)
G3_Mice = 1:length(Mice);
G3_Days = [1 2];

%cd /media/alex/DataDrive1/MoSeqData/CvsS_20180831_MoSeq %7day_preexposure_MoSeq
% cd /media/alex/DataDrive1/MoSeqData/Capoeira/Capoeira_MoSeq
cd /media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq
%cd(['/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/Behavior/Standard setup/' ...
%    'CombineAnalysis/Just-in-case files/Dataset_20190723/'])
load('MoSeqDataFrame_75.mat')
%load('GeneralUsage.mat') %%%%%%%%%%%%%%
cmap=jet(100);
fps=30;
Syllablebinedge=[-6,-0.5:1:99.5];
frameCutoff = 18000;

% load('NearObj_ts.mat')

disp('section 1')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic

% Generate Transition Matrix

for miceiter=1:length(Mice)

    for dayiter=1:length(Mice(miceiter).ExpDay)

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
        frameCutoff = length(MoSeqDataFrame.labels{MSidindex}); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(frameCutoff)
        Labels=double(MoSeqDataFrame.labels{MSidindex}(1:frameCutoff));
        labellen=length(Labels);

        % Calculate Empirical Bigram Transition Matrix (self transition included)
        % BiMatrixSum is the sum of Bi-transition, exclude window with 'none' type
        Mice(miceiter).ExpDay(dayiter).BiMatrix=zeros(100,100);
        Mice(miceiter).ExpDay(dayiter).BiMatrixSum=0;
        
        for frameiter=1:labellen-1

            readwindow=Labels(frameiter:frameiter+1);
            if sum(readwindow<0)>0
                continue
            end

            Mice(miceiter).ExpDay(dayiter).BiMatrix(readwindow(1)+1,readwindow(2)+1) = ...
                1 + Mice(miceiter).ExpDay(dayiter).BiMatrix(readwindow(1)+1,readwindow(2)+1);
            Mice(miceiter).ExpDay(dayiter).BiMatrixSum=Mice(miceiter).ExpDay(dayiter).BiMatrixSum+1;

        end

        % Calculate syllable usage count, usage(1) 'none' syllable, usage(2) number of syllable 0 is used...
        % usagesum is the sum of all frames except 'none' type
        Mice(miceiter).ExpDay(dayiter).usage=histcounts(Labels,Syllablebinedge);
        Mice(miceiter).ExpDay(dayiter).usagesum=sum(Mice(miceiter).ExpDay(dayiter).usage(2:end));
        
    end

end

% Calculate General Usage of all experiments
GUsage=zeros(1,101);
for miceiter=1:length(Mice)
    for dayiter=1:length(Mice(miceiter).ExpDay)
        GUsage = GUsage + Mice(miceiter).ExpDay(dayiter).usage;   
    end
end
PGUsage=GUsage./sum(GUsage);

[GSortedusage,GSortedusageindex]=sort(PGUsage,'descend');

AccGSortedusage=zeros(1,101);
for usageiter=1:length(GSortedusage)
    AccGSortedusage(usageiter)=sum(GSortedusage(1:usageiter));
end

% Calculate Usage of Group1, Group2, and Group3
G1Usage=zeros(1,101);
for miceiter=G1_Mice
    for dayiter=G1_Days
        G1Usage = G1Usage + Mice(miceiter).ExpDay(dayiter).usage; 
    end
end
PG1Usage=G1Usage./sum(G1Usage);

G2Usage=zeros(1,101);
for miceiter=G2_Mice
    for dayiter=G2_Days
        if(dayiter>length(Mice(miceiter).ExpDay))
            continue
        else
            G2Usage = G2Usage + Mice(miceiter).ExpDay(dayiter).usage;
        end
    end
end
PG2Usage=G2Usage./sum(G2Usage);

G3Usage=zeros(1,101);
for miceiter=G3_Mice
    for dayiter=G3_Days
        G3Usage = G3Usage + Mice(miceiter).ExpDay(dayiter).usage;
    end
end
PG3Usage=G3Usage./sum(G3Usage);

G2vsG1usage=(PG2Usage-PG1Usage)./(PG2Usage+PG1Usage);
[G2vsG1Sortedusage,G2vsG1Sortedusageindex]=sort(G2vsG1usage,'descend');

% Calculate General Bigram Transition Matrix of all experiments
GBM=zeros(100,100);
for miceiter=1:length(Mice)
    for dayiter=1:length(Mice(miceiter).ExpDay)
        GBM = GBM + Mice(miceiter).ExpDay(dayiter).BiMatrix;
    end
end
InterGBM=GBM-GBM.*diag(ones(1,100));     % subtract self transition
PGBM=InterGBM./sum(sum(InterGBM));       % Normalization


% Calculate Bigram Transition Matrix of Group1
G1BM=zeros(100,100);
for miceiter=G1_Mice
    for dayiter=G1_Days
        G1BM = G1BM + Mice(miceiter).ExpDay(dayiter).BiMatrix;
    end
end
InterG1BM=G1BM-G1BM.*diag(ones(1,100));     % subtract self transition
PG1BM=InterG1BM./sum(sum(InterG1BM));       % Normalization

% Calculate Bigram Transition Matrix of Group2
G2BM=zeros(100,100);
for miceiter=G2_Mice
    for dayiter=G2_Days
        if(dayiter>length(Mice(miceiter).ExpDay))
            continue
        else
            G2BM = G2BM + Mice(miceiter).ExpDay(dayiter).BiMatrix;
        end
    end
end
InterG2BM=G2BM-G2BM.*diag(ones(1,100));     % subtract self transition
PG2BM=InterG2BM./sum(sum(InterG2BM));       % Normalization

disp('section 2')

if(0)
    clear MoSeqDataFrame
    save('GeneralAnalysis_Dataset_20190723_30min')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Making plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
X=0:100;
SyllablesX=-1:99;
fsize=16;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_GeneralUsage=figure(1);
plot(X,GSortedusage,'b.-', 'LineWidth',1.5)
hold on
line([0 100], [0.01 0.01], 'color', 'k', 'linestyle', '--')
title(['General Syllable Usage (sorted by usage, ', setName, ')'],'FontSize',fsize, ...
    'Interpreter', 'none')
ylabel('Fractional usage','FontSize',fsize)
xlabel('Syllable Number','FontSize',fsize)
xticks(X);
xticklabels(SyllablesX(GSortedusageindex));
xlim([0 60])
ylim([0 0.055])
set(Plot_GeneralUsage, 'Position', [46 353 1870 450])
%saveas(Plot_GeneralUsage, [setName '_generalUsage_showCutoff.tif'])

%% just plot top 20 most used syllables
Plot_GeneralUsage_Top20=figure(2);
hold on
plot(X(1:20),GSortedusage(1:20),'LineWidth',1.5, 'Marker', 'o', 'LineStyle', '-')
line([X(1) X(20)], [1/length(SyllablesX) 1/length(SyllablesX)], 'Color', 'b', 'LineStyle', '--')
title(['Syllable Usage Across All Mice/Days (', setName, ')'],'FontSize',fsize, ...
    'Interpreter', 'none')
ylim([0 max(PGUsage)+0.001])
ylabel('Fractional Usage','FontSize',fsize)
xlabel('Syllable Number','FontSize',fsize)
xticks(X(1:20));
xticklabels(SyllablesX(GSortedusageindex(1:20)));
set(Plot_GeneralUsage_Top20, 'Position', [340 386 716 450])
hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clc

amt_explained = 0.9;
syl_amt_explained = find(AccGSortedusage>=amt_explained,1)-1;

Plot_AccGeneralUsage=figure(3);
hold on
plot(AccGSortedusage,'k','LineWidth',1.5)
line([0 length(X)], [amt_explained amt_explained], ...
    'Color', 'k', 'LineStyle', '--')
line([syl_amt_explained syl_amt_explained], [0 1], ...
    'Color', 'k', 'LineStyle', '--')
text(syl_amt_explained+2, amt_explained-0.05, ...
    ['n=' num2str(syl_amt_explained)], 'FontSize', 20)

title(['Accumulated Syllable Usage (', setName, ')'],'FontSize',fsize,...
    'Interpreter', 'none')
ylabel('Fractional Usage','FontSize',fsize)
xlabel('Syllable Rank','FontSize',fsize)
xlim([0 length(X)])
ylim([0 1.01])
set(gca, 'FontSize', 16)
set(gcf, 'Position', [320 110 700 500])
% xticks(X);

if(0)
    saveas(Plot_AccGeneralUsage, [setName '_AccSyllableUsage_sorted_allFrames.tif'])
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_UsageCompare=figure(4);
plot(X,PG1Usage,'LineWidth',1.5)
hold on
plot(X,PG2Usage,'LineWidth',1.5)

legend('Contextual Novelty','Stimulus Novelty')
title(['Syllable Usage Comparison of Contextual/Stimulus Novely Mice (', setName, ')'],...
    'FontSize',fsize,'Interpreter', 'none')
ylabel('Percentage','FontSize',fsize)
xlabel('Syllables','FontSize',fsize)
xticks(X);
xticklabels(SyllablesX);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot_UsageCompare_sorted=figure(5);
plot(X,PG1Usage(G2vsG1Sortedusageindex),'LineWidth',1.5)
hold on
plot(X,PG2Usage(G2vsG1Sortedusageindex),'LineWidth',1.5)
hold on
plot(X,PG3Usage(G2vsG1Sortedusageindex),'LineWidth',1,'Color','Black')

legend({'Contextual Novelty','Stimulus Novelty','Habituattion'},'FontSize',fsize)
title(['Syllable Usage Comparison of Contextual/Stimulus Novely Mice (', setName, ...
    ') (Sorted by stimulus novelty enrichment)'],'FontSize',fsize,'Interpreter', 'none')
ylabel('Fractional usage','FontSize',fsize)
xlabel('Syllables','FontSize',fsize)
xticks(X);
xticklabels(SyllablesX(G2vsG1Sortedusageindex));
%saveas(Plot_UsageCompare_sorted, [setName '_UsageCompare_sorted.tif'])

%% MoSeq statistical analysis

% step 1: identify syllables with >1% usage across all mice+days,
%         pool with syllables that are used >1% within radius of object
GsortedCutoff = GSortedusageindex(GSortedusage>=0.01);
% SyllablesX(GSortedusageindex(GSortedusage>=0.01));

cropFrames = 900;
disp(['Cropped frames: ', num2str(cropFrames)])
addLabels = [];
for miceiter=1:length(Mice)
    
    for dayiter=1:length(Mice(miceiter).ExpDay)
        disp([Mice(miceiter).name, ': ', Mice(miceiter).ExpDay(dayiter).date])
        
        % find MSid index
        MSidindex=1;
        for indexiter=1:size(MoSeqDataFrame.session_uuid,1)
            if strcmp(MoSeqDataFrame.session_uuid(indexiter,:),...
                      Mice(miceiter).ExpDay(dayiter).MSid)
                break
            end
            MSidindex=MSidindex+1;
            if MSidindex==size(MoSeqDataFrame.session_uuid,1)+1
                error('MSid not found');
            end
        end

        frameCutoff = length(MoSeqDataFrame.labels{MSidindex}); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Labels=double(MoSeqDataFrame.labels{MSidindex}(1:frameCutoff));
        
        %cd(['./vids/', Mice(miceiter).name])
        cd('./Capoeira_MoSeq')
        cd(Mice(miceiter).name)
        cd(Mice(miceiter).ExpDay(dayiter).date)
        
        sessionName = dir('session*');
        cd(sessionName.name)
        
        allDepth_ts = dir('*depth_ts*');
        
        % load depthts
        filename = allDepth_ts(1).name; %%%%%%%%%%
        delimiter = ' ';
        formatSpec = '%f%[^\n\r]';%'%*q%f%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, ...
            'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
        fclose(fileID);
        depthts = [dataArray{1}];
        
        nearObj_ts = Mice(miceiter).ExpDay(dayiter).NearObj_ts;
        nearObj_index = zeros(length(nearObj_ts), 1);
        for nearObjiter=1:length(nearObj_ts)
            nearObj_index(nearObjiter,1)=find(depthts>nearObj_ts(nearObjiter),1);
        end

        nearObj_index = nearObj_index-cropFrames; % cropped frames during moseq-extract
        nearObj_index(nearObj_index<=0)=[];
        nearObj_index(nearObj_index>frameCutoff)=[];
        if(mod(length(nearObj_index),2))
            disp('odd nearObj_index length')
            nearObj_index = nearObj_index(1:end-1);
        end
        
        for depthiter = 1:2:length(nearObj_index)
            toAdd = Labels(nearObj_index(depthiter):nearObj_index(depthiter+1));
            addLabels = [addLabels toAdd];
        end
        
        cd ../../../..
    end
end



addLabelsUsage=histcounts(addLabels,Syllablebinedge);
PaddLabelsUsage=addLabelsUsage./sum(addLabelsUsage);

[PaddLabelsUsageSort,PaddLabelsUsageSortIndex]=sort(PaddLabelsUsage,'descend');
addLabelsCutoff = PaddLabelsUsageSortIndex(PaddLabelsUsageSort>=0.01);
% SyllablesX(PaddLabelsUsageSortIndex(PaddLabelsUsageSort>=0.01));

% include syllables expressed >1% of time within specified radius of obj
% GsortedCutoffCombine = unique([GsortedCutoff addLabelsCutoff]);
GsortedCutoffCombine = GsortedCutoff; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% crop sorted vector to include only syllables >1%
isIn = ones(length(G2vsG1Sortedusageindex),1);
for sortediter = 1:length(G2vsG1Sortedusageindex)
    if(isempty(find(G2vsG1Sortedusageindex(sortediter)==GsortedCutoffCombine,1)))
        isIn(sortediter,1) = 0;
    end
end
G2vsG1Sortedusageindex_crop = G2vsG1Sortedusageindex(isIn==1);

%% plot usage just for syllables with >1% expression
PGUsage_crop = PGUsage(GsortedCutoffCombine);
[GSortedusage_crop,GSortedusageindex_crop]=sort(PGUsage_crop,'descend');
X_crop = 0:length(GSortedusage_crop)-1;

close all
Plot_GeneralUsage=figure(1);
plot(X_crop,GSortedusage_crop,'b.-', 'LineWidth',1.5)
hold on
line([0 length(GSortedusage_crop)], [0.01 0.01], 'color', 'k', 'linestyle', '--')
title(['General Syllable Usage (syllables >1% in general or within radius, ', setName, ')'],...
    'FontSize',fsize, 'Interpreter', 'none')
ylabel('Fractional usage','FontSize',fsize)
xlabel('Syllable Number','FontSize',fsize)
xticks(X_crop);
xticklabels(SyllablesX(GsortedCutoffCombine(GSortedusageindex_crop)));
%xlim([0 60])
ylim([0 0.055])
set(Plot_GeneralUsage, 'Position', [46 353 1870 450])
%saveas(Plot_GeneralUsage, [setName '_generalUsage_crop.tif'])

%% plot syllable enrichment without errorbar
close all
X_cutoff = 0:length(G2vsG1Sortedusageindex_crop)-1;

Plot_UsageCompare_cutoff=figure(6);
plot(X_cutoff,PG1Usage(G2vsG1Sortedusageindex_crop),'LineWidth',1.5)
hold on
plot(X_cutoff,PG2Usage(G2vsG1Sortedusageindex_crop),'LineWidth',1.5)
plot(X_cutoff,PG3Usage(G2vsG1Sortedusageindex_crop),'LineWidth',1,'Color','Black')

legend({'Contextual Novelty','Stimulus Novelty','Habituation'},'FontSize',fsize)
title(['Syllable Usage Comparison of Contextual/Stimulus Novely Mice (', setName, ...
    ') (Sorted by stimulus novelty enrichment)'],'FontSize',fsize,'Interpreter', 'none')
ylabel('Fractional usage','FontSize',fsize)
xlabel('Syllables','FontSize',fsize)
xticks(X_cutoff);
xticklabels(SyllablesX(G2vsG1Sortedusageindex_crop));
set(Plot_UsageCompare_cutoff, 'Position', [46 353 1870 450])

%saveas(Plot_UsageCompare_cutoff, [setName '_UsageCompare_Crop.tif'])

%% divide syllable usage per mouse for statistical tests (ttest2, bonferroni correction)
G12Usage_split=zeros(length(G3_Mice),101);
for miceiter=G3_Mice
    for dayiter=G1_Days
        if(dayiter>length(Mice(miceiter).ExpDay))
            continue
        else
            G12Usage_split(miceiter,:) = G12Usage_split(miceiter,:) + ...
                Mice(miceiter).ExpDay(dayiter).usage; 
        end
    end
end
PG12Usage_split=G12Usage_split./sum(G12Usage_split,2);

G3Usage_split=zeros(length(G3_Mice),101);
for miceiter=G3_Mice
    for dayiter=G3_Days
        if(dayiter>length(Mice(miceiter).ExpDay))
            continue
        else
            G3Usage_split(miceiter,:) = G3Usage_split(miceiter,:) + ...
                Mice(miceiter).ExpDay(dayiter).usage;
        end
    end
end
PG3Usage_split=G3Usage_split./sum(G3Usage_split,2);

GUsage_splitHN = [G3Usage_split; G12Usage_split];

%test = [GUsage_splitH(1:8,2) GUsage_splitN(1:8,2)];
%[p, tb1] = anova2(test,4);

%% ttest with Bonferroni correction
PG12Usage_split_crop = PG12Usage_split(:,G2vsG1Sortedusageindex_crop);

[h, p] = ttest2(PG12Usage_split_crop(G1_Mice,:), ...
                PG12Usage_split_crop(G2_Mice,:), ...
                'Vartype','unequal');

% simple bonferonni 
alpha = 0.05;
m     = length(p);
bonf  = p<(alpha/m);

[p_sorted, p_index] = sort(p);

p_index_adj = SyllablesX(G2vsG1Sortedusageindex_crop(p_index));

%% plot syllable enrichment (cropped) with errorbar
close all
X_cutoff = 0:length(G2vsG1Sortedusageindex_crop)-1;

mean_PG1Usage_split_crop = mean(PG12Usage_split_crop(G1_Mice,:),1);
% equivalent to PG1Usage(G2vsG1Sortedusageindex_crop)
std_PG1Usage_split_crop = std(PG12Usage_split_crop(G1_Mice,:),0,1);

mean_PG2Usage_split_crop = mean(PG12Usage_split_crop(G2_Mice,:),1);
std_PG2Usage_split_crop = std(PG12Usage_split_crop(G2_Mice,:),0,1);


Plot_UsageCompare_cutoff_errorbar=figure(7);
errorbar(X_cutoff,mean_PG1Usage_split_crop,std_PG1Usage_split_crop,'LineWidth',1.5)
hold on
errorbar(X_cutoff,mean_PG2Usage_split_crop,std_PG2Usage_split_crop,'LineWidth',1.5)
plot(X_cutoff,PG3Usage(G2vsG1Sortedusageindex_crop),'LineWidth',1,'Color','Black')

legend({'Contextual Novelty','Stimulus Novelty','Habituation'},'FontSize',fsize)
title(['Syllable Usage Comparison of Contextual/Stimulus Novely Mice (', setName, ...
    ') (Sorted by stimulus novelty enrichment)'],'FontSize',fsize,'Interpreter', 'none')
ylabel('Fractional usage','FontSize',fsize)
xlabel('Syllables','FontSize',fsize)
xticks(X_cutoff);
xticklabels(SyllablesX(G2vsG1Sortedusageindex_crop));
set(Plot_UsageCompare_cutoff_errorbar, 'Position', [46 353 1870 450])

%saveas(Plot_UsageCompare_cutoff_errorbar, 'UsageCompare_Crop_wError.tif')

%% random shuffling to assess significance
[m,n] = size(GUsage_splitN_crop);
rng('default')

shufflemat = zeros(m,n);
for permiter = 1:n
    shufflemat(:,permiter) = randperm(8)';
end
b = GUsage_splitN_crop(shufflemat);
diffMean = mean(b(G1_Mice,:))-mean(b(G2_Mice,:));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pause
close all
clearvars Plot_UsageCompare_sorted Plot_UsageCompare Plot_AccGeneralUsage Plot_GeneralUsage
save('GeneralUsage.mat')
clear



