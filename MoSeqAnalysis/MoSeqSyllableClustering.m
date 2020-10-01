clear
clc
close all
%%
cd /media/alex/DataDrive1/MoSeqData/CvsS_20180831_MoSeq

%load('MoSeqDataFrame_syllabledissorted.mat');
load('MoSeqDataFrame_curr.mat')
%modelLabels = csvread('CvsS_modelLabel_sortedByUsage.txt');

for syllableiter=1:100
    nodename{syllableiter}=num2str(syllableiter-1);
    %nodename{syllableiter}=num2str(modelLabels(syllableiter));
end

CThreshold=0.7;
NodeNum=100;
LWidth=2;
FSize=20;

syllable_dis_vct=squareform(MoSeqDataFrame.syllable_dis);
syllable_linkage=linkage(syllable_dis_vct,'average');
leafOrder = optimalleaforder(syllable_linkage,syllable_dis_vct,'Criteria','group');

%%
clc
close all

Plot_SyllableDis=figure;
set(Plot_SyllableDis, 'Position', [180 430 1500 450])

[SyllableDenG,T,outperm]=dendrogram(syllable_linkage,NodeNum,...
    'reorder',leafOrder,'ColorThreshold',...
    CThreshold*max(syllable_linkage(:,3)),'Labels',nodename);
title('Syllable Distance Clustering','FontSize',FSize);
xlabel('Syllable Number Rank (Sorted by usage)','FontSize',FSize);
ylabel('Syllable Distance','FontSize',FSize)
set(SyllableDenG,'LineWidth',LWidth);

if(0)
    saveas(Plot_SyllableDis, 'CvsS_syllableClustering_test.tif')
end
