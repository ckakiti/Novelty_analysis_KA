clear
clc
close all

cd /media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq

load('MoSeqDataFrame.mat')
Syllablebinedge=[-6,-0.5:1:99.5];

Mice_Index_path='/media/alex/DataDrive1/MoSeqData/Dataset_20190723/MiceIndex/MiceIndex_Capoeira.m';
run(Mice_Index_path);

MarkerSize=5;
Fsize=20;


%%
% IntSyllable=[45 20 9 33 14 17 64 62 60 95 85 94 36 50 55 11 0 10 16 2];%[42 94 39 38 70 13];
% IntSyllable=[36 37 74 18 24];
IntSyllable = [15 41 72 1 4 23 63 75 31 52 44 33 49 47 8];

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
saveas(SyllableDis, 'SyllablePosDis.tif')
