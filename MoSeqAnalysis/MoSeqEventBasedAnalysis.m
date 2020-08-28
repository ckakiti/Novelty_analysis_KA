%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

cd /media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq

load('MoSeqDataFrame.mat')
Syllablebinedge=[-6,-0.5:1:99.5];

MouseSet = 'Dataset_20190723'; 
% Mice_Index_path='/Users/yuxie/Dropbox/YuXie/CvsS_180831/CvsS_180831_MoSeq/Mice_Index.m';
Mice_Index_path=['/media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq/MiceIndex.mat'];%' MouseSet '.m'];
load(Mice_Index_path);
%run(Mice_Index_path);

mksize=80;
fsize=20;

fps=30;
starttime=0;         % in minutes  Pending: in the future, this should depend on the starting frame of each experiments
timeseg=0:2:10;      % in minutes
timeseg=starttime+timeseg;
frameseg=timeseg.*60.*fps;
frameCutoff = 18000;

Analysis_Mice=1:length(Mice);
Analysis_Day=1:6;     
Plot_SingleDay=3;    % first novelty day
Plot_MultiDay=3:6;   % all novelty days
cmap=cool(2*length(Plot_MultiDay));

IntSyl=9; %28 %71 %9            % Interesting Syllable

detectCond = cat(1, Mice.novelty);
Cnum = find(detectCond=='C');
Snum = find(detectCond=='S');
cond2Color = [0.5 0.0 0.5];
cond1Color = [1.0 0.5 0.0];

CVSXTick={'Contextual' 'Stimulus'};
CVSX=1:2;            % number of novelty conditions

for segiter=1:length(timeseg)-1
    IntersessionXTick{segiter}=[num2str(timeseg(segiter)) '-' num2str(timeseg(segiter+1))];
end
IntersessionX=1:length(timeseg)-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for miceiter=Analysis_Mice

    for dayiter=Analysis_Day
        if(miceiter==12 && dayiter==6)
            disp('skip this mouse/day')
            continue
        end
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

%         frameCutoff = length(MoSeqDataFrame.labels{MSidindex}); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(frameCutoff)
        Labels=double(MoSeqDataFrame.labels{MSidindex}(1:frameCutoff));
        labellen=length(Labels);

        Mice(miceiter).ExpDay(dayiter).IntSIndex=find(Labels==IntSyl);
        Mice(miceiter).ExpDay(dayiter).syltime=histcounts(Mice(miceiter).ExpDay(dayiter).IntSIndex,frameseg)./fps;
        Mice(miceiter).ExpDay(dayiter).syltime_total=sum(Mice(miceiter).ExpDay(dayiter).syltime);
        Mice(miceiter).ExpDay(dayiter).syltime_per_minute=Mice(miceiter).ExpDay(dayiter).syltime./(timeseg(2:end)-timeseg(1:end-1));
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(['/media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq/'])%' MouseSet '_MoSeq'])

syltime_total=zeros(length(Mice),1,length(Mice(1).ExpDay));
syltime_per_minute=zeros(length(Mice),length(timeseg)-1,length(Mice(1).ExpDay));
for miceiter=Analysis_Mice
    for dayiter=Analysis_Day
        if(miceiter==12 && dayiter==6)
            disp('skip this mouse/day')
            continue
        end
        syltime_total(miceiter,:,dayiter)=Mice(miceiter).ExpDay(dayiter).syltime_total;
        syltime_per_minute(miceiter,:,dayiter)=Mice(miceiter).ExpDay(dayiter).syltime_per_minute;
    end
end

disp('section 1')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bar plot of total syllable usage (s) for each mouse, divided by novelty condition
close all

% CVSX=1:3;
% CVSXTick={'6OHDA1' '6OHDA2' 'saline'};
% Cnum1 = [2 8 9];
% Cnum2 = [3 5 6];

Plot_total_interaction_time=figure;

Plotdata_avg=syltime_total(:,:,Plot_SingleDay);
csavg=[mean(Plotdata_avg(Cnum,:)) mean(Plotdata_avg(Snum,:))];
csstd=[std(Plotdata_avg(Cnum,:)) std(Plotdata_avg(Snum,:))];
% csavg=[mean(Plotdata(Cnum1,:)) mean(Plotdata(Cnum2,:)) mean(Plotdata(Snum,:))];
% csstd=[std(Plotdata(Cnum1,:)) std(Plotdata(Cnum2,:)) std(Plotdata(Snum,:))];

b = bar(CVSX,csavg);
b.FaceColor  = 'Flat';
b.CData(1,:) = cond2Color;
% b.CData(2,:) = cond2Color;
% b.CData(3,:) = cond1Color;
b.CData(2,:) = cond1Color;

hold on
errorbar(CVSX,csavg,csstd,'.','LineWidth',2,'Color','black')
hold on

scatter(ones(length(Cnum),1),Plotdata_avg(Cnum,1),...
    mksize,'filled','d','MarkerFaceColor','k')
scatter(repmat(2, length(Snum), 1),Plotdata_avg(Snum,1),...
    mksize,'filled','d','MarkerFaceColor','k')
%for miceiter=1:length(Mice)
%    if Mice(miceiter).novelty == 'C'
%        scatter(1,Plotdata_avg(miceiter),mksize,'filled','d')
%        hold on 
%    elseif Mice(miceiter).novelty == 'S'
%        scatter(2,Plotdata_avg(miceiter),mksize,'filled','d')
%        hold on
%    else
%        error('Novelty class not defined');
%    end
%end

% legend(['C vs S','SD',{Mice.name}])
set(Plot_total_interaction_time, 'position', [0 0 500 800]);
set(gca,'FontSize',fsize)
% ylim([0 50])
title(['Average syllable ' num2str(IntSyl) ' usage'])%on N1 (' MouseSet ')'])
%title({'Average syllable usage:','cautious approach'})
ylabel(['Syllable Usage (s)'])
xticks(CVSX);
xticklabels(CVSXTick);
xtickangle(45)
ylim([0 60])
set(gca,'YTick',[0 15 30 45])

if(0)
%     saveas(Plot_total_interaction_time, [MouseSet '_MoSeq_avgUsage_syl', num2str(IntSyl), '.tif'])
    save('Dataset_20190723_eventBased_syl9avg_10min.mat')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot syllable usage across time bins within 1 day, divide by novelty condition
close all

Plot_interaction_time=figure;
set(Plot_interaction_time,'position',[680 495 710 455])

Plotdata_avg=syltime_per_minute(:,:,Plot_SingleDay);

intersessioncavg=mean(Plotdata_avg(Cnum,:));
intersessioncstd=std(Plotdata_avg(Cnum,:));

intersessionsavg=mean(Plotdata_avg(Snum,:));
intersessionsstd=std(Plotdata_avg(Snum,:));

errorbar(IntersessionX,intersessioncavg,intersessioncstd,'Color',cond2Color,'LineWidth',2)
hold on
errorbar(IntersessionX,intersessionsavg,intersessionsstd,'Color',cond1Color,'LineWidth',2)

% hold on
% plot(IntersessionX, Plotdata(Cnum,:), 'Color', cond2Color, 'LineWidth',2)
% plot(IntersessionX, Plotdata(Snum,:), 'Color', cond1Color, 'LineWidth',2)

% for miceiter=1:length(Mice)
%     if Mice(miceiter).novelty == 'C'
%         scatter(IntersessionX,Plotdata(miceiter,:),'Marker','o','MarkerFaceColor','b','LineWidth',1)
%         hold on
%     elseif Mice(miceiter).novelty == 'S'
%         scatter(IntersessionX,Plotdata(miceiter,:),'Marker','o','MarkerFaceColor','r','LineWidth',1)
%         hold on
%     else
%         error('Novelty class not defined');
%     end
% end
xlim([0.5 length(timeseg)-0.5]);
ylim([-1 7])

set(gca,'FontSize',fsize)
% title(['Syllable ' num2str(IntSyl) ' per minute N1 (' MouseSet ')'])
title(['Time course of syllable ' num2str(IntSyl) ' usage'])
% title({'Time course of syllable usage:','cautious approach'})
% legend(['Contextual','Stimulus',{Mice.name}])
% legend({Mice.name},'location','north')
xlabel('Time Segments (min)')
%ylabel(['Syllable ' num2str(IntSyl) ' Usage (s)'])
ylabel('Syllable Usage (s)')
xticks(IntersessionX);
xticklabels(IntersessionXTick);
set(gca,'YTick',[0 2 4 6 8 10])

if(0)
%     saveas(Plot_interaction_time, [MouseSet '_MoSeq_usageAcrossSession_syl', num2str(IntSyl), '.tif'])
%     saveas(Plot_interaction_time, [MouseSet '_MoSeq_usageAcrossSession_syl', num2str(IntSyl), '_err.tif'])
    save('Dataset_20190723_eventBased_syl9acrossDay_10min.mat')
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot average syllable usage per novelty condition across novelty days

Plot_interaction_time_multiday=figure;

Plotdata_bins=syltime_per_minute(:,:,Plot_MultiDay);

intersessioncavg=mean(Plotdata_bins(Cnum,:,:),1);
intersessioncstd=std(Plotdata_bins(Cnum,:,:),1);

intersessionsavg=mean(Plotdata_bins(Snum,:,:),1);
intersessionsstd=std(Plotdata_bins(Snum,:,:),1);

for dayiter=1:length(Plot_MultiDay)
    plot(IntersessionX,intersessioncavg(:,:,dayiter),'LineWidth',1.5,'Color',cmap(dayiter,:))
    hold on
end
for dayiter=1:length(Plot_MultiDay)
    plot(IntersessionX,intersessionsavg(:,:,dayiter),'LineWidth',1.5,'Color',cmap(end-dayiter,:))
    hold on
end
legend('C-N1','C-N2','C-N3','C-N4','S-N1','S-N2','S-N3','S-N4');
xlim([0.5 length(timeseg)-0.5]);
xticks(IntersessionX);
xticklabels(IntersessionXTick);
xlabel('Time Segments (min)','FontSize',fsize)
ylabel(['Syllable ' num2str(IntSyl) ' Usage Time (s)'],'FontSize',fsize)
title(['Syllable ' num2str(IntSyl) ' per min on all novelty days (Capoeira)'],'FontSize',fsize)

% saveas(Plot_interaction_time_multiday, [MouseSet '_MoSeq_usageAcrossDays_syl', num2str(IntSyl), '.tif'])

%% statistical analysis

% for avg syl expression; non-parametric 2-sample test
[curr_p, curr_h, curr_stats] = ...
        ranksum(Plotdata_avg(Snum,1), Plotdata_avg(Cnum,1))

% for plot across days; non-parametric 2-sample test with bonf-holm correction
clear stat_summary curr_p curr_h curr_stats
for biniter = 1:size(Plotdata_bins,2)
    [curr_p, curr_h, curr_stats] = ...
        ranksum(Plotdata_bins(Snum,biniter), Plotdata_bins(Cnum,biniter));
    stat_summary(biniter).p = curr_p;
    stat_summary(biniter).h = curr_h;
end
[corrected_p, corrected_h] = bonf_holm(cat(1,stat_summary.p),0.05)

