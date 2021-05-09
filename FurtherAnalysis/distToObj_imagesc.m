%% imagesc plot for distToObj (combine groups - based on Mitsuko code)
%    index file generated from analy_novelty_multi_sessions_2012_ka.m

clear
close all
clc
pause(0.5)

group = {'StandardSetup_combine','StandardLesion_combine'};%,'Iku_6OHDA_DLC'};
folder_base = '/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/';

% load index files for each group
Mice_combine = [];
for group_iter = 1:length(group)
    cd(strcat(folder_base, group{group_iter}));
    miceIndex_matFile = dir('MiceIndex.mat');
    miceIndex_name = miceIndex_matFile.name;
    
    load(miceIndex_name)
    Mice_combine = [Mice_combine Mice];
end

detectCond = cat(1, Mice_combine.novelty);
cond2 = find((detectCond=='C')|(detectCond=='l'));
cond1 = find((detectCond=='S')|(detectCond=='s'));
which_cond = cond1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% standardize lengths of sessions and # of sessions
min_session_num     = min(cellfun(@length,{Mice_combine.noseLog}));
start_frame         = 30; %1 %30
crop_session_length = start_frame:(25*15*60+start_frame-1);
for miceiter = 1:length(Mice_combine)
    % crop number of sessions to be equal across mice
    Mice_combine(miceiter).noseLog = Mice_combine(miceiter).noseLog(1:min_session_num);
    Mice_combine(miceiter).tailLog = Mice_combine(miceiter).tailLog(1:min_session_num);
    
    % crop length of each session
    uncroppedNoseLog = squeeze(struct2cell(Mice_combine(miceiter).noseLog));
    uncroppedTailLog = squeeze(struct2cell(Mice_combine(miceiter).tailLog));
    croppedNoseLog = cellfun(@(s) s(crop_session_length), uncroppedNoseLog, 'UniformOutput',false);
    croppedTailLog = cellfun(@(s) s(crop_session_length), uncroppedTailLog, 'UniformOutput',false);
    
    Mice_combine(miceiter).noseLog = croppedNoseLog;
    Mice_combine(miceiter).tailLog = croppedTailLog;
end

%% imagesc of nose distToObj
%  (plot stimulus novelty and saline mice together)
close all
clc

noseLog_combine = cell2mat([Mice_combine(which_cond).noseLog])';
tailLog_combine = cell2mat([Mice_combine(which_cond).tailLog])';

which_stat = 'tail'; % nose %tail
if(strcmp(which_stat,'nose'))
    statLog_combine = cell2mat([Mice_combine(which_cond).noseLog])';
elseif(strcmp(which_stat,'tail'))
    statLog_combine = cell2mat([Mice_combine(which_cond).tailLog])';
end

fig1 = figure(1);
set(fig1, 'Position', [170 450 1500 450])

imagesc(-statLog_combine)
% imagesc(noseLog_combine-tailLog_combine)
caxis([-log10(60) 1])
colormap;
clrbar = colorbar('Ticks',[-log10(60) -1 0 1],...
                  'TickLabels',{'60','10','1','0.1'});

% imagesc(10.^(statLog_combine))
% caxis([1 60])
% clrmap = colormap;
% colormap(flipud(clrmap))
% clrbar = colorbar('Direction','reverse');

% yticks      [-1     0    1   log10(60)])
% yticklabels {'0.1','1','10','60'})

xlabel('Time (min)')
ylabel('Mouse')
title(['Distance from object (cm): ' which_stat])
x_ticks = length(crop_session_length)*[1 2 3 4 5];
set(gca,'FontSize',16,'tickdir','out',...'YDir','normal',
    'XTick',x_ticks,...
    'XTickLabel',strsplit(num2str(x_ticks/15/60)))
    %, 'TickLength',2*(get(gca,'TickLength')))

if(0)
    cd(strcat(folder_base, 'StandardSetup_stimSaline'))
    saveas(fig1, ['distToObj_imagesc_' which_stat '.tif'])
end