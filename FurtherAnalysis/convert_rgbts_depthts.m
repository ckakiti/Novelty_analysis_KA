%% convert pokes from rgb ts to depth ts frame (from bout_analysis.m)
clear
close all
clc

cd('/media/alex/DataDrive1/MoSeqData/CvsS_20180831_MoSeq/')
% poke_labels=csvread('Alcohol_181130_01_poke_manual_test.csv',1);
%'Capoeira_poke_labels_N1.csv',1,3);
load('MiceIndex_isIn_nose_8cm.mat')

% cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Machines_DLC/')
% load('PokesApproaches.mat')

% cd('/media/alex/DataDrive1/MoSeqData/Superset/')

% curr_labels = [];
for mousei=1:length(Mice)
%     curr_rgb = Mice(mousei).Pokes_Day3'; % identify pokes for N1
%     curr_rgb = Mice(mousei).Approach_Day3';
%     curr_rgb = Mice(mousei).ExpDay; %%%%%%%%%%%%%%%%%%%%%
%     curr_rgb = poke_labels;
    
    cd(Mice(mousei).name)
    
    % if using conf novelty file structure (1 day, 12 sessions/day)
%     cd 181130  % day
%     files = dir('session*');
%     files = {files.name}';
%     cd(files{1}) % session 1
     
    files    = dir;
    whichDir = [files.isdir];
    dayDir  = files(whichDir);
    dayDir  = {dayDir.name};
    dayDir(ismember(dayDir,{'.','..','temp'})) = [];

    for dayiter = 1:length(dayDir)
        cd(dayDir{dayiter})
        curr_session = dir('session*');
        cd(curr_session.name)
        
        rgbts_file   = dir('*rgb_ts.txt');
        depthts_file = dir('*depth_ts.txt');
        % Capoeira: 190413_01
        % Hiking: 190622
        % Chess: 190712
        % Machines: 190906
        % Planets: 200102

        rgbts   = load(rgbts_file.name);
        depthts = load(depthts_file.name);
        
        curr_rgb = Mice(mousei).ExpDay(dayiter).isIn_rgb;
        curr_depth = zeros(length(curr_rgb),1);
        for ts_iter = 1:length(curr_rgb)
            if(curr_rgb(ts_iter)>length(rgbts))
                curr_rgb(ts_iter) = [];
                curr_depth(ts_iter) = [];
            elseif(rgbts(curr_rgb(ts_iter),1)>depthts(end,1))
                curr_depth(ts_iter,1) = length(depthts);
            else
                curr_depth(ts_iter,1) = find(depthts>rgbts(curr_rgb(ts_iter),1),1);
            end
        end
        
        Mice(mousei).ExpDay(dayiter).isIn_depth = curr_depth;
%         curr_AddLabel  = [repmat(mousei,length(curr_rgb),1) curr_rgb curr_depth];
%         curr_labels    = [curr_labels; curr_AddLabel];
        
        cd ../..
    end
    
    cd .. 
end
disp('end')

if(0)
    % csvwrite('***_poke_labels_N1.csv', curr_labels)
    
    cd /media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq/
    save('MiceIndex_isIn_head_depth.mat', 'Mice')
end

%% create csv file indicating which frames mouse is near obj (depthts)
clear
clc

cd /media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq/
load('MiceIndex_isIn_head_depth.mat')

isIn_MSid  = [];
isIn_depth = [];
for mouseiter = 1:length(Mice)
    currSessions = Mice(mouseiter).ExpDay;
    
    for dayiter = 1:length(currSessions)
        currMSid  = currSessions(dayiter).MSid;
        currDepth = currSessions(dayiter).isIn_depth;
        
        isIn_MSid  = [isIn_MSid; repmat(currMSid,length(currDepth),1)];
        isIn_depth = [isIn_depth; currDepth];
        
        if(mouseiter==2 & dayiter==4)
            %pause
            disp(length(isIn_depth)-length(currDepth))
            disp(length(isIn_depth))
        end
    end
end
isIn_MSid = string(isIn_MSid);

if(0)
    csvwrite('***_isIn_depth.csv', isIn_depth)
    writetable('***_isIn_MSid.txt', isIn_MSid_table)
    
    fid = fopen('***_isIn_MSid.txt','w');
    fprintf(fid, '%s\n', isIn_MSid);
    fclose(fid);
end
