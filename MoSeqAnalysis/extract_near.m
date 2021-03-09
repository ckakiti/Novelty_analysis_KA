%% extract when mouse is within r radius of obj, integrate into MiceIndex.mat
%  (later, convert rgb_ts to depth_ts with convert_rgbts_depthts.m)
clear
clc

% set_names = {'Capoeira_DLC', 'Chess_DLC', 'Hiking_DLC'};
set_names = {'CvsS_180831_DLC'};

% load MiceIndex file (contains mice, condition, date, MSid) 
% cd /media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq
cd /media/alex/DataDrive1/MoSeqData/CvsS_20180831_MoSeq
load('MiceIndex.mat') % generated from extract_uuid.m

for setiter =1:length(set_names)
    cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/'...
        set_names{setiter}])
    
    files    = dir;
    whichDir = [files.isdir];
    nameDir  = files(whichDir);
    nameDir  = {nameDir.name};
    nameDir(ismember(nameDir,{'.','..','temp'})) = [];
    if(length(nameDir)~=6)
        disp('Incorrect number of mice in folder (~=6)')
    end

    for mouseiter = 1:length(nameDir)
        cd(nameDir{mouseiter})
%         cd Analyzed_Data_1obj_head/
        cd Analyzed_Data/nose_8cm %%%
        session_files = dir('*rgb_Converted.mat');
        session_names = cat(1,session_files.name);
        
        for sessioniter = 1:size(session_names,1)
            clear Labels isIn_curr
            load(session_names(sessioniter,:))
            
            isIn_curr = find(Labels(:,21));
            
            currMousePos = find(strcmp({Mice.name}, ...
                                        nameDir{mouseiter}));
            find_date = strfind(session_names(sessioniter,:),'_');
            currDate = session_names(sessioniter,find_date(1)+1:find_date(2)-1);
            currSessionPos = find(strcmp({Mice(currMousePos).ExpDay.date},currDate));
            
            Mice(currMousePos).ExpDay(currSessionPos).isIn_rgb = isIn_curr;
        end
        
        cd ../..
        cd .. %%%
    end
end
disp('end')

if(0)
    cd /media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq/
    save('MiceIndex_isIn_head.mat', 'Mice')
end