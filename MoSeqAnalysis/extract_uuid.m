%% DO THIS BEFORE RUNNING CODE BELOW
% automated way to extract all uuids from moseq2-index.yaml file

% open moseq2-index.yaml in Visual Studio
% CTRL F, enter this into box (path:\r?$)|(.*00.h5)|(.*yaml)|(\s uuid.*)
% CTRL + SHIFT + L to select all occurrances found
% CTRL + C --> CTRL + V into new file and name it moseq2-index_extract
% run below code
%%
clear
close all
clc

cd /media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq
run('MiceIndex_blank')

% identify basic information (name, date, uuid) within file
extract_file = importdata('moseq2-index_extract');
extract_reshape = reshape(extract_file, 4, [])';
extract_reshape = extract_reshape(:,[2 4]);

extract_separate = cell(size(extract_reshape,1),3);
for extractiter = 1:length(extract_reshape)
    curr_string = extract_reshape{extractiter,1};
    session_pos = strfind(curr_string,'session');
    
    % date
    curr_date = curr_string(session_pos-7:session_pos-2);
    extract_separate(extractiter,2) = {curr_date};
    
    % mouse name
    where_slash = strfind(curr_string, '/');
    
    % change curr_name depending on path format
    curr_name = curr_string(where_slash(1)+1:where_slash(2)-1); %%%%%%%%%%%%
%     curr_name = curr_string(5:where_slash(1)-1);                %%%%%%%%%%%%
    
    extract_separate(extractiter,1) = {curr_name};
    
    % uuid string
    uuid_curr = extract_reshape{extractiter,2};
    uuid_crop = uuid_curr(9:end);
    extract_separate(extractiter,3) = {uuid_crop};
end

% integrate basic information into MiceIndex
for mouseiter = 1:length(Mice)
    currMouse = Mice(mouseiter).name;
    
    whereMouse = cellfun(@(s) contains(currMouse,s), extract_separate(:,1));
    whereMouse = find(whereMouse);
    curr_days  = cat(1,extract_separate{whereMouse,2});
    curr_uuids = cat(1,extract_separate{whereMouse,3});
    days_cell  = cellstr(curr_days);
        
    % sort days
    translate_day = datetime(curr_days, 'InputFormat', 'yyyyMMdd');
    sorted_days   = datestr(sort(translate_day),'yymmdd');
    
    for daysiter = 1:length(whereMouse)
        % add date
        Mice(mouseiter).ExpDay(daysiter).date = sorted_days(daysiter,:);
        
        % find index of date and add corresponding uuid
        curr_idx = find(strcmp(sorted_days(daysiter,:),days_cell));
        Mice(mouseiter).ExpDay(daysiter).MSid = curr_uuids(curr_idx,:);
    end
end

if(0)
    save('MiceIndex.mat', 'Mice')
end

