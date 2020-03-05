%% extract all uuids from moseq2-index.yaml file
clear
close all
clc

cd /media/alex/DataDrive1/MoSeqData/Planets/Planets_MoSeq
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
%     curr_name = curr_string(where_slash(1)+1:where_slash(2)-1); 
    curr_name = curr_string(5:where_slash(1)-1);
    
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

