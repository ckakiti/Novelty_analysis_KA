clear
clc

cd('/media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq')
load('MiceIndex.mat')

%% find mouse/day

%curr_MSid = '0bcee6ae-8f86-416c-81c5-a727c9268eee';
curr_MSid = 'e9e81cda-4b31-4f40-85ac-fbed1608bbbb';
%curr_MSid = '47538c1f-1566-4f29-9d60-03bc5839dbf7';
for miceiter=1:length(Mice)
    for dayiter=1:length(Mice(miceiter).ExpDay)
        compare_MSid = Mice(miceiter).ExpDay(dayiter).MSid;
        if(strcmp(curr_MSid,compare_MSid))
            disp('found')
            disp(['Mouse: ' num2str(miceiter) '; Day: ' num2str(dayiter)])
        end
    end
end
%% find MSid index
miceiter = 12;
day_iter = 5;

MSidindex=1;
for indexiter=1:size(MoSeqDataFrame.session_uuid,1)
    if strcmp(MoSeqDataFrame.session_uuid(indexiter,:),Mice(miceiter).ExpDay(day_iter).MSid)
        break
    end
    MSidindex=MSidindex+1;
    if MSidindex==size(MoSeqDataFrame.session_uuid,1)+1
        error('MSid not found');
    end
end
disp(MSidindex)