%% create extract .m file with each uuid, 
%       name, and novelty condition (.mat -> .m)

clear
clc
close all

setName = 'CvsS_180831_MoSeq';
whichDay = 3; % 3=N1

% Mice_Index_path = ['/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/Behavior/Standard setup/' ...
%    'CombineAnalysis/Just-in-case files/Dataset_20190723/MiceIndex.mat'];
Mice_Index_path = '/media/alex/DataDrive1/MoSeqData/CvsS_20180831_MoSeq/MiceIndex.mat';
load(Mice_Index_path)

Mice_uuid_N1 = struct('name', {Mice.name}, ...
                      'novelty', {Mice.novelty});
for mouseiter = 1:length(Mice)
    Mice_uuid_N1(mouseiter).uuid = Mice(mouseiter).ExpDay(whichDay).MSid;
end

Mice_uuid_N1_table = cell2table({Mice_uuid_N1.uuid}');
if(0)
    writetable(Mice_uuid_N1_table, 'MiceIndex_combine3_uuidN1.csv', ...
        'WriteVariableNames', 0)
end
