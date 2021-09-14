%% combine information from cumulative orientation mat files

clear
close all
clc

path_to_index_files = ['/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/' ...
    'Behavior/others/Standard setup/CombineAnalysis/Dataset_20190723'];
cd(path_to_index_files)
run('MiceIndex_combine3')

path_to_mat_files = ['/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/' ...
    'Behavior/others/Standard setup/CombineAnalysis/Just-in-case files'];
cd(path_to_mat_files)
set_names = {'Capoeira_DLC','Chess_DLC','Hiking_DLC'};

orientCumul_all = [];
for seti = 1:length(set_names)
    cd(set_names{seti})
    
    load('TimeStatistic_orient_cumul.mat') % generated in TimeStatistic_trim.m
    orientCumul_all = [orientCumul_all; Time_angle_cumul];
    
    cd ..
end

if(0)
    cd(path_to_index_files)
    save('Dataset_20190723_orientCumul.mat', 'orientCumul_all')
    
    disp('saved')
end