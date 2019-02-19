% in subfolder when looking for .mat files, prsumed there is an Arena_Obj_Pos.mat file
% omit the first .mat file which is Arena_Obj_Pos.mat. 
% If there are other .mat files in subfolder will cause error

clear
clc

Config_NovAna

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/7day_preexposure_combine %CvsS_180831_DLC/
%cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/NewHope-ROTJ/')

start_min=0.005;
end_min=start_min+10;
startframe=start_min.*60.*fps;
endframe=end_min.*60.*fps;

folderpath = cd;
folderd = dir(folderpath);
isub = [folderd(:).isdir];
foldernames = {folderd(isub).name}'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foldernames(ismember(foldernames,{'.','..','incomplete'})) = []; %%%%%%%%%%%%%%%%%%%%%%%
folderlen=length(foldernames);

for folderi=1:folderlen
    cd(foldernames{folderi});
    cd Analyzed_Data/body %%%%%%%%%%%%%%%%

    subpath=cd;
    PathRoot=[subpath '/'];
    filelist=dir([foldernames{folderi}(4:end), '*.mat']); % 'session01*.mat']; PathRoot, '*.mat']); %%%%%%%%%%%%%%%%%%
    flen = length(filelist);

    for filei = 1:flen
        filename = filelist(filei).name;
        load(filename)
        Time_distance(filei)=sum(Labels(startframe:endframe,21));
        Time_angle(filei)=sum(Labels(startframe:endframe,23));
        clearvars -except filelist flen filei Time_angle Time_distance ...
            startframe endframe folderi folderlen foldernames write_pointer ...
            TableName Time_distance_all Time_angle_all DisOrAng
    end
    Time_distance_all(folderi,:)=Time_distance;
    Time_angle_all(folderi,:)=Time_angle;
    clearvars Time_angle Time_distance
    cd ..
    cd .. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd ..
    DisOrAng{folderi,1}='Distatnce';
    DisOrAng{folderlen+folderi,1}='Angle';


end

Time_distance_all=Time_distance_all./(endframe-startframe);
Time_angle_all=Time_angle_all./(endframe-startframe);
Table=table(cat(1,foldernames,foldernames),DisOrAng,[Time_distance_all;Time_angle_all]);
writetable(Table,'TimeStatistic_body.csv');
clear all