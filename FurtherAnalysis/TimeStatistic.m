% in subfolder when looking for .mat files, prsumed there is an Arena_Obj_Pos.mat file
% omit the first .mat file which is Arena_Obj_Pos.mat. 
% If there are other .mat files in subfolder will cause error

clear
close all
clc

Config_NovAna

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/7day_preexposure_combine
%cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/NewHope-ROTJ/')

start_min=0.005;
end_min=start_min+10;
startframe=start_min.*60.*fps;
endframe=end_min.*60.*fps;

folderpath = cd;
folderd = dir(folderpath);
isub = [folderd(:).isdir];
foldernames = {folderd(isub).name}'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foldernames(ismember(foldernames,{'.','..','incomplete','Trajectories_H2N1'})) = []; %%%%%%%%%%%%%%%%%%%%%%%
folderlen=length(foldernames);

for folderi=1:folderlen
    cd(foldernames{folderi});
    cd Analyzed_Data/
    load('Arena_Obj_Pos.mat','arena')
    cd ./body %%%%%%%%%%%%%%%%

    subpath=cd;
    PathRoot=[subpath '/'];
    filelist=dir(['*rgb.mat']); % 'session01*.mat']; PathRoot, '*.mat']); foldernames{folderi}(4:end) %%%%%%%%%%%%%%%%%%
    flen = length(filelist);

    for filei = 1:flen
        filename = filelist(filei).name;
        load(filename)
        
        Time_distance(filei) = sum(Labels(startframe:endframe,21));
        Time_angle(filei)    = sum(Labels(startframe:endframe,23));
        
        arenaCorners = [arena(filei,1) arena(filei,3) arena(filei,3) arena(filei,1) arena(filei,1); ...
                        arena(filei,2) arena(filei,2) arena(filei,4) arena(filei,4) arena(filei,2)];
        arenaCenter  = [(arena(filei,1)+arena(filei,3))/2 (arena(filei,2)+arena(filei,4))/2];
        
        scale = 190;
        L  = linspace(0,2*pi,5) + pi/4; % select angles and thus number of sides for polygon,
                                        %   rotate so edges are parallel to arena edges
        xv = scale*cos(L)' + arenaCenter(1); % x points
        yv = scale*sin(L)' + arenaCenter(2); % y points
        
        xq = Labels(startframe:endframe,14);
        yq = Labels(startframe:endframe,15);
        
        [in,on] = inpolygon(xq,yq,xv,yv); % detect points that are inside defined polygon region
        
        % determine fraction of points outside polygon region
        Time_periphery(filei) = numel(xq(~in))/numel(xq);
        
        clearvars -except filelist flen filei Time_angle Time_distance Time_periphery ...
            startframe endframe folderi folderlen foldernames write_pointer arena ...
            TableName Time_distance_all Time_angle_all Time_periphery_all DisOrAng Periph
    end
    
    Time_distance_all(folderi,:)=Time_distance;
    Time_angle_all(folderi,:)=Time_angle;
    Time_periphery_all(folderi,:)=Time_periphery;
    clearvars Time_angle Time_distance Time_periphery
    cd ..
    cd .. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd ..
    DisOrAng{folderi,1}='Distance';
    DisOrAng{folderlen+folderi,1}='Angle';
    Periph{folderi,1}='Periphery';

end

Time_distance_all = Time_distance_all./(endframe-startframe);
Time_angle_all    = Time_angle_all./(endframe-startframe);

Table  = table(cat(1,foldernames,foldernames),DisOrAng,[Time_distance_all;Time_angle_all]);
Table2 = table(cat(1,foldernames), Periph, Time_periphery_all);

%writetable(Table,'TimeStatistic_body.csv');
%writetable(Table2,'TimeStatistic_body_periph.csv');

clear all