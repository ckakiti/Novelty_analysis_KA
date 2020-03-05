% in subfolder when looking for .mat files, prsumed there is an Arena_Obj_Pos.mat file
clear
clc

Config_NovAna % ppc = 42/6.3;

currMouse = 'Queixada';
cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Capoeira_DLC/' currMouse '/'])

start_min=(10/60);
end_min=start_min+10;
fps=15; %%%%%%%%%%%%%%%%%%%%%%
startframe=start_min.*60.*fps;
endframe=end_min.*60.*fps;

folderpath = cd;
folderd = dir(folderpath);
isub = [folderd(:).isdir];
foldernames = {folderd(isub).name}';
foldernames(ismember(foldernames,{'.','..',...
    'Analyzed_Data','Analyzed_Data_1obj','RetrainSep17','temp'})) = [];
folderlen=length(foldernames);

% squareSize = 4; % must be a factor of both vidWidth and vidHeight
% roundTargets = squareSize:squareSize:max(video_xlen, video_ywid);

% for folderi=1:folderlen
%     cd(foldernames{folderi});
    cd('Analyzed_Data_4obj')

    subpath=cd;
    PathRoot=[subpath '/'];
    filelist=dir([PathRoot, '*rgb_Converted.mat']); %[currMouse '*.mat']]); %foldernames{folderi}, '*.mat']);
    flen = length(filelist);

    for filei = 1:flen
        filename = filelist(filei).name;
        load(filename)
%         load('Arena_Obj_Pos_4obj.mat')
        
        if(endframe>length(Labels))
            endframe = length(Labels);
        end
        
        Time_distance(filei,1)=sum(Labels(startframe:endframe,21));
        Time_distance(filei,2)=sum(Labels(startframe:endframe,22));
        Time_distance(filei,3)=sum(Labels(startframe:endframe,23));
        Time_distance(filei,4)=sum(Labels(startframe:endframe,24));
        
        diffs = hypot(diff(Labels(:,14)), diff(Labels(:,15)));
        totalDist(filei,1) = sum(diffs)/ppc;
        
%         roundedPawnArea  = interp1(roundTargets, roundTargets, [Labels(:,14) Labels(:,15)], 'nearest');
%         
%         if(sum(isnan(roundedArea(:)))>0)
%             for replace_nan = 1:length(roundedArea)
%                 if(isnan(roundedArea(replace_nan,1)))
%                     roundedArea(replace_nan,1) = roundedArea(replace_nan-1,1);
%                 end
%                 if(isnan(roundedArea(replace_nan,2)))
%                     roundedArea(replace_nan,2) = roundedArea(replace_nan-1,2);
%                 end
%             end
%         end
%         
%         linearInd    = sub2ind([video_xlen video_ywid], roundedArea(:,1), roundedArea(:,2)); %%%%%%
%         arenaTargets = interp1(roundTargets, roundTargets, arena(filei,:), 'nearest');
%         arenaArea    = (arenaTargets(1,4)-arenaTargets(1,2))*...
%                        (arenaTargets(1,3)-arenaTargets(1,1))/16;
%         runArea(filei,1) = length(unique(linearInd))/arenaArea;
        
        clearvars -except filelist flen filei startframe endframe Time_distance totalDist ... runArea  ...
            folderi folderlen foldernames write_pointer TableName Time_distance_all DisOrAng ...
            ppc arena currMouse video_ywid video_xlen %roundTargets
    end
%     Time_distance_all(folderi,:)=Time_distance;
    cd ..

% end

for filei = 1:flen
    rownames{filei, 1} = ['Session' num2str(filei)];
end

Time_distance=Time_distance./(endframe-startframe);
Table=table(rownames, ...
    Time_distance(:,1), Time_distance(:,2), Time_distance(:,3), Time_distance(:,4), ...
    totalDist);%, runArea);
Table.Properties.VariableNames = {'Session' 'Obj1' 'Obj2' 'Obj3' 'Obj4', 'totalDist'};%, 'runArea'};

%% saving
% cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Rim_KO_DLC/' currMouse])
cd Analyzed_Data_4obj
writetable(Table,'TimeStatistic_4obj.csv');
clear all