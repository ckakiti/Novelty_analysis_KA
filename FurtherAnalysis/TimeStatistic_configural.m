% in subfolder when looking for .mat files, prsumed there is an Arena_Obj_Pos.mat file
clear
clc

Config_NovAna % ppc = 42/6.3;

cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Configural_set/'])
% currSet = 'Configural_set';
% cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' currSet '/'])

start_min=1;
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

% smoothing window for SAP analysis
Swindow=30;

% maximum body length (pixels) for SAP, anything above counted as error
bodyLen_cutoff = 80;

% squareSize = 4; % must be a factor of both vidWidth and vidHeight
% roundTargets = squareSize:squareSize:max(video_xlen, video_ywid);

%%
for folderi=1:folderlen
    cd(foldernames{folderi});
    disp(foldernames{folderi})
    
%     cd('Analyzed_Data_4obj')
    cd('Analyzed_Data_10cm')

    subpath=cd;
    PathRoot=[subpath '/'];
    filelist=dir([PathRoot, '*Converted.mat']); 
        %'*rgb_Converted.mat']); %[currMouse '*.mat']]); %foldernames{folderi}, '*.mat']);
    flen = length(filelist);

    for filei = 1:flen
        filename = filelist(filei).name;
        load(filename)
%         load('Arena_Obj_Pos_4obj.mat')
        
        if(endframe>length(Labels))
            endframe = length(Labels);
        end
        
        LabelsCut=Labels(startframe:endframe,:);
        
        Time_distance(filei,1)=sum(LabelsCut(:,21));
        Time_distance(filei,2)=sum(LabelsCut(:,22));
        Time_distance(filei,3)=sum(LabelsCut(:,23));
        Time_distance(filei,4)=sum(LabelsCut(:,24));
        
        diffs = hypot(diff(LabelsCut(:,14)), diff(LabelsCut(:,15)));
        diffs(isnan(diffs))=[];
        totalDistCut(1,filei) = sum(diffs)/ppc;
        
        
        % Analyse body length and SAP
        smooth_headX = smoothdata( LabelsCut(:,14), 'rloess', Swindow);%%%%
        smooth_headY = smoothdata( LabelsCut(:,15), 'rloess', Swindow);%%%%
        smooth_tailX = smoothdata( LabelsCut(:,11), 'rloess', Swindow);
        smooth_tailY = smoothdata( LabelsCut(:,12), 'rloess', Swindow);
        curr_velX  = LabelsCut(:,25);
        curr_velY  = LabelsCut(:,26);

        bodyLen = sqrt( (smooth_tailX-smooth_headX).^2 + (smooth_tailY-smooth_headY).^2 ); % pixels
        bodyLen(bodyLen>bodyLen_cutoff)=0; % body lengths >80px are likely to be errors

        vel        = sqrt( curr_velX.^2 + curr_velY.^2 ); % head velocity, cm/s
        stretches  = find(bodyLen>(mean(bodyLen)+std(bodyLen)));
        slows      = find(vel<5);
        sap        = intersect(stretches,slows);
%         sap_in_rad = intersect(sap, find(LabelsCut(:,21)==1));
        
        diffSap      = diff(sap);
%         diffSapInRad = diff(sap_in_rad);
        sap_frame{filei,1} = sap(find(diffSap>1)+1);
%         sap_num      = sum(diffSap>1); % total SAPs (any location)
%         sap_num_in_rad = sum(diffSapInRad>1); % SAPs near object
%         sap_num_in_rad(filei) = sum(diffSapInRad>1)/Time_distance(filei); % (norm)
        sap_dist{filei,1}  = LabelsCut(diffSap>1,17);
        nose_dist_sap{filei,1} = smooth_headX(diffSap>1);
        nose_dist_sap{filei,2} = smooth_headY(diffSap>1);
        nose_dist_all{filei,1} = LabelsCut(:,14);
        nose_dist_all{filei,2} = LabelsCut(:,15);
        
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

%         clearvars -except filelist flen filei startframe endframe Time_distance totalDist ... runArea  ...
%             folderi folderlen foldernames write_pointer TableName Time_distance_all DisOrAng ...
%             ppc arena currMouse video_ywid video_xlen Swindow bodyLen_cutoff dist_of_sap ...
%             sap_num sap_num_in_rad sap_frame sap_dist nose_dist %roundTargets
    end
    
    dist_of_sap(folderi).name       = foldernames{folderi};
    dist_of_sap(folderi).sapFrame   = sap_frame;
    dist_of_sap(folderi).sapDist    = sap_dist;
    dist_of_sap(folderi).nosePosSAP = nose_dist_sap;
    dist_of_sap(folderi).nosePosAll = nose_dist_all;
    
    if(strcmp(foldernames{folderi}, 'Aldehyde'))
        Time_distance_all(2:12,:,folderi) = Time_distance;
        TotalDistCut_all(folderi,2:12)    = totalDistCut;
    else
        Time_distance_all(:,:,folderi) = Time_distance;
        TotalDistCut_all(folderi,:)    = totalDistCut;
    end
        
    clearvars Time_angle Time_distance Time_periphery totalDistCut ...
        fracArea sap_num sap_num_in_rad sap_frame sap_dist nose_dist
    
    cd ../..

end

for filei = 1:flen
    rownames{filei, 1} = ['Session' num2str(filei)];
end

Time_distance_all=Time_distance_all./(endframe-startframe);
% Table=table(rownames, ...
%     Time_distance_all(:,1), Time_distance_all(:,2), ...
%     Time_distance_all(:,3), Time_distance_all(:,4), ...
%     totalDistCut_all);%, runArea);
% Table.Properties.VariableNames = {'Session' 'Obj1' 'Obj2' 'Obj3' 'Obj4', 'totalDist'};%, 'runArea'};

%% saving
cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Configural_set/')
save('Configural_set_SAP_plus_dist_10min_smooth.mat', 'dist_of_sap')
save('Configural_set_timeDist_totalDistCut.mat', 'Time_distance_all', 'TotalDistCut_all')

% cd Analyzed_Data_4obj
% writetable(Table,'TimeStatistic_4obj.csv');
% clear all