% Analyze bouts of interaction with object

clear
close all
clc

Config_NovAna
cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Planets_DLC/')
run('MiceIndex')

% only analyze first 10 minutes of video
start_min=0;%0.5;
end_min=start_min+10.25;
% startframe=start_min.*60.*fps;
% endframe=end_min.*60.*fps;
% startframe=round(startframe);
% endframe=round(endframe);
startframe=Dis_ts_frame;
endframe=Dis_te_frame;

% for smoothing
Swindow=40;
Xtime=startframe:endframe;

files    = dir;
whichDir = [files.isdir];
nameDir  = files(whichDir);
nameDir  = {nameDir.name};
nameDir(ismember(nameDir,{'.','..','Retrain_Sep17','temp'})) = [];
days = 6;

shift_time = input('Exclude first 5 min for some mice? 0/1: ');
if(shift_time==1)
    whichFiles_shift = [4 5 6];
    disp(['Mice where lego is placed in arena at 5min: ' num2str(whichFiles_shift)])
end

%% number of bouts and bout length
boutNum = zeros(length(nameDir),days);
boutLen = zeros(length(nameDir),days);

for mousei = 1:length(nameDir)
    cd([nameDir{1,mousei}, '/Analyzed_Data_1obj_12cm_tail'])
    load('Arena_Obj_Pos.mat')
    
    matFiles = dir('*rgb_Converted.mat'); %'*rgb.mat'); %'*Converted.mat');
    matFiles = cat(1, matFiles.name);
    
    
    % exclude first 5 min of recording (if lego is placed in arena after 5 min)
    if(shift_time==1 && any(mousei==whichFiles_shift))
        timeShift = 10100; %11000-900
        disp(timeShift)
        
        startframe = Dis_ts_frame+timeShift;
        endframe   = Dis_te_frame+timeShift;
        Xtime      = startframe:endframe;
        
        disp(nameDir{1,mousei})
        disp(startframe)
    else
        startframe = Dis_ts_frame;
        endframe   = Dis_te_frame;
        Xtime      = startframe:endframe;
    end
    
    
    for filei = 1:size(matFiles,1)     %%%% which day to analyze %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load(matFiles(filei,:))

        if(startframe<=0)
            startframe=1;
            Xtime=startframe:endframe;
        end
        if(endframe>size(Labels,1))
            endframe = size(Labels,1);
            Xtime=startframe:endframe;
        end
        
        DisX=Labels(Xtime,2); % pix, nose
        DisY=Labels(Xtime,3); % pix, nose
        Dis=sqrt((obj_center(filei,1)-DisX).^2+(obj_center(filei,2)-DisY).^2)/ppc; % cm

        SDis=smoothdata(Dis,'rloess',Swindow);
        SDis_inRad = double(SDis<radius_cm);
        
        DisMin=islocalmin(SDis);
        DisMax=islocalmax(SDis);
        findDisMax = find(DisMax);
        
        
        % find frames where mouse crosses within radius specified by Config_NovAna
%         crossTmp  = crossing(Labels(startframe:endframe,21), [], 0.5);
        crossTmp  = crossing(SDis_inRad, [], 0.5);
        cross_In  = crossTmp(1:2:end);
        cross_Out = crossTmp(2:2:end);
        
        % make sure pairs of crossIn/crossOut are correct
        %  (make sure first crossing point is indeed crossing in, make sure
        %   there's an even number of crossings)
        vel_SDis = diff(SDis);
        
        if(vel_SDis(cross_In(1))>0)
            disp('incorrect cross in')
            cross_In(1) = [];
            cross_Out_new = cross_In;
            cross_In_new  = cross_Out;
            cross_Out     = cross_Out_new;
            cross_In      = cross_In_new;
        end
        
        if(length(cross_Out)>length(cross_In))
            disp('cross out > cross in')
            cross_Out(end) = [];
        elseif(length(cross_In)>length(cross_Out))
            disp('cross in > cross out')
            cross_In(end) = [];
        end
        
        % find local minimum within each bout (aka pokes)
        pokesTemp=[];
        approachTemp=[];
        for crossIter = 1:length(cross_In)
            localMinTemp = find(DisMin(cross_In(crossIter):cross_Out(crossIter)));
            if(length(localMinTemp)==1)
                pokesTemp = [pokesTemp cross_In(crossIter)+localMinTemp];
            else
                [minTemp, whichMin] = min(SDis(cross_In(crossIter)+localMinTemp));
                pokesTemp = [pokesTemp cross_In(crossIter)+localMinTemp(whichMin)];
            end
            
            findApproach = find(DisMax(1:cross_In(crossIter)),1,'last');
            approachTemp = [approachTemp findApproach];
        end
        
        Mice(mousei).(['Pokes_Day' num2str(filei)])=pokesTemp+startframe; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mice(mousei).(['Approach_Day' num2str(filei)])=approachTemp+startframe; %%%%%%%%%%%%%%%%%%%%%%%
        
        if(~isempty(cross_In))
            boutNum(mousei,filei) = length(cross_In);
            boutLen(mousei,filei) = mean(cross_Out-cross_In)/fps; % seconds
        end
    end
    
    cd ../..
% end
disp('end')

% % sanity check plotting
if(0)
    load('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/DRILLS_DLC/DRILLS_poke_labels_N1.csv')
    whichLabels = find(DRILLS_poke_labels_N1(:,1)==mousei);
    poke_labels = DRILLS_poke_labels_N1(whichLabels,2);
    
    figure;
    hold on
    plot(obj_center(mousei,1), obj_center(mousei,2), 'r*')
    for int = 1:length(poke_labels)
        plot(DisX(poke_labels(int)), DisY(poke_labels(int)), 'k.')
    end
    
    
    %
    cd(Mice(mousei).name)
    cd('Analyzed_Data_1obj_10min')
    load('Arena_Obj_Pos.mat')
    close all
    
    boutPlot = figure(1);
    set(gcf, 'Position', [46 351 1872 450])
    plot(SDis)
    hold on
    line([0 endframe], [radius_cm radius_cm])
    scatter(cross_In, repmat(radius_cm, 1, length(cross_In)))
    scatter(cross_Out, repmat(radius_cm, 1, length(cross_Out)))
    scatter(pokesTemp, SDis(pokesTemp))
    scatter(approachTemp, SDis(approachTemp))
    ylim([0 60])
    title([Mice(mousei).name ' Day ' num2str(filei)])
%     scatter(Mice(1).Pokes_Day3, SDis(Mice(1).Pokes_Day3))
%     scatter(Mice(1).Approach_Day3, SDis(Mice(1).Approach_Day3))
%     scatter(find(DisMin),SDis(DisMin))


    correctPokingLabels = Mice(mousei).(['Approach_Day' num2str(filei)])...
        <Mice(mousei).(['Pokes_Day' num2str(filei)]);
    currPokes = Mice(mousei).(['Pokes_Day' num2str(filei)])(correctPokingLabels);
    currApproach = Mice(mousei).(['Approach_Day' num2str(filei)])(correctPokingLabels);

    boutTracePlot = figure(2);
    hold on
    for int = 1:length(currPokes)
        currDisX = DisX(currApproach(int):currPokes(int));
        currDisY = DisY(currApproach(int):currPokes(int));
        plot(currDisX, currDisY)%, 'k')
    end
    plot(obj_center(1,1), obj_center(1,2), 'r*')
    
    th = 0:pi/50:2*pi;
    x  = obj_center(1,1);
    y  = obj_center(1,2);
    xunit = radius_cm * ppc * cos(th) + x; %radius
    yunit = radius_cm * ppc * sin(th) + y; %radius
    plot(xunit, yunit)
    xlim([0 video_xlen])
    ylim([0 video_ywid])
    
    axis square
    set(gca, 'YDir', 'reverse')
    title(['Threshold: ', num2str(radius_cm), 'cm'])
        
    pause
    
    saveas(boutPlot,[Mice(mousei).name '_Day' num2str(filei) '_bout_analysis'])
    saveas(boutTracePlot,[Mice(mousei).name '_Day' num2str(filei) '_bout_traces'])
    
%     csvwrite('boutStart',cross_In)
%     csvwrite('boutEnd',cross_Out)
%     csvwrite('NoveltyResponse_Queens_N1_poke', pokesTemp+startframe-1) %%%%%%%%%% change name %%%%%%%%%%%%%%
%     csvwrite('NoveltyResponse_Queens_N1_boutStart', approachTemp+startframe-1) %% change name %%%%%%%%%%%%%%
    
    % figure(2)
    % hist(cross_Out-cross_In,100)
    cd ../..
end

end
close all

% saving
numOrLen = cell(2*length(nameDir),1);
numOrLen(1:length(nameDir)) = {'boutNum'};
numOrLen(length(nameDir)+1:end) = {'boutLen'};
TableLabels = [cat(2,nameDir,nameDir)', numOrLen];
TableInfo   = [boutNum;boutLen];
Table       = [cell2table(TableLabels) array2table(TableInfo)];
disp('end')
if(0)
    writetable(Table,'boutAnalysis_tail.csv');
    save('PokesApproaches', 'Mice')
end
%clear all

%% convert pokes from rgb ts to depth ts frame
clear
close all
clc

% cd('/media/alex/DataDrive1/MoSeqData/Dataset_20190723/MoSeq/Capoeira_MoSeq')
% poke_labels=csvread('Capoeira_poke_labels_N1.csv',1,3);

cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Machines_DLC/')
load('PokesApproaches.mat')

files    = dir;
whichDir = [files.isdir];
nameDir  = files(whichDir);
nameDir  = {nameDir.name};
nameDir(ismember(nameDir,{'.','..','temp'})) = [];

curr_labels = [];
for mousei=1:length(Mice)
%     curr_rgb = Mice(mousei).Pokes_Day3'; % identify pokes for N1
    curr_rgb = Mice(mousei).Approach_Day3';
    
    cd(nameDir{mousei})
    
    rgbts_file   = dir('*190906_rgb_ts.txt');
    depthts_file = dir('*190906_depth_ts.txt');
    % Capoeira: 190413_01
    % Hiking: 190622
    % Chess: 190712
    % Machines: 190906
    % Planets: 200102
    
    rgbts   = load(rgbts_file.name);
    depthts = load(depthts_file.name);
    
    curr_depth = zeros(length(curr_rgb),1);
    for ts_iter = 1:length(curr_rgb)
        curr_depth(ts_iter,1) = find(depthts>rgbts(curr_rgb(ts_iter)),1);
    end
    
    curr_AddLabel  = [repmat(mousei,length(curr_rgb),1) curr_rgb curr_depth];
    curr_labels    = [curr_labels; curr_AddLabel];
    
    cd .. 
end

% csvwrite('***_poke_labels_N1.csv', curr_labels)


