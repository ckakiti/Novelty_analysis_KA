% in subfolder when looking for .mat files, prsumed there is an Arena_Obj_Pos.mat file
% omit the first .mat file which is Arena_Obj_Pos.mat. 
% If there are other .mat files in subfolder will cause error

clear
close all
clc

Config_NovAna

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Planets_DLC/
%cd('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Holidays/')

% start_min=0.5;
% end_min=start_min+10;
% startframe=uint16(start_min.*60.*fps);
% endframe=uint16(end_min.*60.*fps);

durTotal = 30; % duration of analysis (min)
disp(['Duration of analysis: ' num2str(durTotal) 'min'])

Dis_ts_frame=500;
Dis_te_frame=durTotal.*60.*fps+Dis_ts_frame;
startframe  = Dis_ts_frame;
endframe    = Dis_te_frame;

folderpath = cd;
folderd = dir(folderpath);
isub = [folderd(:).isdir];
foldernames = {folderd(isub).name}'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foldernames(ismember(foldernames,{'.','..','Retrain_Sep17','temp'})) = []; %%%%%%%%%%%%%%%%%%%%%%%
folderlen=length(foldernames);

% for orientation analysis (smoothing)
Swindow=40;

disp('next')

%%
shift_time = input('Exclude first 5 min for some mice? 0/1: ');
if(shift_time==1)
    whichFiles_shift = [4 5 6];
    disp(['Mice where lego is placed in arena at 5min: ' num2str(whichFiles_shift)])
end

for folderi=1:folderlen
    cd(foldernames{folderi});
    disp(foldernames{folderi})
%     poke_file = dir('*poke');
%     disp(poke_file.name)
%     pokes = csvread(poke_file.name);
    
    cd Analyzed_Data_1obj_12cm_tail
    load('Arena_Obj_Pos.mat','arena')
%     cd ./tail %%%%%%%%%%%%%%%%

    subpath=cd;
    PathRoot=[subpath '/'];
    filelist=dir('*rgb_Converted.mat'); %*rgb_Converted.mat'); %'session01*.mat']; PathRoot, '*.mat']); foldernames{folderi}(4:end) %%%%%%%%%%%%%%%%%%
    flen = length(filelist);
    
    for filei = 1:flen
        filename = filelist(filei).name;
        load(filename)
        
        % exclude first 5 min of recording (if lego is placed in arena after 5 min)
        if(shift_time==1 && any(folderi==whichFiles_shift))
            timeShift = 5000;
            disp(folderi)
            disp(timeShift)
            
            startframe = Dis_ts_frame+timeShift;
            endframe   = min(Dis_te_frame+timeShift, size(Labels,1));
        else
            startframe = Dis_ts_frame;
            endframe   = min(Dis_te_frame, size(Labels,1));
        end

        LabelsCut=Labels(startframe:endframe,:);
        
        % Analyze time spent near obj (within specified radius)
        Time_distance(filei) = sum(LabelsCut(:,21));
        
        
        % Analyze time spent oriented to obj (after smoothing)
        orientSmooth = smoothdata(LabelsCut(:,22),'rloess',Swindow);
        isOrient = intersect(find(orientSmooth<=angle_radius), find(orientSmooth>=(-angle_radius)));
        Time_angle(filei) = length(isOrient);
        
%         % sanity orientation histogram
%         histCurr = figure(1);
%         histogram(abs(orientSmooth), 0:15:180)
%         title([foldernames{folderi} ': Day' num2str(filei)])
%         xlabel('abs(angle)')
%         ylim([0 2500])
%         saveas(histCurr,['OrientHist_' [foldernames{folderi} '_Day' num2str(filei)] '.tif'])
        

%         % scatter plot of body length as a function of distance to object
        bodyLen = sqrt( (LabelsCut(:,11)-LabelsCut(:,2)).^2 + ...
                        (LabelsCut(:,12)-LabelsCut(:,3)).^2 )/ppc;
% 
%         bodyLenFigure=figure(3);
%         scatter(LabelsCut(:,17), bodyLen, 20,'filled');
%         title(['Len vs Dist ' foldernames{folderi}, ': Day', num2str(filei)],...
%             'Interpreter', 'none');
%         ylim([0 15])
%         xlabel('Distance to object (cm)')
%         ylabel('Body length')
%         set(bodyLenFigure, 'position', [0 0 1200 900]);

        
%         % plot body length aligned to each poke
%         pokeMat = zeros(length(pokes), fps*2);
%         for pokeiter = 1:length(pokes)
%             pokeMat(pokeiter,:) = bodyLen(pokes(pokeiter)-fps:pokes(pokeiter)+(fps-1));
%         end
        
%         close all        
%         bodyLenPokeFigure=figure(4);
%         hold on
%         plot(pokeMat', 'Color', [0.8 0.8 0.8])
%         plot(mean(pokeMat,1), 'k')
%         
%         title(['Body length aligned to pokes (' foldernames{folderi}, ': Day', num2str(filei) ')'],...
%             'Interpreter', 'none');
%         xlabel('Time (s)')
%         ylabel('Body length (cm)')
%         ylim([0 15])
%         curr_xtick = get(gca, 'xtick');
%         new_xtick = arrayfun(@num2str, curr_xtick-fps, 'Uniform', false);
%         set(gca,'xticklabel', new_xtick)
%         set(gca, 'FontSize', 16)
%         set(bodyLenPokeFigure, 'position', [0 0 800 700]);
%         saveas(bodyLenPokeFigure,['BodyLenPoke_' [foldernames{folderi} '_Day' num2str(filei)] '.tif'])
        
        
        % Analyze time spent in periphery
        arenaCorners = [arena(filei,1) arena(filei,3) arena(filei,3) arena(filei,1) arena(filei,1); ...
                        arena(filei,2) arena(filei,2) arena(filei,4) arena(filei,4) arena(filei,2)];
        arenaCenter  = [(arena(filei,1)+arena(filei,3))/2 (arena(filei,2)+arena(filei,4))/2];
        
        scale = 190;
        L  = linspace(0,2*pi,5) + pi/4; % select angles and thus number of sides for polygon,
                                        %   rotate so edges are parallel to arena edges
        xv = scale*cos(L)' + arenaCenter(1); % x points
        yv = scale*sin(L)' + arenaCenter(2); % y points
        
        xq = LabelsCut(:,14);
        yq = LabelsCut(:,15);
        
        [in,on] = inpolygon(xq,yq,xv,yv); % detect points that are inside defined polygon region
        
        % determine fraction of points outside polygon region
        Time_periphery(filei) = numel(xq(~in))/numel(xq);
        
        
        % Analyze total distance covered (based on nose pos - more reliable than center of mass)
        centerCalcX = (LabelsCut(:,2));% + LabelsCut(:,5) + ...
%                        LabelsCut(:,8) + LabelsCut(:,11))/4;
        centerCalcY = (LabelsCut(:,3));% + LabelsCut(:,6) + ...
%                        LabelsCut(:,9) + LabelsCut(:,12))/4;
        centerCalc = [centerCalcX centerCalcY]/ppc; % convert pixel to cm
        
        Swindow=40;
        SDis=smoothdata(centerCalc,'rloess',Swindow);
%         SDis_pix=smoothdata([centerCalcX centerCalcY],'rloess',Swindow);
        
        cumDist = [0; cumsum( sqrt( diff(SDis(:,1)).^2+diff(SDis(:,2)).^2 ) )];
        totalDistCut(filei) = cumDist(end); % cm
                   
        % Analyze total area covered
%         squareSize   = 4;    % must be a factor of both video_ywid and video_xlen
%         roundTargets = squareSize:squareSize:max(video_xlen,video_ywid);
%         roundedArea  = interp1(roundTargets,roundTargets,[SDis_pix(:,1) SDis_pix(:,2)],'nearest');
%         linearInd    = sub2ind([video_xlen video_ywid], roundedArea(:,1), roundedArea(:,2));
%         runArea      = length(unique(linearInd));
%         
%         allCoords = zeros(round(arena(1,3)), round(arena(1,4)));
%         allCoords(round(arena(1,1)):end, round(arena(1,2)):end) = 1;
%         [allCoordsX, allCoordsY] = find(allCoords==1);
%         roundedArena = interp1(roundTargets,roundTargets,[allCoordsX, allCoordsY],'nearest');
%         linearIndArena = sub2ind([video_xlen video_ywid], roundedArena(:,1), roundedArena(:,2));
%         arenaArea = length(unique(linearIndArena));
%                 
%         fracArea(filei)  = runArea/arenaArea;
        
%         clearvars -except filelist flen filei Time_angle Time_distance Time_periphery ...
%             startframe endframe folderi folderlen foldernames write_pointer arena ...
%             TableName Time_distance_all Time_angle_all Time_periphery_all TotalDistCut_all ...
%             DisOrAng Periph TotalDist totalDistCut fracArea ppc video_xlen video_ywid ...
%             Swindow angle_radius
    end
    
    Time_distance_all(folderi,:)=Time_distance;
    Time_angle_all(folderi,:)=Time_angle;
    Time_periphery_all(folderi,:)=Time_periphery;
    TotalDistCut_all(folderi,:)=totalDistCut;
%     FracArea_all(folderi,:)=fracArea;
    clearvars Time_angle Time_distance Time_periphery totalDistCut fracArea
    cd ..
%     cd .. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd ..
    DisOrAng{folderi,1}='Distance';
    DisOrAng{folderlen+folderi,1}='Angle';
    Periph{folderi,1}='Periphery';
    TotalDist{folderi,1}='DistanceRun';
%     TotalArea{folderi,1}='AreaCovered';

end

Time_distance_all = Time_distance_all./double(endframe-startframe);
Time_angle_all    = Time_angle_all./double(endframe-startframe);

Table  = table(cat(1,foldernames,foldernames),DisOrAng,[Time_distance_all;Time_angle_all]);
Table2 = table(cat(1,foldernames), Periph, Time_periphery_all);
Table3 = table(cat(1,foldernames), TotalDist, TotalDistCut_all);

if(0)
    writetable(Table,'TimeStatistic.csv');
    writetable(Table2,'TimeStatistic_body_periph.csv');
    writetable(Table3,'TimeStatistic_nose_totalDistCut.csv');
    disp('saved')
end

%clear all
disp('end')