% in subfolder when looking for .mat files, prsumed there is an Arena_Obj_Pos.mat file
% omit the first .mat file which is Arena_Obj_Pos.mat. 
% If there are other .mat files in subfolder will cause error

clear
close all
clc

Config_NovAna_Ghana

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/StandardLesion_combine/

% start_min=0.5;
% end_min=start_min+10;
% startframe=uint16(start_min.*60.*fps);
% endframe=uint16(end_min.*60.*fps);

durTotal = 10; % duration of analysis (min)
disp(['Duration of analysis: ' num2str(durTotal) 'min'])

fps=15; %%%%%%%%%%%%%%%%%%%%%%%%%%%

Dis_ts_frame=500;
Dis_te_frame=durTotal.*60.*fps+Dis_ts_frame;
startframe  = Dis_ts_frame;
endframe    = Dis_te_frame;

folderpath = cd;
folderd = dir(folderpath);
isub = [folderd(:).isdir];
foldernames = {folderd(isub).name}'; 
foldernames(ismember(foldernames,{'.','..','Retrain_Sep17','temp','Trajectories_H2N1'})) = []; 
folderlen=length(foldernames);

% smoothing window for SAP+orientation analysis
Swindow=30;
Swindow_orient=40;

% maximum body length (pixels), anything above counted as error
bodyLen_cutoff = 80;

disp('next')

%%
% shift_time = input('Exclude first 5 min for some mice? 0/1: ');
shift_time = 0;
if(shift_time==1)
    whichFiles_shift = [4 5 6];
    disp(['Mice where lego is placed in arena at 5min: ' num2str(whichFiles_shift)])
else
    shift_time = 0;
end

for folderi=1:folderlen
    cd(foldernames{folderi});
    disp(foldernames{folderi})
%     poke_file = dir('*poke');
%     disp(poke_file.name)
%     pokes = csvread(poke_file.name);
    dist_of_sap(folderi).name = foldernames{folderi};
    
    cd Analyzed_Data_1obj_12cm_tail %%%%%%%%%%%%%%%%%%%%%%%%%
    load('Arena_Obj_Pos.mat','arena')
%     cd ./nose_8cm

    subpath=cd;
    PathRoot=[subpath '/'];
%     filelist=dir('*rgb.mat'); 
%     filelist=dir('session01*.mat');
    filelist=dir('*rgb_Converted*.mat'); %PathRoot, '*.mat']); foldernames{folderi}(4:end) 
    flen = length(filelist);
    
    for filei = 1:flen
        disp(filei)
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
        orientSmooth = smoothdata(LabelsCut(:,22),'rloess',Swindow_orient);
        isOrient = intersect(find(orientSmooth<=angle_radius), find(orientSmooth>=(-angle_radius)));
        Time_angle(filei) = length(isOrient);
        
%         % sanity orientation histogram
%         histCurr = figure(1);
%         histogram(abs(orientSmooth), 0:15:180)
%         title([foldernames{folderi} ': Day' num2str(filei)])
%         xlabel('abs(angle)')
%         ylim([0 2500])
%         saveas(histCurr,['OrientHist_' [foldernames{folderi} '_Day' num2str(filei)] '.tif'])
        
        % Analyse body length and SAP
        smooth_headX = smoothdata( LabelsCut(:,2), 'rloess', Swindow);
        smooth_headY = smoothdata( LabelsCut(:,3), 'rloess', Swindow);
        smooth_tailX = smoothdata( LabelsCut(:,11), 'rloess', Swindow);
        smooth_tailY = smoothdata( LabelsCut(:,12), 'rloess', Swindow);
        curr_velX  = LabelsCut(:,24);
        curr_velY  = LabelsCut(:,25);

        bodyLen = sqrt( (smooth_tailX-smooth_headX).^2 + (smooth_tailY-smooth_headY).^2 ); % pixels
        bodyLen(bodyLen>bodyLen_cutoff)=0; % body lengths >80px are likely to be errors

        vel        = sqrt( curr_velX.^2 + curr_velY.^2 ); % head velocity, cm/s
        stretches  = find(bodyLen>(mean(bodyLen)+std(bodyLen)));
        slows      = find(vel<5);
        sap        = intersect(stretches,slows);
        sap_in_rad = intersect(sap, find(LabelsCut(:,21)==1));
        
        diffSap      = diff(sap);
        diffSapInRad = diff(sap_in_rad);
        sap_frame{filei,1}    = sap(find(diffSap>1)+1);
        sap_num(filei)        = sum(diffSap>1); % total SAPs (any location)
        sap_num_in_rad(filei) = sum(diffSapInRad>1); % SAPs near object
%         sap_num_in_rad(filei) = sum(diffSapInRad>1)/Time_distance(filei); % (norm)
        sap_dist{filei,1}     = LabelsCut(diffSap>1,17);
        nose_dist{filei,1}    = smooth_headX(diffSap>1);
        nose_dist{filei,2}    = smooth_headY(diffSap>1);
%         sap_frame{filei,1}    = sap;
%         sap_num(filei)        = length(sap); % total SAPs (any location)
%         sap_num_in_rad(filei) = length(sap_in_rad)/Time_distance(filei); % (norm)
%         sap_dist{filei,1}     = LabelsCut(sap,17);
%         nose_dist{filei,1}    = smooth_headX(sap);
%         nose_dist{filei,2}    = smooth_headY(sap);
        
%         load('Arena_Obj_Pos.mat')
%         disp(['radius: ' num2str(radius_cm)])
%         
%         % scatter plot of SAP as a function of distance
%         SAP_figure = figure(1);
% %         plot(LabelsCut(:,2), LabelsCut(:,3),'k.')
%         hold on
%         plot(LabelsCut(sap,2), LabelsCut(sap,3),'b.')
%         
%         th = 0:pi/50:2*pi;
%         x  = obj_center(filei,1);
%         y  = obj_center(filei,2);
%         xunit = radius_cm * ppc * cos(th) + x;
%         yunit = radius_cm * ppc * sin(th) + y;
% %         plot(xunit, yunit,'r-','linewidth',3)
%         xlim([50 500])
%         set(gca, 'YDir', 'reverse')
%         
%         pause(0.5)
                
        % scatter plot of body length as a function of distance to object
%         bodyLen = sqrt( (LabelsCut(:,11)-LabelsCut(:,2)).^2 + ...
%                         (LabelsCut(:,12)-LabelsCut(:,3)).^2 )/ppc;
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
        
        
        % Analyze total distance covered (based on nose/tail avg)
        centerCalcX = (LabelsCut(:,2) + LabelsCut(:,11))/2;
        centerCalcY = (LabelsCut(:,3) + LabelsCut(:,12))/2;
        centerCalc = [centerCalcX centerCalcY]/ppc; % convert pixel to cm
        
        SDis=smoothdata(centerCalc,'rloess',Swindow);
        
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
    
    Time_distance_all(folderi,:)  = Time_distance;
    Time_angle_all(folderi,:)     = Time_angle;
    Time_periphery_all(folderi,:) = Time_periphery;
    TotalDistCut_all(folderi,:)   = totalDistCut;
    SAP_num_all(folderi,:)        = sap_num;
    SAP_in_rad_all(folderi,:)     = sap_num_in_rad;
    dist_of_sap(folderi).sapFrame = sap_frame;
    dist_of_sap(folderi).sapDist  = sap_dist;
    dist_of_sap(folderi).nosePos  = nose_dist;
%     FracArea_all(folderi,:)=fracArea;
    clearvars Time_angle Time_distance Time_periphery totalDistCut ...
        fracArea sap_num sap_num_in_rad sap_frame sap_dist nose_dist
    cd ..
%     cd .. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd ..
    DisOrAng{folderi,1}='Distance';
    DisOrAng{folderlen+folderi,1}='Angle';
    Periph{folderi,1}='Periphery';
    TotalDist{folderi,1}='DistanceRun';
    SAPNumNear{folderi,1}='SAPnum';
    SAPNumNear{folderlen+folderi,1}='SAPnear';
%     TotalArea{folderi,1}='AreaCovered';

end

Time_distance_all = Time_distance_all./double(endframe-startframe);
Time_angle_all    = Time_angle_all./double(endframe-startframe);

Table  = table(cat(1,foldernames,foldernames),DisOrAng,[Time_distance_all;Time_angle_all]);
Table2 = table(cat(1,foldernames), Periph, Time_periphery_all);
Table3 = table(cat(1,foldernames), TotalDist, TotalDistCut_all);
Table4 = table(cat(1,foldernames,foldernames),SAPNumNear,[SAP_num_all;SAP_in_rad_all]);

if(0)
    writetable(Table,'TimeStatistic.csv');
    writetable(Table2,'TimeStatistic_body_periph.csv');
    writetable(Table3,'TimeStatistic_nose_totalDistCut.csv');
    writetable(Table4,'TimeStatistic_SAP_norm.csv');
    save('SAP_plus_dist.mat', 'dist_of_sap')
    
    disp('saved')
end

%clear all
disp('end')