%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

Config_NovAna

fps=25; % 15 rgb; % 25; 30 depth; 

startframe=1;%1;%900;
endframe=(fps*60*10)+startframe;%15000;%18000;
Swindow=40;         % Smooth Window Size %40;

% mouseSet = '7day_preexposure_combine';
mouseSet = 'Chess_DLC';
% mouseSet = 'Hiking_DLC';
% mouseSet = 'Chess_DLC';

%run(['/media/alex/DataDrive1/MoSeqData/' mouseSet '/' mouseSet '_MoSeq/Mice_Index.m'])
cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' mouseSet])% '_DLC'])

run('MiceIndex')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation (Richard's code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AnalysisDay=3; % 3 = first novelty day

for miceiter=1:length(Mice)
    cd(Mice(miceiter).name);
    
    %pokes = load(['NoveltyResponse_' Mice(miceiter).name '_N3_poke']);
    
    cd('Analyzed_Data'); %%%
    
    pathname = cd; %%%
    PathRoot=[pathname '/']; %%%
    filelist=dir([PathRoot,'*rgb.mat']); %%%%
    flen = length(filelist); %%%
    load(filelist(AnalysisDay+0).name, 'Labels'); %%%
    disp(['radius_cm: ' num2str(radius_cm)])
    DisThreshold=radius_cm;   % Distance threshold % 10; 60; %radius;
    
    load('Arena_Obj_Pos.mat')
    Labels17 = Labels(:,17);     % head distance from object
    %Labels17 = sqrt((obj_center(1,1)-Labels(:,2)).^2+(obj_center(1,2)-Labels(:,3)).^2)/ppc;
    bodyLen = sqrt( (Labels(:,11)-Labels(:,2)).^2 + ...
                    (Labels(:,12)-Labels(:,3)).^2 )/ppc;
                
                
    % Distance
    Xtime=startframe:endframe;
    Dis=Labels17(startframe:endframe,1); %%%
    SDis=smoothdata(Labels17(startframe:endframe,1),'rloess',Swindow); %%%

    SDis_whole=smoothdata(Labels17(:,1),'rloess',Swindow); %%%
    DisMin=islocalmin(SDis_whole);
    DisMax=islocalmax(SDis_whole);
    
    findDisMax = find(DisMax);

    % Derivative of Distance
    DisVelocity=diff(Labels17(startframe:endframe,1)');
%     DisVelocity=diff(Dis);
    DisVelocity=[0 DisVelocity];
    SDisVelocity=smoothdata(DisVelocity,'rloess',Swindow);
 
    Mice(miceiter).PokingLabels=[];
    Mice(miceiter).BoutApproach=[];
 
    
    if(0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    windowiter=startframe;
    while windowiter<endframe
        % mouse crosses into radius around object
%         if SDis_whole(windowiter)<DisThreshold
        if Labels17(windowiter)<DisThreshold
            
            minfound=0;
            for wpointer=1:length(Labels17)
                % within bout, find local minimum (poke)
                if minfound==0
                    if DisMin(windowiter+wpointer)==1
                        %disp(['minfound: ' num2str(windowiter)])
                        minfound=1;
                        Mice(miceiter).PokingLabels=[Mice(miceiter).PokingLabels windowiter+wpointer];
                        
                        % find beginning of bout (last local maximum before poke)
                        prevBoutApproachAll = find(findDisMax<(windowiter+wpointer),1,'last');
                        prevBoutApproach = findDisMax(prevBoutApproachAll);
                        Mice(miceiter).BoutApproach=[Mice(miceiter).BoutApproach prevBoutApproach];
                    end
                end
                
                % start over when mouse leaves radius
%                 if SDis_whole(windowiter+wpointer)>DisThreshold
                if Labels17(windowiter+wpointer)>DisThreshold
                    break
                end
            end
            windowiter=windowiter+wpointer;
        end
        windowiter=windowiter+1;
    end
    
%     if(0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DisPlot=figure(1);
        plot(Xtime,SDis)
        hold on
        %     plot(Xtime/fps,Labels(startframe:endframe,18),'LineWidth',1.5)
        %     plot(Xtime/fps,Labels(startframe:endframe,19),'LineWidth',1.5)
        %     plot(Xtime/fps,Labels(startframe:endframe,20),'LineWidth',1.5)
        scatter(Mice(miceiter).PokingLabels,SDis_whole(Mice(miceiter).PokingLabels))
        scatter(Mice(miceiter).BoutApproach,SDis_whole(Mice(miceiter).BoutApproach))
        line([startframe endframe], [DisThreshold DisThreshold], 'color', 'k')
        %     line([poke_byHand-900 poke_byHand-900]', repmat([0 350]', 1, length(poke_byHand)), 'color', 'k')
        set(DisPlot, 'Position', [44 296 1871 505])
        title(['Threshold: ', num2str(DisThreshold), 'cm'])
        ylim([0 60])
        
        
        
        
        xPos = Labels(:,2);
        yPos = Labels(:,3);
        correctPokingLabels = Mice(miceiter).PokingLabels(Mice(miceiter).BoutApproach(1)...
            <Mice(miceiter).PokingLabels);
        
        boutPlot = figure(2);
        hold on
        for int = 1:length(Mice(miceiter).BoutApproach)
            plot(xPos(Mice(miceiter).BoutApproach(int):correctPokingLabels(int)), ...
                yPos(Mice(miceiter).BoutApproach(int):correctPokingLabels(int)))%, 'k')
        end
        plot(obj_center(1,1), obj_center(1,2), 'r*')
        
        th = 0:pi/50:2*pi;
        x  = obj_center(1,1);
        y  = obj_center(1,2);
        xunit = DisThreshold*ppc * cos(th) + x; %radius
        yunit = DisThreshold*ppc * sin(th) + y; %radius
        plot(xunit, yunit)
        xlim([0 video_xlen])
        ylim([0 video_ywid])
        
        axis square
        set(gca, 'YDir', 'reverse')
        title(['Threshold: ', num2str(DisThreshold), 'cm'])
        
        
        bodyLenFigure=figure(3);
        scatter(Labels17(startframe:endframe), bodyLen(startframe:endframe), 20,'filled');
        %     scatter(Labels17(startframe:endframe), DisVelocity, 20,'filled');
        
        title(['Len vs Dist ' Mice(miceiter).name, ' Day', num2str(AnalysisDay), ': '...
            num2str(round((endframe-startframe)/fpm)) 'min'],...
            'Interpreter', 'none');
        xlabel('Distance to object (cm)')
        ylabel('Body length')
        set(bodyLenFigure, 'position', [0 0 1200 900]);
        
        pause
        
        saveas(DisPlot, ['FindingPokeDTW_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '.tif'])
        saveas(boutPlot, ['FindingPokeDTW_boutTrace_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '.tif'])
        saveas(bodyLenFigure, ['BodyLen_vs_DistToObj_', ...
            Mice(miceiter).name, '_Day', num2str(AnalysisDay), '.tif'])
        
        csvwrite(['NoveltyResponse_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '_poke'], ...
            Mice(miceiter).PokingLabels)
        csvwrite(['NoveltyResponse_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '_boutApproach'], ...
            Mice(miceiter).BoutApproach)
        csvwrite(['Bouts_10min_', Mice(miceiter).name, '_Day', num2str(AnalysisDay)], boutsReshape)
    end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pokeMat = zeros(length(pokes), fps*2);
    for pokeiter = 1:length(pokes)
        pokeMat(pokeiter,:) = bodyLen(pokes(pokeiter)-fps:pokes(pokeiter)+(fps-1));
    end
    
    close all
    bodyLenPokeFigure=figure(4);
    hold on
    plot(pokeMat', 'Color', [0.8 0.8 0.8])
    plot(mean(pokeMat,1), 'k')
    
    title(['BodyLen aligned to poke (' Mice(miceiter).name, ' Day', num2str(AnalysisDay), ': '...
        num2str(round((endframe-startframe)/fpm)) 'min)'],...
        'Interpreter', 'none');
    xlabel('Time (s)')
    ylabel('Body length (cm)')
    ylim([4 10])
    set(gca, 'FontSize', 16)
    set(bodyLenPokeFigure, 'position', [0 0 800 700]);
        
    cd ../..
    
    saveas(bodyLenPokeFigure, ['BodyLen_alignPoke_', ...
        Mice(miceiter).name, '_Day', num2str(AnalysisDay), '.tif'])
    
    close all
    
%     DisMin_cut = DisMin(startframe:endframe);
%     DisMax_cut = DisMax(startframe:endframe);
%    DisPlot2=figure(2);
%    plot(Xtime,SDis, 'k-')
%    hold on
%    plot(Xtime(DisMin_cut), SDis(DisMin_cut), 'r.')
%    plot(Xtime(DisMax_cut), SDis(DisMax_cut), 'g.')
%     plot(Xtime/fps,Labels(startframe:endframe,18),'LineWidth',1.5)
%
%     DisPatchY=[0 70 70 0];
%     for actiter=1:length(Mice(miceiter).act(:,1))
%         hold on
%         patch([Mice(miceiter).act(actiter,1) Mice(miceiter).act(actiter,1) ...
%                Mice(miceiter).act(actiter,3) Mice(miceiter).act(actiter,3)]/fps,DisPatchY,[0.6 0.4 0.9],...
%         'FaceAlpha',0.3, 'EdgeColor','none')
%     end
%     xlim([375 445])
%     ylabel('Distance (cm)')
%     xlabel('time (s)')

    
%     VelPlot=figure(2);
%     plot(Xtime/fps,SDisVelocity*fps);
%     hold on
%     plot(Xtime/fps,zeros(1,length(Xtime)),'LineWidth',1.5);
%     DisPatchY=[-30 50 50 -30];
% 
%     for actiter=1:length(Mice(miceiter).act(:,1))
%         hold on
%         patch([Mice(miceiter).act(actiter,1) Mice(miceiter).act(actiter,1) ...
%                Mice(miceiter).act(actiter,3) Mice(miceiter).act(actiter,3)]/fps,DisPatchY,[0.6 0.4 0.9],...
%         'FaceAlpha',0.3, 'EdgeColor','none')
%     end
%     xlim([375 445])
%     ylabel('Derivative of Distance cm/s')
%     xlabel('time (s)')
end
disp('end')


%% alternative bout analysis
AnalysisDay=3; % 3 = first novelty day

for miceiter=1:length(Mice)
    cd(Mice(miceiter).name);
    cd('Analyzed_Data');
    
    pathname = cd;
    PathRoot=[pathname '/'];
    filelist=dir([PathRoot,'*rgb.mat']);
    flen = length(filelist);
    load(filelist(AnalysisDay+0).name, 'Labels');
    disp(['radius_cm: ' num2str(radius_cm)])
    DisThreshold=radius_cm;   % Distance threshold
    
    load('Arena_Obj_Pos.mat')
    Labels17 = Labels(:,17);  % head distance from object
    
    
    bodyLen = sqrt( (Labels(:,11)-Labels(:,2)).^2 + ...
                    (Labels(:,12)-Labels(:,3)).^2 )/ppc;
                
    % Distance
    Xtime=startframe:endframe;
    Dis=Labels17(startframe:endframe,1);
    SDis=smoothdata(Labels17(startframe:endframe,1),'rloess',Swindow);

    SDis_whole=smoothdata(Labels17(:,1),'rloess',Swindow);
    DisMin=islocalmin(SDis_whole);
    DisMax=islocalmax(SDis_whole);

    Mice(miceiter).PokingLabels=[];
    Mice(miceiter).BoutApproach=[];
    
    % Bout analysis
    crossTmp  = crossing(Labels(startframe:endframe,21), [], 0.5);
    if(mod(length(crossTmp),2))
        crossTmp = crossTmp(1:end-1);
    end
    cross_In  = crossTmp(1:2:end);
    cross_Out = crossTmp(2:2:end);
    
    for crossiter = 1:length(cross_In)
        pokecurr = find(DisMin(cross_In(crossiter):cross_Out(crossiter)),1);
        
        if(isempty(pokecurr))
            disp('poke not found')
            continue
        end
        
        approachcurr = find(DisMax(1:cross_In(crossiter)),1,'last');
                
        Mice(miceiter).PokingLabels=[Mice(miceiter).PokingLabels pokecurr+cross_In(crossiter)];
        Mice(miceiter).BoutApproach=[Mice(miceiter).BoutApproach approachcurr];
    end
    
    close all
    
    % Scatter plot of distance to object w/ pokes
    DisPlot=figure(1);
    plot(Xtime,Dis)
    hold on
    scatter(Mice(miceiter).PokingLabels,Dis(Mice(miceiter).PokingLabels))
    scatter(Mice(miceiter).BoutApproach,Dis(Mice(miceiter).BoutApproach))
    line([startframe endframe], [DisThreshold DisThreshold], 'color', 'k')
    line([cross_In; cross_In], repmat([-5 0]', 1, length(cross_In)))
    set(DisPlot, 'Position', [44 296 1871 505])
    title(['Threshold: ', num2str(DisThreshold), 'cm'])
    ylim([0 60])
    
    
    xPos = Labels(:,2);
    yPos = Labels(:,3);
    correctPokingLabels = Mice(miceiter).PokingLabels(Mice(miceiter).BoutApproach(1)...
        <Mice(miceiter).PokingLabels);
    
    % Traces of each approach to poke
    boutPlot = figure(2);
    hold on
    for int = 1:length(Mice(miceiter).BoutApproach)
        plot(xPos(Mice(miceiter).BoutApproach(int):correctPokingLabels(int)), ...
            yPos(Mice(miceiter).BoutApproach(int):correctPokingLabels(int)))%, 'k')
    end
    plot(obj_center(1,1), obj_center(1,2), 'r*')
    
    th = 0:pi/50:2*pi;
    x  = obj_center(1,1);
    y  = obj_center(1,2);
    xunit = DisThreshold*ppc * cos(th) + x; %radius
    yunit = DisThreshold*ppc * sin(th) + y; %radius
    plot(xunit, yunit)
    xlim([0 video_xlen])
    ylim([0 video_ywid])
    
    axis square
    set(gca, 'YDir', 'reverse')
    title(['Threshold: ', num2str(DisThreshold), 'cm'])
    
    if(0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        bodyLenFigure=figure(3);
        scatter(Labels17(startframe:endframe), bodyLen(startframe:endframe), 20,'filled');
        %     scatter(Labels17(startframe:endframe), DisVelocity, 20,'filled');
        
        title(['Len vs Dist ' Mice(miceiter).name, ' Day', num2str(AnalysisDay), ': '...
            num2str(round((endframe-startframe)/fpm)) 'min'],...
            'Interpreter', 'none');
        xlabel('Distance to object (cm)')
        ylabel('Body length')
        set(bodyLenFigure, 'position', [0 0 1200 900]);
        
        pause
        
        saveas(DisPlot, ['FindingPokeDTW_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '.tif'])
        saveas(boutPlot, ['FindingPokeDTW_boutTrace_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '.tif'])
        saveas(bodyLenFigure, ['BodyLen_vs_DistToObj_', ...
            Mice(miceiter).name, '_Day', num2str(AnalysisDay), '.tif'])
        
        csvwrite(['NoveltyResponse_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '_poke'], ...
            Mice(miceiter).PokingLabels)
        csvwrite(['NoveltyResponse_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '_boutApproach'], ...
            Mice(miceiter).BoutApproach)
        csvwrite(['Bouts_10min_', Mice(miceiter).name, '_Day', num2str(AnalysisDay)], boutsReshape)

    
    pokeMat = zeros(length(pokes), fps*2);
    for pokeiter = 1:length(pokes)
        pokeMat(pokeiter,:) = bodyLen(pokes(pokeiter)-fps:pokes(pokeiter)+(fps-1));
    end
    
    close all
    bodyLenPokeFigure=figure(4);
    hold on
    plot(pokeMat', 'Color', [0.8 0.8 0.8])
    plot(mean(pokeMat,1), 'k')
    
    title(['BodyLen aligned to poke (' Mice(miceiter).name, ' Day', num2str(AnalysisDay), ': '...
        num2str(round((endframe-startframe)/fpm)) 'min)'],...
        'Interpreter', 'none');
    xlabel('Time (s)')
    ylabel('Body length (cm)')
    ylim([4 10])
    set(gca, 'FontSize', 16)
    set(bodyLenPokeFigure, 'position', [0 0 800 700]);
        
    cd ../..
    
    saveas(bodyLenPokeFigure, ['BodyLen_alignPoke_', ...
        Mice(miceiter).name, '_Day', num2str(AnalysisDay), '.tif'])
    
    close all
    
    end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
disp('end')


%% Saving

saveas(DisPlot, ['FindingPokeDTW_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '.tif'])
saveas(boutPlot, ['FindingPokeDTW_boutTrace_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '.tif'])

cd ..

csvwrite(['NoveltyResponse_', Mice(miceiter).name, '_N', num2str(AnalysisDay), '_poke'], ...
    Mice(miceiter).PokingLabels)
csvwrite(['NoveltyResponse_', Mice(miceiter).name, '_N', num2str(AnalysisDay), '_boutApproach'], ...
    Mice(miceiter).BoutApproach)

cd ..

close all
disp('saved')

%save(['NoveltyResponse_Day' num2str(AnalysisDay)])
%poking_labels = extractfield(Mice,'PokingLabels');
%csvwrite([num2str(mouseSet), '_poke_labels_Day', num2str(AnalysisDay), '_auto'], poking_labels)

%% add back in datanum, index, and act fields
load('NoveltyResponse_Day3.mat')

Mice(1).datanum=41;
Mice(2).datanum=29;
Mice(3).datanum=21;
Mice(4).datanum=45;
Mice(5).datanum=19;
Mice(6).datanum=26;
Mice(7).datanum=20;
Mice(8).datanum=30;
AllActLabels=csvread('CvsS_poke_labels_N1_auto')';
    
Mice(1).index=1;
for miceiter=2:length(Mice)
    Mice(miceiter).index=Mice(miceiter-1).index+Mice(miceiter-1).datanum;
    endindex=Mice(miceiter).index+Mice(miceiter).datanum;
end

if endindex ~= length(AllActLabels)+1
   error('index calculation error');
end


for miceiter=1:length(Mice)
   Mice(miceiter).act=AllActLabels(Mice(miceiter).index:Mice(miceiter).index+Mice(miceiter).datanum-1,:);
end

save('NoveltyResponse_Day3.mat')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization (sanity check)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

% load MoSeq position information
cd /media/alex/DataDrive1/MoSeqData/MSFP_Test/DataSet180922
load('MoSeqDataFrame.mat')

% load depth_ts and rgb_ts for conversion
cd session_20180924185040/
delimiter  = ' ';
formatSpec = '%*q%f%[^\n\r]';

filename = 'depth_ts.txt';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
    'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
depthts = [dataArray{1}];

filename = 'rgb_ts.txt';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
    'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
rgbts = [dataArray{1}];

% load DLC tracking and object position
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/MSFP_Test/MSFP/Analyzed_Data
load('MSFP_Test_180924_rgb.mat')
load('Arena_Obj_Pos.mat')
ppc = 42/6.3;

% extract distance from object for DLC and MoSeq
uuid_subset        = ismember(MoSeqDataFrame.uuid, MoSeqDataFrame.session_uuid(1,:), 'rows');
session_uuid_start = [find(uuid_subset,1,'first'); find(~uuid_subset,1,'first')];

xPos = MoSeqDataFrame.centroid_x_mm(session_uuid_start(2):end);
yPos = MoSeqDataFrame.centroid_y_mm(session_uuid_start(2):end);

xyPos_MoSeq = sqrt((obj_center(2,1)-xPos).^2+(obj_center(2,2)-yPos).^2);
xyPos_DLC   = Labels(:,17);

% zscore to better compare y-values from each array
xPos(isnan(xPos))=0;
yPos(isnan(yPos))=0;
xPos_zscore = zscore(xPos);
yPos_zscore = zscore(yPos);
xLabels_zscore = zscore(Labels(:,14));
yLabels_zscore = zscore(Labels(:,15));

xPos_rescale = rescale(xPos, min(Labels(1:15000,14))-5, max(Labels(1:15000,14))+12);
yPos_rescale = rescale(yPos, min(Labels(1:15000,15))-5, max(Labels(1:15000,15)));

xyPos_MoSeq(isnan(xyPos_MoSeq))=0;
xyPos_MoSeq_zscore = zscore(xyPos_MoSeq);
xyPos_DLC_zscore = zscore(-xyPos_DLC); % -1


% translate automatically identified pokes from 25fps (DLC) to 30fps (MoSeq)
poke_byHand = csvread('NoveltyResponse_byHand_korleki.csv'); % depth_ts units
poke_auto   = load('NoveltyResponse_MSFP_N1')';              % DLC output units
rescaling   = length(xyPos_DLC)/length(depthts);             % (DLC output units / depth actual units)
poke_auto   = round(poke_auto/rescaling);                    % depth_ts units
%rescaling   = length(xyPos_MoSeq)/length(xyPos_DLC);
%x_new = (1:length(xyPos_DLC))*rescaling;

poke_adj = zeros(length(poke_auto),1);
for i = 1:length(poke_auto)
    [d, ix] = min(abs(rgbts(poke_auto(i))-depthts));
    poke_adj(i,1) = ix;
end

% define circle around obj center
th = 0:pi/50:2*pi;
xunit = radius * cos(th) + obj_center(2,1);
yunit = radius * sin(th) + obj_center(2,2);

%%
close all

figure(1)
hold on
%plot(xPos_zscore(1:18000),yPos_zscore(1:18000), 'k-')
%plot(xLabels_zscore(1:15000), yLabels_zscore(1:15000), 'r-')
%plot(xPos_zscore(poke_byHand), yPos_zscore(poke_byHand), 'g*')
%plot(xPos_zscore(poke_adj), yPos_zscore(poke_adj), 'c*')

plot(Labels(1:15000,14), Labels(1:15000,15), 'r-')
plot(xPos_rescale(1:18000), yPos_rescale(1:18000), 'k-')
plot(Labels(poke_auto,14), Labels(poke_auto,15), 'g*')

plot(xunit, yunit, 'r')
legend({'DLC', 'MoSeq', 'poke auto'})
set(gca,'ydir','reverse')


axis square


fig2 = figure(2);
set(fig2, 'Position', [46 26 1865 450])
hold on
%plot(xyPos_MoSeq, 'k-')
%plot(xyPos_MoSeq_zscore, 'k-')
%plot(x_new, xyPos_DLC_zscore, 'r.-')
%plot(poke_byHand, xyPos_MoSeq_zscore(poke_byHand), 'g*')
%plot(poke_adj, xyPos_DLC_zscore(poke_adj), 'c*')
% line([poke_byHand poke_byHand]', repmat([5 45], length(poke_byHand), 1)', ... %[-2.5 2]
%     'color', 'g')

plot(1:length(xyPos_DLC), xyPos_DLC, 'r-')
plot(poke_auto, xyPos_DLC(poke_auto), 'g*')
%line([poke_auto poke_auto]', repmat([5 45], length(poke_auto), 1)', ... %[-2.5 2]
%    'color', 'g')


xlim([0 20000])
legend({'MoSeq', 'DLC', 'poke_byHand', 'poke_auto'})


%% convert DLC output units to rgb_ts units
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry_DLC/

for miceiter = 1:length(Mice)
    cd(['./' Mice(miceiter).name '/Analyzed_Data'])
    filename = dir('NoveltyResponse_*_Day3');
    poke_auto = load(filename.name);
    
    cd ..
    rgbts_list = dir('*rgb_ts.txt');
    rgbts = load(rgbts_list(3).name);
    rgbtsLen = length(rgbts);
    
    csv_list = dir('*csv');
    DLC_output = csvread(csv_list(3).name, 3, 0);
    DLC_out_len = length(DLC_output);
    
    rescaling = rgbtsLen/DLC_out_len;
    poke_adj = round(poke_auto*rescaling);
    
    cd Analyzed_Data
    csvwrite(['NoveltyResponse_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), ...
        '_rgbAdj'], poke_adj)
    
    cd ../..
end

%% consolidate bout information (if each day is saved separately)
days = 5;
boutsNum       = zeros(length(Mice),days);
boutsLenMean   = zeros(length(Mice),days);
boutsLenMedian = zeros(length(Mice),days);    

for miceiter=1:length(Mice)
    cd(Mice(miceiter).name);
    cd('Analyzed_Data'); %%%
    %cd head
    
    boutFiles = dir('Bouts*');
    for dayiter = 1:length(boutFiles)
        boutFileCurr = load(boutFiles(dayiter).name);
        
        boutsNumCurr = size(boutFileCurr,2);
        boutsLenCurr = boutFileCurr(2,:)-boutFileCurr(1,:);
        boutsNum(miceiter,dayiter) = boutsNumCurr;
        boutsLenMean(miceiter,dayiter) = mean(boutsLenCurr)/fps;
        boutsLenMedian(miceiter,dayiter) = median(boutsLenCurr)/fps;
    end
    
    cd ../../
end

csvwrite([mouseSet, '_boutNum'],boutsNum)
csvwrite([mouseSet, '_avgBoutLen'],boutsLenMean)
csvwrite([mouseSet, '_medBoutLen'],boutsLenMedian)
