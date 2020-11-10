%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: use bout_analysis.m script instead ~CKA 190813

clear
close all
clc

Config_NovAna

fps=15; % 15 rgb; % 25 unconverted; 30 depth; 

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
