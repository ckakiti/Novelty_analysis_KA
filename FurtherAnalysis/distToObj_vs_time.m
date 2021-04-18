%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate temporal development of mouse behaviors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

% config_file = 'Config_NovAna_NewHope_ROTJ.m';
% config_file = 'Config_NovAna_CvsS.m';
config_file = 'Config_NovAna_Ghana.m';
% config_file = 'Config_NovAna_combine3.m';
run(config_file)

fps=15; % 15 rgb; % 25 unconverted; 30 depth; 

startframe=1;%1;%900;
endframe=(fps*60*30)+startframe;%15000;%18000;
Swindow=40;         % Smooth Window Size %40;

%mouseSet = 'NewHope-ROTJ';
%mouseSet = 'CvsS_180831_DLC';
%mouseSet = 'Capoeira_DLC';
%mouseSet = 'Hiking_DLC';
%mouseSet = 'Chess_DLC';
%mouseSet = 'Planets_DLC';
%mouseSet = 'Iku_6OHDA_DLC';
mouseSet = 'Montana_DLC';
%mouseSet = 'StandardSetup_combine';

% groups   = {'Capoeira_DLC', 'Hiking_DLC', 'Chess_DLC'};

setfolder = strcat('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/',mouseSet);
cd(setfolder)
%run('MiceIndex_NewHope_ROTJ')
run('MiceIndex_Montana')
%run('MiceIndex_combine3')

plot_kind = input('Plot nose, tail, or both? n/t/b: ', 's');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%group_n = 1;
%groupfolder = strcat('/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/',groups{group_n});

AnalysisDay=2; % 3 = first novelty day

for miceiter=1:length(Mice)
    clc
%     cd(groupfolder)
    cd(setfolder)
    cd(Mice(miceiter).name)
    
    % load head information
    cd('Analyzed_Data_1obj_8cm_nose')
%     cd('Analyzed_Data_1obj_head')
    pathname = cd;
    PathRoot=[pathname '/'];
    filelist=dir([PathRoot,'*rgb_Converted*.mat']); 
%     filelist=dir([PathRoot,'session01*.mat']); 
    flen = length(filelist);
    load(filelist(AnalysisDay).name, 'Labels');
    Labels_head = Labels;
    clear Labels
    cd ..
    
    % load tail information
    cd('Analyzed_Data_1obj_12cm_tail')
%     cd('Analyzed_Data_1obj_tail')
    pathname = cd;
    PathRoot=[pathname '/'];
    filelist=dir([PathRoot,'*rgb_Converted*.mat']); 
%     filelist=dir([PathRoot,'session01*.mat']); 
    flen = length(filelist);
    load(filelist(AnalysisDay).name, 'Labels');
    Labels_tail = Labels;
    clear Labels
    cd ..
    
    %disp(['radius_cm: ' num2str(radius_cm)])
    %DisThreshold=radius_cm;   % Distance threshold % 10; 60; %radius;
    %load('Arena_Obj_Pos.mat')
    
    headDist = Labels_head(:,17); % head distance from obj
    tailDist = Labels_tail(:,17); % tail distance from obj    
    
    % smooth trajectories for plotting
    if(length(Labels_head)<endframe)
        Xtime=startframe:length(Labels_head);
    else
        Xtime=startframe:endframe;
    end
    %Dis=headDist(Xtime,1);
    SDis_head=smoothdata(headDist(Xtime,1),'rloess',Swindow);
    SDis_tail=smoothdata(tailDist(Xtime,1),'rloess',Swindow);
    bodyLen = abs(SDis_head-SDis_tail);
    bodyLen(bodyLen>10)=10;
    
    close all
    
    if(strcmp(plot_kind,'n'))
        DisPlot=figure(1);
        hold on
        plot(Xtime/fps,log10(SDis_head),'r')
%         scatter(Xtime/fps, (SDis_head), 10, bodyLen)
        line([startframe/fps endframe/fps], log10([8 8]), 'color', [0.9 0.8 0.8])
        
        set(DisPlot, 'Position', [44 400 1871 505])
        set(gca, 'YDir', 'reverse', 'FontSize', 18)
        title(['Distance from object: ' ...
            Mice(miceiter).name ' Day ' num2str(AnalysisDay)], ...
            'Interpreter', 'none')
        legend('nose')
        xlabel('Time (s)')
        ylabel('Distance (log(cm))')
%         ylabel('Distance (cm)')
        xlim([0 length(Xtime)/fps])
        ylim(log10([0.1 60]))
%         ylim([0 60])
%         colorbar
        
        %if(0)
        pause
        cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' mouseSet '/temp/distToObj'])
%         saveas(DisPlot, ['distToObj_vs_time_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '_nose_log.tif'])
        saveas(gcf, ['distToObj_vs_time_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '_nose_log_color.tif'])
        %end
    elseif(strcmp(plot_kind,'t'))
        DisPlot=figure(1);
        hold on
        plot(Xtime/fps,log10(SDis_tail),'k')
        line([startframe/fps endframe/fps], log10([12 12]), 'color', [0.7 0.7 0.7])
        
        set(DisPlot, 'Position', [44 350 1871 505])       
        set(gca, 'YDir', 'reverse', 'FontSize', 18)
        title(['Distance from object: ' ...
            Mice(miceiter).name ' Day ' num2str(AnalysisDay)], ...
            'Interpreter', 'none')
        legend('tail')
        xlabel('Time (s)')
        ylabel('Distance (log(cm))')
        xlim([0 length(Xtime)/fps])
        ylim(log10([0.1 60]))
        
        %if(0)
        pause
        cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' mouseSet '/temp/distToObj'])
        saveas(DisPlot, ['distToObj_vs_time_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '_tail_log.tif'])
        %end
    elseif(strcmp(plot_kind,'b'))
        DisPlot=figure(1);
        hold on
        %plot(Xtime/fps,SDis_head,'r')
        %plot(Xtime/fps,SDis_tail,'k')
        %line([startframe/fps endframe/fps], [8 8], 'color', [0.9 0.8 0.8])
        %line([startframe/fps endframe/fps], [12 12], 'color', [0.7 0.7 0.7])
        
        plot(Xtime/fps/60,log10(SDis_head),'r','LineWidth',1.5)
        plot(Xtime/fps/60,log10(SDis_tail),'k','LineWidth',1.5)
%         line([startframe/fps endframe/fps], log10([8 8]), 'color', [0.9 0.8 0.8])
%         line([startframe/fps endframe/fps], log10([12 12]), 'color', [0.7 0.7 0.7])
        
        set(DisPlot, 'Position', [44 450 1871 505])
        set(gca, 'YDir', 'reverse', 'FontSize', 18, 'tickdir','out')
        title(['Distance from object: ' ...
            Mice(miceiter).name ' Day ' num2str(AnalysisDay)], ...
            'Interpreter', 'none')
        legend({'nose', 'tail'})
        %legend('nose')
        %legend('tail')
        xlabel('Time (min)')
        ylabel('Distance (cm)')
%         ylabel('Distance (log(cm))')
        yticks([-1 0 1 log10(60)])
        yticklabels({'0.1','1','10','60'})
        xlim([0 length(Xtime)/fps/60])
        ylim(log10([0.1 70]))
        %ylim([0 60])
        
%         if(0)
        pause(0.1)
        cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' mouseSet '/temp/distToObj'])
        saveas(DisPlot, ['distToObj_vs_time_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '_log_edit.tif'])
        %saveas(DisPlot, ['distToObj_vs_time_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '_nose_log.tif'])
        %saveas(DisPlot, ['distToObj_vs_time_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '_tail_log.tif'])
%         end
    end
    
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
    headDist = Labels(:,17);  % head distance from object
    
    
    bodyLen = sqrt( (Labels(:,11)-Labels(:,2)).^2 + ...
                    (Labels(:,12)-Labels(:,3)).^2 )/ppc;
                
    % Distance
    Xtime=startframe:endframe;
    Dis=headDist(startframe:endframe,1);
    SDis_head=smoothdata(headDist(startframe:endframe,1),'rloess',Swindow);

    SDis_whole=smoothdata(headDist(:,1),'rloess',Swindow);
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
