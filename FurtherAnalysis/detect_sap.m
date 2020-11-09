%%
clear
close all
clc

Config_NovAna
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/

curr_set = 'Capoeira_DLC';
curr_day = 3;  % N1
durTotal = 10; % duration of analysis (min)
disp(['Duration of analysis: ' num2str(durTotal) 'min'])

Dis_ts_frame = 500;
Dis_te_frame = durTotal.*60.*fps+Dis_ts_frame;
startframe   = Dis_ts_frame;
endframe     = Dis_te_frame;

cd(curr_set)
%poke_file = dir('*poke_labels_N1.csv');
%pokes     = csvread(poke_file.name);

folderpath = cd;
folderd = dir(folderpath);
isub = [folderd(:).isdir];
foldernames = {folderd(isub).name}'; 
foldernames(ismember(foldernames,{'.','..','Retrain_Sep17','temp'})) = []; 
folderlen=length(foldernames);

disp('next')

%%
cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' curr_set])

% for smoothing
Swindow=30;

% maximum body length (pixels), anything above counted as error
bodyLen_cutoff = 80;

for miceiter = 1:folderlen
    close all
    clc
    
    cd(foldernames{miceiter})
    cd Analyzed_Data_1obj_tail
    %cd Analyzed_Data_10cm
    
    mat_files = dir('*Converted.mat'); %'Converted.mat') %'rgb.mat')
    load(mat_files(curr_day).name)
    LabelsCut = Labels(startframe:endframe,:);
    
    smooth_headX = smoothdata( LabelsCut(:,2), 'rloess', Swindow);
    smooth_headY = smoothdata( LabelsCut(:,3), 'rloess', Swindow);
    smooth_tailX = smoothdata( LabelsCut(:,11), 'rloess', Swindow);
    smooth_tailY = smoothdata( LabelsCut(:,12), 'rloess', Swindow);
    curr_velX  = LabelsCut(:,24);
    curr_velY  = LabelsCut(:,25);
    
    bodyLen    = sqrt( (smooth_tailX-smooth_headX).^2 + (smooth_tailY-smooth_headY).^2 ); % pixels
    bodyLen(bodyLen>bodyLen_cutoff)=0; % body lengths >80px are likely to be errors
    
    vel        = sqrt( curr_velX.^2 + curr_velY.^2 ); % head velocity, cm/s
    stretches  = find(bodyLen>(mean(bodyLen)+std(bodyLen)));
    slows      = find(vel<5);
    sap        = intersect(stretches,slows);
    sap_in_rad = intersect(sap, find(LabelsCut(:,21)==1));
    
    diffSap      = diff(sap);
    diffSapInRad = diff(sap_in_rad);
    sap_frame    = sap(find(diffSap>1)+1);
    sap_num        = sum(diffSap>1); % total SAPs (any location)
    sap_num_in_rad = sum(diffSapInRad>1); % SAPs near object
    sap_dist      = LabelsCut(diffSap>1,17);
    nose_distX    = smooth_headX(diffSap>1);
    nose_distY    = smooth_headY(diffSap>1);
        
%     curr_pokes = find(pokes(:,1)==miceiter);
    
    % plot trajectories for nose and tail side-by-side
%     fig1 = figure(1);
%     set(gcf, 'Position', [100 300 1100 450])
%     subplot(1,2,1)
%     plot(curr_headX, curr_headY, 'k.-')
%     axis square
%     subplot(1,2,2)
%     plot(curr_tailX, curr_tailY, 'k.-')
%     axis square
    
    % plot body length (distance between nose and tail), superimpose automatically calculated poke 
%     fig2 = figure(2);
%     hold on
%     plot(bodyLen, 'k.-')
%     line([pokes(curr_pokes,2) pokes(curr_pokes,2)]', repmat([0 max(bodyLen(:))], length(curr_pokes), 1)')
    
    % histogram of body lengths (sanity check)
    fig3 = figure(3);
    set(gcf, 'Position', [300 500 560 420])
    hold on
    hist(bodyLen,100)
    sap_threshold = mean(bodyLen)+std(bodyLen);
    line([sap_threshold sap_threshold]', (fig3.CurrentAxes.YLim)', 'Color', 'r')
    title(['Body lengths + sap threshold: ' foldernames{miceiter} '_N1_' num2str(durTotal) 'min'], 'Interpreter', 'none')
    set(gca, 'FontSize', 14)
    
    % spatial location of SAPs (superimposed on smoothed x/y position)
    fig4 = figure(4);
    set(gcf, 'Position', [900 500 560 420])
    hold on
    plot(smooth_headX, smooth_headY, 'k.-')
%     plot(smooth_headX(sap), smooth_headY(sap), 'r*')
    plot(nose_distX, nose_distY, 'r*')
    title(['Location of SAP: ' foldernames{miceiter} '_N1_' num2str(durTotal) 'min'], 'Interpreter', 'none')
    xlim([50 500])
    ylim([0 400])
    set(gca, 'YDir', 'reverse', 'FontSize', 14)
    
    pause
    if(0)
        csvwrite([curr_set '_' foldernames{miceiter} '_N1_' num2str(durTotal) 'min_sap.csv'], sap)
        saveas(fig3, [curr_set '_' foldernames{miceiter} '_N1_' num2str(durTotal) 'min_bodyLen_hist.tif'])
        saveas(fig4, [curr_set '_' foldernames{miceiter} '_N1_' num2str(durTotal) 'min_sapPos.tif'])
    end
    
    cd ../..
end

%% plot SAP from .mat file (already created from TimeStatistic_configural)
clear
close all
clc
pause(0.5)
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/

curr_set   = 'Configural_set';
curr_mouse = 'Arya';
fps=15;
disp(['Mouse: ' curr_mouse])
disp(['fps=' num2str(fps)])

cd(curr_set)
load('Configural_set_SAP_plus_dist_10min_smooth.mat')
curr_mouse_loc = find(strcmp({dist_of_sap.name},curr_mouse));

cd(curr_mouse)
cd('Analyzed_Data_10cm')

for session_iter = 1:length(dist_of_sap(1).sapFrame)
    if(strcmp(curr_mouse, 'Aldehyde') && session_iter==1)
        disp('skip Aldehyde session 1')
        continue
    end
    
    % % histogram of body lengths (sanity check)
    % fig3 = figure(3);
    % set(gcf, 'Position', [300 500 560 420])
    % hold on
    % hist(bodyLen,100)
    % sap_threshold = mean(bodyLen)+std(bodyLen);
    % line([sap_threshold sap_threshold]', (fig3.CurrentAxes.YLim)', 'Color', 'r')
    % title(['Body lengths + sap threshold: ' foldernames{miceiter} '_N1_' num2str(durTotal) 'min'], 'Interpreter', 'none')
    % set(gca, 'FontSize', 14)
    
    curr_posX = dist_of_sap(curr_mouse_loc).nosePosAll{session_iter,1};
    curr_posY = dist_of_sap(curr_mouse_loc).nosePosAll{session_iter,2};
    curr_sapX = dist_of_sap(curr_mouse_loc).nosePosSAP{session_iter,1};
    curr_sapY = dist_of_sap(curr_mouse_loc).nosePosSAP{session_iter,2};
    durTotal  = floor(length(curr_posX)/60/fps);
    
    % spatial location of SAPs (superimposed on smoothed x/y position)
    fig1 = figure(1);
    set(gcf, 'Position', [900 500 560 420])
    hold on
    plot(curr_posX, curr_posY, 'k.-')
    plot(curr_sapX, curr_sapY, 'r*')
    title(['Location of SAP: ' curr_mouse '_S' num2str(session_iter) ...
        '_' num2str(durTotal) 'min'], 'Interpreter', 'none')
    xlim([50 500])
    ylim([0 400])
    set(gca, 'YDir', 'reverse', 'FontSize', 14)
    
%     if(0)
%         csvwrite([curr_set '_' curr_mouse '_N1_' num2str(durTotal) 'min_sap.csv'], sap)
%         saveas(fig3, [curr_set '_' curr_mouse '_N1_' num2str(durTotal) 'min_bodyLen_hist.tif'])
        saveas(fig1, [curr_mouse '_S' num2str(session_iter) '_' ...
            num2str(durTotal) 'min_sapPos.tif'])
%     end
    clf
end
disp('end')

%% concatenate structures of SAP distance from object for all mouse sets
clear
close all
clc

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Configural_set
matFiles = dir('*SAP_plus_dist_10min_incl.mat'); % generated at end of TimeStatistic

dist_of_sap_all = [];
for fileiter = 1:length(matFiles)
    load(matFiles(fileiter).name)
    
    dist_of_sap_all = [dist_of_sap_all dist_of_sap];
end
if(0)
    save('SAP_plus_dist_all_crop_incl.mat','dist_of_sap_all')
end

