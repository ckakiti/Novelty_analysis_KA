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
poke_file = dir('*poke_labels_N1.csv');
pokes     = csvread(poke_file.name);

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
    
    mat_files = dir('*Converted.mat');
    load(mat_files(curr_day).name)
    
    curr_headX = smoothdata( Labels(startframe:endframe,2), 'rloess', Swindow);
    curr_headY = smoothdata( Labels(startframe:endframe,3), 'rloess', Swindow);
    curr_tailX = smoothdata( Labels(startframe:endframe,11), 'rloess', Swindow);
    curr_tailY = smoothdata( Labels(startframe:endframe,12), 'rloess', Swindow);
    curr_velX  = Labels(startframe:endframe,24);
    curr_velY  = Labels(startframe:endframe,25);
    
    bodyLen    = sqrt( (curr_tailX-curr_headX).^2 + (curr_tailY-curr_headY).^2 ); % pixels
    bodyLen(bodyLen>bodyLen_cutoff)=0; % body lengths >80px are likely to be errors
    
    vel       = sqrt( curr_velX.^2 + curr_velY.^2 ); % head velocity, cm/s
    stretches = find(bodyLen>(mean(bodyLen)+std(bodyLen)));
    slows     = find(vel<5);
    sap       = intersect(stretches,slows);
    
    curr_pokes = find(pokes(:,1)==miceiter);
    
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
    plot(curr_headX, curr_headY, 'k.-')
    plot(curr_headX(sap), curr_headY(sap), 'r*')
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

%% concatenate structures of SAP distance from object for all mouse sets
clear
close all
clc

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/StandardSetup_combine
matFiles = dir('*SAP_plus_dist.mat');

dist_of_sap_all = [];
for fileiter = 1:length(matFiles)
    load(matFiles(fileiter).name)
    
    dist_of_sap_all = [dist_of_sap_all dist_of_sap];
end
if(0)
    save('SAP_plus_dist_all.mat','dist_of_sap_all')
end
