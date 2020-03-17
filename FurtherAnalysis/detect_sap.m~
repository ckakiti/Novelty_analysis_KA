%%
clear
close all
clc

Config_NovAna
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Capoeira_DLC/

poke_file = dir('*poke_labels_N1.csv');
pokes     = csvread(poke_file.name);


curr_day = 3;  % N1
durTotal = 10; % duration of analysis (min)
disp(['Duration of analysis: ' num2str(durTotal) 'min'])

Dis_ts_frame = 500;
Dis_te_frame = durTotal.*60.*fps+Dis_ts_frame;
startframe   = 1;    %Dis_ts_frame;
endframe     = 9000; %Dis_te_frame;


folderpath = cd;
folderd = dir(folderpath);
isub = [folderd(:).isdir];
foldernames = {folderd(isub).name}'; 
foldernames(ismember(foldernames,{'.','..','Retrain_Sep17','temp'})) = []; 
folderlen=length(foldernames);


disp('next')

%%
close all
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Capoeira_DLC/

% for smoothing
Swindow=30;

for miceiter = 1%1:folderlen
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
    bodyLen(bodyLen>80)=0; % body lengths >80px are likely to be errors
    
    vel       = sqrt( curr_velX.^2 + curr_velY.^2 ); % head velocity, cm/s
    stretches = find(bodyLen>(mean(bodyLen)+std(bodyLen)));
    slows     = find(vel<5);
    sap       = intersect(stretches,slows);
    
    curr_pokes = find(pokes(:,1)==miceiter);
    
    % plot trajectories for nose and tail side-by-side
    fig1 = figure(1);
    set(gcf, 'Position', [100 300 1100 450])
    subplot(1,2,1)
    plot(curr_headX, curr_headY, 'k.-')
    axis square
    
    subplot(1,2,2)
    plot(curr_tailX, curr_tailY, 'k.-')
    axis square
    
    % plot body length (distance between nose and tail), superimpose automatically calculated poke 
    fig2 = figure(2);
    hold on
    plot(bodyLen, 'k.-')
    line([pokes(curr_pokes,2) pokes(curr_pokes,2)]', repmat([0 max(bodyLen(:))], length(curr_pokes), 1)')
    
    % histogram of body lengths (sanity check)
    figure(3)
    hist(bodyLen,100)
    
    pause
    cd ../..
end


