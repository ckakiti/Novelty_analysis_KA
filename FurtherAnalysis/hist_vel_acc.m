clear
close all
clc

% parameters
fps = 15;
ppc = 42/6.3;
cutoff = fps*60*10;

currSet = 'Rim_KO_DLC';
mice = {'01', '02', '03', '04', '05', '06'};

% for plotting
fig1 = figure(1);
set(fig1, 'Position', [46 1 1726 973])
suptitle(['Velocity (avg every ' num2str(fps) ' frames)'])

fig2 = figure(2);
set(fig2, 'Position', [46 1 1726 973])
suptitle(['Acceleration (avg every ' num2str(fps) ' frames)'])
%%
for mousei = 1:length(mice)
    currMouse = mice{mousei};
    currFile  = 2; %H2
    cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' currSet '/' currMouse '/Analyzed_Data_1obj'])
    
    posFile = dir('*rgb.mat');
    DLCPos = load(posFile(currFile).name);
    load('Arena_Obj_Pos.mat');

    
    xPos_body = DLCPos.Labels(1:cutoff,2);
    yPos_body = DLCPos.Labels(1:cutoff,3);
    xPos_smooth = smooth(xPos_body,30);
    yPos_smooth = smooth(yPos_body,30);
    
    % first average x and y position every 15 frames (1s)
    %  then calculate velocity (difference in distance / dt)
    posWin = fps;
    modPos = mod(length(xPos_smooth),posWin);
    xPos_smooth_cut = xPos_smooth(1:end-modPos);
    yPos_smooth_cut = yPos_smooth(1:end-modPos);
    
    reshape_xPos = reshape(xPos_smooth_cut, posWin, []);
    reshape_yPos = reshape(yPos_smooth_cut, posWin, []);
    avg_xPos = mean(reshape_xPos,1);
    avg_yPos = mean(reshape_yPos,1);    
    
    % calculating velocity+acceleration from averaged xPos/yPos
    diff_xPos = diff(avg_xPos);
    diff_yPos = diff(avg_yPos);
    dt = posWin/fps;
    vel = sqrt( (diff_xPos .* diff_xPos) + (diff_yPos .* diff_yPos) )/(dt*ppc); %units = cm/s
    
    vel_smooth = smooth(vel, 4);  
    acc = diff(vel_smooth); %units = cm/s^2
    
    
    % plotting velocity
    figure(1)
    subplot(2,3,mousei)
    histogram(vel_smooth, 0:1:20)
    title(['Mouse' currMouse ' Day' num2str(currFile)])
    xlabel('Velocity (cm/s)')
    xlim([0 20])
    ylim([0 150])
    
    % plotting acceleration
    figure(2)
    subplot(2,3,mousei)
    histogram(acc, -50:1:50)
    title(['Mouse' currMouse ' Day' num2str(currFile)])
    xlabel('Acceleration (cm/s^2)')
    xlim([-15 15])
%     ylim([0 300])
end

%% saving
cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' currSet])

saveas(fig1, ['histVel_' currSet])
saveas(fig1, ['histVel_' currSet '.tif'])
saveas(fig2, ['histAcc_' currSet])
saveas(fig2, ['histAcc_' currSet '.tif'])

close all