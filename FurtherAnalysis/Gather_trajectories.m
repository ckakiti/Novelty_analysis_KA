% Plot trajectories from one mouse in one image (use plots from Analysis.m

clear
close all
clc


currMouse = 'AllMice_N1_tail'; %'AllMice_N1_nose';

cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Planets_DLC/'])

run('MiceIndex.m')
detectCond = cat(1, Mice.novelty);
cond2 = find(detectCond=='C');
cond1 = find(detectCond=='S');
names = {Mice.name};

cd temp
cd(currMouse);%, '/190305/'])
% cd allFiles
% cd Analyzed_Data_1obj_12cm_tail/;

trajs = dir('*Trajectory*.tif');%'*Trajectory*.tif'); % 'Heatmap*.tif');
trajsOrder = [trajs(cond2'); trajs(cond1')];
%sessionNames = {'1' '2' '3' '4' '5' '6' '7' '8' 'N1' 'N2' 'N3' 'N4'};
% sessionNames = {'H1' 'H2' 'N1' 'N2' 'N3' 'N4' 'N5' 'N6' 'N7' 'N8' 'N9' 'N10'};
sessionNames = [names(cond2) names(cond1)];

fig1 = figure(1);
% suptitle(currMouse)
% suptitle([currMouse(1:7), ' ' currMouse(9:end)])
suptitle('Trajectories N1: Planets (tail)')
% set(fig1, 'Position', [200 1 1500 900]); %DRILLS
set(fig1, 'Position', [200 1 1860 811]); %NewHope-ROTJ

% suptitle([currMouse(1:2), ' ' currMouse(4:end)])
% set(fig1, 'Position', [200 1 1500 500]);
hold on

for i = 1:length(trajs)
%     image = imread(trajs(i).name);
    image = imread(trajsOrder(i).name);
    imageCropped = imcrop(image, [500 200 950 950]); % for novelty arena    
%     imageCropped = imcrop(image, [650 400 750 850]); % for Mitsuko arena
%     imageCropped = imcrop(image, [550 150 950 1050]);  % NewHope-ROTJ(tail)
%     imageCropped = imcrop(image, [600 150 1150 900]);  % NewHope-ROTJ(head)
    
    if(i<=2)
        subplot(3, 3, i)
        imshow(imageCropped, 'InitialMagnification', 60)
%         title(sessionNames(i))
        title(sessionNames(i), 'Interpreter', 'none')
    else
        subplot(3, 3, i) %i+1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        imshow(imageCropped, 'InitialMagnification', 50)
        title(sessionNames(i), 'Interpreter', 'none')
    end
end
%%
cd ../
%saveas(fig1, ['Trajectories: ' currMouse])
saveas(fig1, ['Trajectories: ' currMouse, '.tif'])
close all
disp('Saved')