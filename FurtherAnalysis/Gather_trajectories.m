% Plot trajectories from one mouse in one image (use plots from Analysis.m

clear
close all
clc

currMouse = 'MSFP';

cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/MSFP_Test/' ...
    currMouse])%, '/190305/'])
cd Analyzed_Data/;

trajs = dir('Trajectory*.tif'); % dir('Heatmap*.tif');
%sessionNames = {'1' '2' '3' '4' '5' '6' '7' '8' 'N1' 'N2' 'N3' 'N4'};
%sessionNames = {'H1' 'H2' 'N1' 'N2' 'N3' 'N4' 'N5' 'N6' 'N7' 'N8' 'N9' 'N10'};
sessionNames = {'H1' 'N1' 'N2'};

fig1 = figure(1);
%suptitle([currMouse(1:2), ' ' currMouse(4:end)])
suptitle(currMouse)
set(fig1, 'Position', [200 1 1500 500]);%[200 1 1500 500]); %[200 1 1500 500]);
hold on

for i = 1:length(trajs)
    image = imread(trajs(i).name);
    imageCropped = imcrop(image, [500 150 1050 1050]);%[500 200 950 950]);%[550 150 1150 950]);
    
    if(i<=2)
        subplot(1, 3, i)
        imshow(imageCropped, 'InitialMagnification', 60)
        title(sessionNames(i))
    else
        subplot(1, 3, i) %i+1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        imshow(imageCropped, 'InitialMagnification', 50)
        title(sessionNames(i))
    end
end
%%
cd ../../
%saveas(fig1, ['Trajectories: ' currMouse])
saveas(fig1, ['Trajectories: ' currMouse, '.tif'])
%saveas(fig1, ['Trajectories: ' currMouse '_body_H2N1.tif'])
disp('Saved')