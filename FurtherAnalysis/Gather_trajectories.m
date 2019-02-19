% Plot trajectories from one mouse in one image (use plots from Analysis.m

clear
clc

currMouse = 'S6_Wedge';

cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/NewHope-ROTJ/' ...
    currMouse])% '/181223/'])
cd Analyzed_Data;

trajs = dir('Trajectory*.tif'); % dir('Heatmap*.tif');
%sessionNames = {'1' '2' '3' '4' '5' '6' '7' '8' 'N1' 'N2' 'N3' 'N4'};
sessionNames = {'H1' 'H2' 'N1' 'N2' 'N3' 'N4' 'N5' 'N6' 'N7' 'N8' 'N9' 'N10'};

fig1 = figure(1);
suptitle([currMouse(1:2), ' ' currMouse(4:end)]) % currMouse)
set(fig1, 'Position', [200 1 1500 800]);
hold on

for i = 1:8 %length(trajs)
    image = imread(trajs(i).name);
    imageCropped = imcrop(image, [550 150 1150 950]);
    
    if(i<=2)
        subplot(3, 3, i)
        imshow(imageCropped, 'InitialMagnification', 60)
        title(sessionNames(i))
    else
        subplot(3, 3, i+1) %i) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        imshow(imageCropped, 'InitialMagnification', 50)
        title(sessionNames(i))
    end
end

%saveas(fig1, ['Trajectories: ' currMouse])
%saveas(fig1, ['Trajectories: ' currMouse '_body.tif'])
%disp('Saved')