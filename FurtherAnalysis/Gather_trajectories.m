% Plot trajectories from one mouse in one image (use plots from Analysis.m

clear
close all
clc

groupfolder = ['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Parents_6OHDA_combine/'];
cd(groupfolder)

run('MiceIndex_Parents.m')
detectCond = cat(1, Mice.novelty);
cond2 = find(detectCond=='l'); %'C'; %'l';
cond1 = find(detectCond=='s'); %'S'; %'s';
names = {Mice.name};

for mouseiter = 1:length(Mice)
    currMouse = Mice(mouseiter).name;
    % currMouse = 'M01_Accra';%'AllMice_N1_tail'; %'AllMice_N1_nose';
    % cd temp
    cd(currMouse);%, '/190305/'])
    % cd allFiles
    cd Analyzed_Data_1obj_10cm_nose/;
    
    trajs = dir('*Trajectory*.tif');%'*Trajectory*.tif'); % 'Heatmap*.tif');
    % trajsOrder = [trajs(cond2'); trajs(cond1')];
    %sessionNames = {'1' '2' '3' '4' '5' '6' '7' '8' 'N1' 'N2' 'N3' 'N4'};
    sessionNames = {'H1' 'H2' 'N1' 'N2' 'N3' 'N4' 'N5' 'N6' 'N7' 'N8' 'N9' 'N10'};
    % sessionNames = [names(cond2) names(cond1)];
    
    fig1 = figure(1);
    % suptitle(currMouse)
    % suptitle([currMouse(1:7), ' ' currMouse(9:end)])
    % suptitle('Trajectories N1: Planets (tail)')
    % set(fig1, 'Position', [200 1 1500 900]); %DRILLS
    % set(fig1, 'Position', [200 1 1860 811]); %NewHope-ROTJ
    set(fig1, 'Position', [472 60 888 841])
    
    suptitle([currMouse(1:3), ' ' currMouse(5:end)])
    % set(fig1, 'Position', [200 1 1500 500]);
    hold on
    
    for i = 1:length(trajs)
        image = imread(trajs(i).name);
        %     image = imread(trajsOrder(i).name);
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
            subplot(3, 3, i+1) %i+1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            imshow(imageCropped, 'InitialMagnification', 50)
            title(sessionNames(i), 'Interpreter', 'none')
        end
    end
    % %
    cd ../
    %saveas(fig1, ['Trajectories: ' currMouse])
    saveas(fig1, ['Trajectories: ' currMouse, '_nose.tif'])
    close all
    disp('Saved')
    cd(groupfolder)
end