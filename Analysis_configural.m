%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapting Analysis.m code to a configural novelty setup
% CKA 181207
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All the information are saved as folowing
% Labels(:,2) Nose x (pixel)
% Labels(:,3) Nose y (pixel)
% Labels(:,5) Leftear x (pixel)
% Labels(:,6) Leftear y (pixel)
% Labels(:,8) Rightear x (pixel)
% Labels(:,9) Rightear y (pixel)
% Labels(:,11) Tailbase x (pixel)
% Labels(:,12) Tailbase y (pixel)
% Labels(:,14) Head x (pixel) 'average of nose, leftear and rightear'
                             % for center of mass, use average of nose and tail 
% Labels(:,15) Head y (pixel)
% Labels(:,17) Head distance from object 1 (cm)
% Labels(:,18) Head distance from object 2
% Labels(:,19) Head distance from object 3
% Labels(:,20) Head distance from object 4
% Labels(:,21) if distance to object 1 <= radius 1; else 0
% Labels(:,22) if distance to object 2 <= radius 1; else 0
% Labels(:,23) if distance to object 3 <= radius 1; else 0
% Labels(:,24) if distance to object 4 <= radius 1; else 0
% Labels(:,25)
% Labels(:,26)
% Labels(:,27)

% Dis_t_obj_x Time spent in the radius for obj x according to the distance

%***********************************************************
% Initialization
%***********************************************************
clear
clc
cd /home/alex/Programs/Novelty_analysis_KA
Config_NovAna;

% radius_cm = 10; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% radius = radius_cm.*ppc; %%%%%%%%%%%%%%%%%%%%%%%%%

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Rim_KO_DLC/03
pathname = cd;
PathRoot=[pathname '/'];
filelist=dir([PathRoot,'*' videoname_format(end-3:end)]);
flen = length(filelist);
cd Analyzed_Data_4obj;

if isfile('Arena_Obj_Pos_4obj.mat')
    load('Arena_Obj_Pos_4obj.mat', 'obj_center', 'obj', 'arena');
    obj_temp = cell(flen, 1);
    
    for replaceObj = 1:flen
       obj_temp(replaceObj,1) = {obj(replaceObj,:)}; 
    end
    obj = obj_temp;
    
else
    load('Arena_Obj_Pos.mat');
end

cd ..

for fiter =1:flen
    if ~isempty(strfind(filelist(fiter).name,'abeled'))
        filelist(fiter)=[];
    end
end
flen = length(filelist);

%%
for fiter =1:flen%%%%%
    vn = filelist(fiter).name;
    fn=[vn(1:end-4) networkname_format '.csv'];
    disp(['Analyzing: ' fn]);

    Labels = csvread(fn,3,0);
    len = length(Labels(:,1));
    Labels = [Labels zeros(len,14)];

    %***********************************************************
    % Parameters
    %***********************************************************
    vspace=3;      % average frame space to calculate the velocity Must be a odd intiger
    % Plot parameters
    plot_fs = 1;      % Distance/Orientation plot start and end frame
    plot_fe = fpm*10;
    if(plot_fe>len)
        plot_fe = len;
    end
    
    x_length=video_xlen;   % Heatmap x and y axis length (pixels)
    y_length=video_ywid;

    %***********************************************************
    % Calculation
    %***********************************************************

    % Calculate head position
    Labels(:,14)=Labels(:,2); %Labels(:,5)+Labels(:,8))./3;
    Labels(:,15)=Labels(:,3); %Labels(:,6)+Labels(:,9))./3;
    %Labels(:,16)=(Labels(:,4)+Labels(:,7)+Labels(:,10))./3;

    % head distance from object center
    Labels(:,17)=sqrt((obj_center{fiter,1}(1)-Labels(:,14)).^2+(obj_center{fiter,1}(2)-Labels(:,15)).^2)/ppc;
    Labels(:,18)=sqrt((obj_center{fiter,2}(1)-Labels(:,14)).^2+(obj_center{fiter,2}(2)-Labels(:,15)).^2)/ppc;
    Labels(:,19)=sqrt((obj_center{fiter,3}(1)-Labels(:,14)).^2+(obj_center{fiter,3}(2)-Labels(:,15)).^2)/ppc;
    Labels(:,20)=sqrt((obj_center{fiter,4}(1)-Labels(:,14)).^2+(obj_center{fiter,4}(2)-Labels(:,15)).^2)/ppc;

    % Time spent near the obj
    for i=1:len
        if Labels(i,17)<=radius_cm
            Labels(i,21)=1;
        else
            Labels(i,21)=0;
        end
        
        if Labels(i,18)<=radius_cm
            Labels(i,22)=1;
        else
            Labels(i,22)=0;
        end
        
        if Labels(i,19)<=radius_cm
            Labels(i,23)=1;
        else
            Labels(i,23)=0;
        end
        
        if Labels(i,20)<=radius_cm
            Labels(i,24)=1;
        else
            Labels(i,24)=0;
        end
    end

    if(Dis_te_frame>len)
        Dis_te_frame = len;
    end
    
    Dis_t_obj_1 = sum(Labels(Dis_ts_frame:Dis_te_frame,21))./(Dis_te_frame-Dis_ts_frame);
    Dis_t_obj_2 = sum(Labels(Dis_ts_frame:Dis_te_frame,22))./(Dis_te_frame-Dis_ts_frame);
    Dis_t_obj_3 = sum(Labels(Dis_ts_frame:Dis_te_frame,23))./(Dis_te_frame-Dis_ts_frame);
    Dis_t_obj_4 = sum(Labels(Dis_ts_frame:Dis_te_frame,24))./(Dis_te_frame-Dis_ts_frame);
    
    %***********************************************************
    % Plot
    %***********************************************************

    % plot distances
%     Disfigure=figure('visible','off');
%     plot(Labels(plot_fs:plot_fe,1)./fpm,Labels(plot_fs:plot_fe,17),'linewidth',2);
%     hold on
%     plot(Labels(plot_fs:plot_fe,1)./fpm,Labels(plot_fs:plot_fe,18),'linewidth',2);
%     hold on
%     plot(Labels(plot_fs:plot_fe,1)./fpm,Labels(plot_fs:plot_fe,19),'linewidth',2);
%     hold on
%     plot(Labels(plot_fs:plot_fe,1)./fpm,Labels(plot_fs:plot_fe,20),'linewidth',2);
%     title(['Distance (first ' num2str((plot_fe-plot_fs)/fpm) 'min) radius=' num2str(radius_cm) ' cm']);
%     xlabel('time (min)')
%     ylabel('Distance (cm)')
%     set(Disfigure, 'position', [0 0 1000 500]);
            
    % plot orientation
%     Angfigure=figure('visible','off');
%     plot(Labels(plot_fs:plot_fe,1)./fpm,Labels(plot_fs:plot_fe,22),'linewidth',2);
%     title(['Orientation (first ' num2str((plot_fe-plot_fs)/fpm) 'min) radius=' num2str(angle_radius)]);
%     xlabel('time min')
%     ylabel('degree');
%     set(Angfigure, 'position', [0 0 1000 500]);

    % Heatmap
    fov=zeros(x_length,y_length);
    for i=plot_fs:plot_fe
        if round(Labels(i,14))<x_length && round(Labels(i,15))<y_length...
            && round(Labels(i,14))>0 && round(Labels(i,15))>0
            fov(round(Labels(i,14)),round(Labels(i,15)))=fov(round(Labels(i,14)),round(Labels(i,15)))+10;
        end
    end


    %pooling method 1
    pool_size=50;
    pooled_map=zeros(x_length-pool_size+1,y_length-pool_size+1);

    for i=1:x_length-pool_size+1
        for j=1:y_length-pool_size+1
            pooled_map(i,j)=sum(sum(fov(i:i+pool_size-1,j:j+pool_size-1)));
        end
    end

    Hmfigure=figure('visible','off');
    hm=heatmap(pooled_map','GridVisible','off','Colormap',parula,'FontSize',0.01);
    % rectangle('Position',[arena(fiter,1)-round(pool_size/2),...
    %                       arena(fiter,2)-round(pool_size/2),...
    %                       arena(fiter,3)-arena(fiter,1),...
    %                       arena(fiter,4)-arena(fiter,2)])

    
    % Plot trajectory
    Trafigure=figure('visible','off');
    scatter(Labels(plot_fs:plot_fe,14),Labels(plot_fs:plot_fe,15),8,'filled');
    rectangle('Position',[arena(fiter,1),arena(fiter,2),...
                          arena(fiter,3)-arena(fiter,1),...
                          arena(fiter,4)-arena(fiter,2)],'EdgeColor','r','linewidth',4)
    for numOfObj = 1%:4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rectangle('Position', [obj{fiter,numOfObj}(1), obj{fiter,numOfObj}(2), ...
                               obj{fiter,numOfObj}(3)-obj{fiter,numOfObj}(1), ...
                               obj{fiter,numOfObj}(4)-obj{fiter,numOfObj}(2)],...
            'EdgeColor','r','linewidth',4)
    end
    
    hold on
    th = 0:pi/50:2*pi;
    for numOfObj = 1:4
        x  = obj_center{fiter,numOfObj}(1);
        y  = obj_center{fiter,numOfObj}(2);
        xunit = radius * cos(th) + x;
        yunit = radius * sin(th) + y;
        plot(xunit, yunit)
    end
    
    set(gca,'ydir','reverse')
    title(['Trajectory ' vn(1:end-8) ': '...
        num2str(round((plot_fe-plot_fs)/fpm)) 'min, radius=' num2str(radius_cm) 'cm'],...
        'Interpreter', 'none');
    xlim([0 video_xlen+50]);
    ylim([0 video_ywid]);
    set(Trafigure, 'position', [0 0 1200 900]);
    
    % ***********************************************************
    % Save
    % ***********************************************************
    % pause
    cd Analyzed_Data_4obj

%     mkdir([vn(1:end-4) '_Plots'])
%     cd([vn(1:end-4) '_Plots'])

%     saveas(Disfigure,['Distance_' vn(1:end-4) '.png'])
%     saveas(Angfigure,['Orientation_' vn(1:end-4) '.png'])
    saveas(Hmfigure,['Heatmap_' vn(1:end-4) '.tif'])
    saveas(Trafigure,['Trajectory_' vn(1:end-4) '.tif'])

%     cd ..

    save(vn(1:end-4),'Labels','Dis_t_obj_1','Dis_t_obj_2','Dis_t_obj_3','Dis_t_obj_4','radius','radius_cm');
    close all
    clearvars -except arena obj obj_center filelist fiter fpm fps ppc radius radius_cm angle_radius...
                      Dis_ts_frame Dis_te_frame Ang_ts_frame Ang_te_frame video_xlen video_ywid...
                      networkname_format videoname_format

    cd ..
    
end

close all
clear
cd /home/alex/Programs/Novelty_analysis_KA

disp('Done analyzing')