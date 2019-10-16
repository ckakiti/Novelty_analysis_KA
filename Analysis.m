%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pending: 
% straight line metric
% pooling method 2 (current pooling method is not ideal)
% the calculation method for derivative needs more thinking because 
% Labels(:,27) Delta Velosity angle (degree/s)
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
% Labels(:,17) Head distance from object (cm)
% Labels(:,18) upper left corner distance from object center  (cm)
% Labels(:,19) average buttom left or upper right conerner distance from object center  (cm)
% Labels(:,20) buttom right corner distance from object center  (cm)
% Labels(:,21) if distance to object <= radius 1; else 0
% Labels(:,22) orientation related to object [-180,180]
% Labels(:,23) if orientation related to object <= +angle_radius and >= -angle_radius 1; else 0
% Labels(:,24) x velocity (cm/s)
% Labels(:,25) y velocity (cm/s)
% Labels(:,26) velosity angle (degree)

% Dis_t_obj Time spent in the radius according to the distance
% Ang_t_obj Time spent orienting towards the object according to the orientation


%***********************************************************
% Initialization
%***********************************************************
clear
clc
cd /home/alex/Programs/Novelty_analysis_KA
Config_NovAna;

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Machines_DLC/Wheel/
%cd(['/media/alex/DataDrive1/MoSeqData/Iku_photometry2/Iku_photometry2_MoSeq/'
% 'Nashville/190425/session_20190425162005/proc'])
cd Analyzed_Data_1obj;
load('Arena_Obj_Pos.mat');
cd ..
pathname = cd;
PathRoot=[pathname '/'];
filelist=dir([PathRoot,'*' videoname_format(end-3:end)]);
flen = length(filelist);

for fiter =1:flen
    if(flen<=length(filelist))
        if ~isempty(strfind(filelist(fiter).name,'abeled'))
            filelist(fiter)=[];
            fiter=fiter-1;
        end
    end
end
flen = length(filelist);

tic;
for fiter = 1:flen %%%%%%%% 
    vn = filelist(fiter).name;
    fn=[vn(1:end-4) networkname_format '.csv'];
    disp(['Analyzing: ' fn]);

    Labels = csvread(fn,3,0);
    len = length(Labels(:,1));
    Labels = [Labels zeros(len,14)];
    Dis_te_frame = min(Dis_te_frame, size(Labels,1));
    Ang_te_frame = min(Ang_te_frame, size(Labels,1));

    %***********************************************************
    % Parameters
    %***********************************************************
    vspace=3;      % average frame space to calculate the velocity Must be a odd intiger
    % Plot parameters
    plot_fs = Dis_ts_frame;      % Distance/Orientation plot start and end frame
    plot_fe = Dis_te_frame;%min(fpm*10, size(Labels,1));

    x_length=video_xlen;   % Heatmap x and y axis length (pixels)
    y_length=video_ywid;

    %***********************************************************
    % Calculation
    %***********************************************************

    % Calculate head or body position
    Labels(:,14)=(Labels(:,2));%+Labels(:,5)+Labels(:,8))./3; %Labels(:,11))./2;
    Labels(:,15)=(Labels(:,3));%+Labels(:,6)+Labels(:,9))./3; %Labels(:,12))./2;
%    Labels(:,16)=(Labels(:,4)+Labels(:,7)+Labels(:,10))./3;

    % head distance from object center
    Labels(:,17)=sqrt((obj_center(fiter,1)-Labels(:,14)).^2+(obj_center(fiter,2)-Labels(:,15)).^2)/ppc;

    % upper left corner distance from object center
    Labels(:,18)=sqrt((obj_center(fiter,1)-arena(fiter,1)).^2+(obj_center(fiter,2)-arena(fiter,2)).^2)/ppc;
    % bottom left or upper right corner distance from object center
    Labels(:,19)=0.5.*(sqrt((obj_center(fiter,1)-arena(fiter,1)).^2+(obj_center(fiter,2)-arena(fiter,4)).^2)...
                    + sqrt((obj_center(fiter,1)-arena(fiter,3)).^2+(obj_center(fiter,2)-arena(fiter,2)).^2))/ppc;                   
    % bottom right corner distance from object center
    Labels(:,20)=sqrt((obj_center(fiter,1)-arena(fiter,3)).^2+(obj_center(fiter,2)-arena(fiter,4)).^2)/ppc;


    % Time spent near the obj
    for i=1:len
        if Labels(i,17)<=radius_cm
            Labels(i,21)=1;
        else
            Labels(i,21)=0;
        end
    end

    Dis_t_obj = sum(Labels(Dis_ts_frame:Dis_te_frame,21))./(Dis_te_frame-Dis_ts_frame);

    % Calculate Orientation  v as vector h2t head to tail o2t obj to tail
    v_h2t = [(Labels(:,2)-Labels(:,11)),(Labels(:,3)-Labels(:,12))];
    v_o2t = [(obj_center(fiter,1)-Labels(:,11)),(obj_center(fiter,2)-Labels(:,12))];

    for i=1:len
        Labels(i,22)= atan2d(det([v_h2t(i,:);v_o2t(i,:)]),dot(v_h2t(i,:),v_o2t(i,:)));
    end

    for i=1:len
        if Labels(i,22)<=angle_radius && Labels(i,22)>=(-angle_radius)
            Labels(i,23)=1;
        else
            Labels(i,23)=0;
        end
    end

    Ang_t_obj = sum(Labels(Ang_ts_frame:Ang_te_frame,23))./(Ang_te_frame-Ang_ts_frame);


    % Calculate head velocity
    floorfra = floor(vspace/2);
    for i = floorfra+1:len-vspace-floorfra
        Labels(i,24) = ((sum(Labels(i-floorfra:i+floorfra,14))-sum(Labels(i-floorfra+vspace:i+floorfra+vspace,14)))./ppc)./(vspace./fps);
        Labels(i,25) = ((sum(Labels(i-floorfra:i+floorfra,15))-sum(Labels(i-floorfra+vspace:i+floorfra+vspace,15)))./ppc)./(vspace./fps);
    end
    Labels(:,26)=atan2d(Labels(:,25),Labels(:,24));

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
%     Trafigure=figure('visible','on');
    scatter(Labels(plot_fs:plot_fe,14),Labels(plot_fs:plot_fe,15),6,'filled');
%     plot(Labels(plot_fs:plot_fe,14),Labels(plot_fs:plot_fe,15),'k')
    rectangle('Position',[arena(fiter,1),arena(fiter,2),...
                          arena(fiter,3)-arena(fiter,1),...
                          arena(fiter,4)-arena(fiter,2)],'EdgeColor','r','linewidth',8)
    rectangle('Position',[obj(fiter,1),obj(fiter,2),...
                          obj(fiter,3)-obj(fiter,1),obj(fiter,4)-obj(fiter,2)],...
                          'EdgeColor','r','linewidth',2)
    
    hold on
    th = 0:pi/50:2*pi;
    x  = obj_center(fiter,1);
    y  = obj_center(fiter,2);
    xunit = radius * cos(th) + x;
    yunit = radius * sin(th) + y;
    plot(xunit, yunit,'r--','linewidth',3)
    
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
    cd Analyzed_Data_1obj

%     mkdir([vn(1:end-4) '_Plots'])
%     cd([vn(1:end-4) '_Plots'])

%     saveas(Disfigure,['Distance_' vn(1:end-4) '.png'])
%     saveas(Angfigure,['Orientation_' vn(1:end-4) '.png'])
    %cd('./body') %%%
    saveas(Hmfigure,['Heatmap_' vn(1:end-4) '.tif'])
    saveas(Trafigure,['Trajectory_' vn(1:end-4) '.tif'])

%     cd ..

    save(vn(1:end-4),'Labels','Dis_t_obj','Ang_t_obj', 'radius', 'radius_cm');
    
    %cd .. %%%
    
    close all
    clearvars -except arena obj obj_center filelist fiter fpm fps ppc radius radius_cm angle_radius...
                      Dis_ts_frame Dis_te_frame Ang_ts_frame Ang_te_frame video_xlen video_ywid...
                      networkname_format videoname_format

    cd ..
    toc;
end

close all
clear
cd /home/alex/Programs/Novelty_analysis_KA

disp('end')