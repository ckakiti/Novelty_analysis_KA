clear
clc
close all

Config_NovAna

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/DRILLS_DLC/temp/Heatmap_subtract

files1 = dir('*nose.mat');
files2 = dir('*tail.mat');
files1_name = cat(1,files1.name);
files2_name = cat(1,files2.name);

load(files1_name)
Labels1 = Labels;
clear Labels

load(files2_name)
Labels2 = Labels;
clear Labels

%%


%% Heatmap: nose - tail
% Plot parameters
vspace=3;               % average frame space to calculate the velocity Must be a odd intiger
plot_fs = Dis_ts_frame; % Distance/Orientation plot start and end frame
plot_fe = Dis_te_frame;%min(fpm*10, size(Labels,1));

x_length=video_xlen;   % Heatmap x and y axis length (pixels)
y_length=video_ywid;

% Heatmap
fov1=zeros(x_length,y_length);
fov2=zeros(x_length,y_length);
for i=plot_fs:plot_fe
    if round(Labels1(i,14))<x_length && round(Labels1(i,15))<y_length...
            && round(Labels1(i,14))>0 && round(Labels1(i,15))>0
        fov1(round(Labels1(i,14)),round(Labels1(i,15)))=fov1(round(Labels1(i,14)),round(Labels1(i,15)))+10;
    end
    if round(Labels2(i,14))<x_length && round(Labels2(i,15))<y_length...
            && round(Labels2(i,14))>0 && round(Labels2(i,15))>0
        fov2(round(Labels2(i,14)),round(Labels2(i,15)))=fov2(round(Labels2(i,14)),round(Labels2(i,15)))+10;
    end
end

%pooling method 1
pool_size=5; %5-10 %og=50
pooled_map1=zeros(x_length-pool_size+1,y_length-pool_size+1);
pooled_map2=zeros(x_length-pool_size+1,y_length-pool_size+1);

for i=1:x_length-pool_size+1
    for j=1:y_length-pool_size+1
        pooled_map1(i,j)=sum(sum(fov1(i:i+pool_size-1,j:j+pool_size-1)));
        pooled_map2(i,j)=sum(sum(fov2(i:i+pool_size-1,j:j+pool_size-1)));
    end
end


%%
close all

Hmfigure=figure('visible','on');
hm=heatmap((pooled_map1-pooled_map2)','GridVisible','off','Colormap',parula,'FontSize',0.01);
caxis([-1440 760]/4)

if(0)
    saveas(Hmfigure,'Heatmap_subtract_temp.tif')
end