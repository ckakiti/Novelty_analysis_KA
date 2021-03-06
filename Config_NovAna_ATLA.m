networkname_format='DeepCut_resnet50_MoSeqFP_retrainMay14shuffle1_1030000';
videoname_format='XYZ_YYMMDD_rgb.avi'; %mp4; %avi';

video_xlen=512;%520;%512;
video_ywid=424;%420;%424;            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fpm = 900; %1800 %1500 %900 % frames per minute
fps = fpm./60;          % frames per second
ppc = 42/6.3; %123/23; %340/(2*30);   % pixels per cm   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radius_cm = 7; %10 for tea obj  % Time spent around the obj radius (cm)   %%%%%%%%%%%%%
radius = radius_cm.*ppc;% Time spent around the obj radius (pixels)
angle_radius = 15;      % Time orient towards the obj radius (degree)

% Calculation Parameters
Dis_ts_frame=500;          % Time spent around the obf start and end frame
Dis_te_frame=10.*60.*fps+Dis_ts_frame;

Ang_ts_frame=500;         % Time orient towards the obj start and end frame
Ang_te_frame=10.*60.*fps+Ang_ts_frame;


