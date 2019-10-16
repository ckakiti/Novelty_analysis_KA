% networkname_format='DeepCut_resnet50_noveltyMay21shuffle1_700000'; % rubber toy orig
% networkname_format='DeepCut_resnet50_MoSeqNoveltySep12shuffle1_1030000';
networkname_format='DeepCut_resnet50_MoSeqNovelty_RetrainSep17shuffle1_1030000';
% networkname_format='DeepCut_resnet50_MoSeqConfiguralDec04shuffle1_1030000';
% networkname_format='DeepCut_resnet50_MoSeqWithOldArenaDec31shuffle1_1030000';
% networkname_format='DeepCut_resnet50_MoSeqConfiguralUpdateMar04shuffle1_1030000';
% networkname_format='DeepCut_resnet50_MoSeqFPApr24shuffle1_1030000';
% networkname_format='DeepCut_resnet50_MitsukoOperantLearningSep15shuffle1_1030000';
% networkname_format='DeepCut_resnet50_MitsukoOperant_FPJun18shuffle1_1030000';
% networkname_format = 'DeepCut_resnet50_MoSeqGreyMouseJul23shuffle1_1030000';

videoname_format='XYZ_YYMMDD_rgb.avi'; %mp4; %avi';

video_xlen=512;%520;%512;%400
video_ywid=424;%420;%424;%480            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fpm = 900; %5040 %1800 %1500 %900 % frames per minute
fps = fpm./60;          % frames per second
ppc = 42/6.3; %123/23; %340/(2*30);   % pixels per cm   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radius_cm = 8; %10 for tea obj %9 for MoSeq %8 default %6 Mit_FP  % Time spent around the obj radius (cm)   %%%%%%%%%%%%%
radius = radius_cm.*ppc;% Time spent around the obj radius (pixels)
angle_radius = 15;      % Time orient towards the obj radius (degree)

% Calculation Parameters
Dis_ts_frame=500; % Time spent around the obf start and end frame
Dis_te_frame=10.*60.*fps+Dis_ts_frame;

Ang_ts_frame=500; % Time orient towards the obj start and end frame
Ang_te_frame=10.*60.*fps+Ang_ts_frame;