clear
close all
clc

% path_to_your_files = ['/home/alex/Programs/DeepLabCut_new/DeepLabCut/'...
%     'videos/Capoeira_DLC/'];
path_to_your_files = ['/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/' ...
   'Behavior/others/Standard setup/CombineAnalysis/Just-in-case files/Hiking_DLC'];
cd(path_to_your_files)

% Config_NovAna_trim
path_to_configFile = ['/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/' ...
   'Behavior/others/Standard setup/CombineAnalysis/Just-in-case files/Code/'];
run([path_to_configFile 'Config_NovAna.m'])

durTotal = 10; % duration of analysis (min)
disp(['Duration of analysis: ' num2str(durTotal) 'min'])

Dis_ts_frame=500;
Dis_te_frame=durTotal.*60.*fps+Dis_ts_frame;
startframe  = Dis_ts_frame;
endframe    = Dis_te_frame;

folderpath = cd;
folderd = dir(folderpath);
isub = [folderd(:).isdir];
foldernames = {folderd(isub).name}'; 
foldernames(ismember(foldernames,{'.','..','Retrain_Sep17','temp'})) = []; 
folderlen=length(foldernames);

% smoothing window for orientation+SAP analysis
Swindow=30;
Swindow_orient=40;

% maximum body length (pixels), anything above counted as error
bodyLen_cutoff = 80;

disp('next')

%%
nearOrSAP = input('Analyzing time near obj or SAP? 0=timeNear 1=SAP: ');
if(nearOrSAP==0)
    whichAnalyzedFolder = 'Analyzed_Data_1obj_head';
    disp('Analyzing time near / orientation to obj')
    disp(['radius_cm = ' num2str(radius_cm)])
    
    normOrNo = [];
elseif(nearOrSAP==1)
    whichAnalyzedFolder = 'Analyzed_Data_1obj_head'; %8cm=head 12cm=tail
    disp('Analyzing SAP number')
    disp(['radius_cm = ' num2str(radius_cm)])
    
    normOrNo = input('Normalize sapNear by how much time mice are in radius? 0/1: ');
    if(normOrNo==1)
        disp('normalize sapNear')
    else
        disp('no normalization')
    end
else
    error_msg = 'please type 0 or 1';
    error(error_msg)
end

% pause(1)

for folderi=1:folderlen
    cd(foldernames{folderi});
    disp(foldernames{folderi})
    dist_of_sap(folderi).name = foldernames{folderi};
    
    cd(whichAnalyzedFolder)
    load('Arena_Obj_Pos.mat','arena')

    subpath=cd;
    PathRoot=[subpath '/'];
    filelist=dir('*rgb_Converted.mat'); 
    flen = length(filelist);
    
    for filei = 3%1:flen %3
        disp(filei)
        filename = filelist(filei).name;
        load(filename)
        
        % make sure end frame doesn't exceed current Labels size
        endframe = min(Dis_te_frame, size(Labels,1));
        LabelsCut=Labels(startframe:endframe,:);
        
        % Analyze time spent near obj (within specified radius)
        Time_distance(filei) = sum(LabelsCut(:,21));
        
        
        % Analyze time spent oriented to obj (after smoothing)
        orientSmooth = smoothdata(LabelsCut(:,22),'rloess',Swindow_orient);
        isOrient = intersect(find(orientSmooth<=angle_radius), find(orientSmooth>=(-angle_radius)));
        Time_angle(filei) = length(isOrient);
        
        % analyze orientation across session (just for N1)
        orientAllSmooth = smoothdata(Labels(1:(15*60*25),22),'rloess',Swindow_orient);
        isOrientAll = (orientAllSmooth<=angle_radius) & (orientAllSmooth>=(-angle_radius));
        orient_cumsum = isOrientAll;%cumsum(isOrientAll);

        % Analyze body length and SAP
        smooth_headX = smoothdata( LabelsCut(:,2), 'rloess', Swindow);
        smooth_headY = smoothdata( LabelsCut(:,3), 'rloess', Swindow);
        smooth_tailX = smoothdata( LabelsCut(:,11), 'rloess', Swindow);
        smooth_tailY = smoothdata( LabelsCut(:,12), 'rloess', Swindow);
        curr_velX  = LabelsCut(:,24);
        curr_velY  = LabelsCut(:,25);

        bodyLen = sqrt( (smooth_tailX-smooth_headX).^2 + (smooth_tailY-smooth_headY).^2 ); % pixels
        bodyLen(bodyLen>bodyLen_cutoff)=0; % body lengths >80px are likely to be errors

        vel        = sqrt( curr_velX.^2 + curr_velY.^2 ); % head velocity, cm/s
        stretches  = find(bodyLen>(mean(bodyLen)+std(bodyLen)));
        slows      = find(vel<5);
        sap        = intersect(stretches,slows);
        sap_in_rad = intersect(sap, find(LabelsCut(:,21)==1));
        
%         diffSap      = diff(sap); 
%         diffSapInRad = diff(sap_in_rad);
%         sap_frame{filei,1}    = sap(find(diffSap>1)+1); % only want to count a SAP once
%         sap_num(filei)        = sum(diffSap>1); % total SAPs (any location)
%         sap_dist{filei,1}     = LabelsCut(diffSap>1,17);
%         nose_dist{filei,1}    = smooth_headX(diffSap>1);
%         nose_dist{filei,2}    = smooth_headY(diffSap>1);
        sap_frame{filei,1}    = sap;
        sap_num(filei)        = length(sap); % total SAPs (any location)
        sap_dist{filei,1}     = LabelsCut(sap,17);
        sap_pos{filei,1}      = smooth_headX(sap);
        sap_pos{filei,2}      = smooth_headY(sap);
        nose_dist{filei,1}    = smooth_headX;
        nose_dist{filei,2}    = smooth_headY;
        
        if(normOrNo) % normalized
%             sap_num_in_rad(filei) = sum(diffSapInRad>1)/Time_distance(filei); 
            sap_num_in_rad(filei) = length(sap_in_rad)/Time_distance(filei); % (norm)
        else         % not normalized
%             sap_num_in_rad(filei) = sum(diffSapInRad>1);
            sap_num_in_rad(filei) = length(sap_in_rad);
        end
    end
    
    Time_distance_all(folderi,:)  = Time_distance;
    Time_angle_all(folderi,:)     = Time_angle;
    Time_angle_cumul(folderi,:)   = orient_cumsum';
    SAP_num_all(folderi,:)        = sap_num;
    SAP_in_rad_all(folderi,:)     = sap_num_in_rad;
    dist_of_sap(folderi).sapFrame = sap_frame;
    dist_of_sap(folderi).sapDist  = sap_dist;
    dist_of_sap(folderi).sapPos   = sap_pos;
    dist_of_sap(folderi).nosePos  = nose_dist;
    
    clearvars Time_angle Time_distance ...
        sap_num sap_num_in_rad sap_frame sap_dist sap_pos nose_dist isOrientAll
    cd ../..
    DisOrAng{folderi,1}='Distance';
    DisOrAng{folderlen+folderi,1}='Angle';
    SAPNumNear{folderi,1}='SAPnum';
    SAPNumNear{folderlen+folderi,1}='SAPnear';

end

Time_distance_all = Time_distance_all./double(endframe-startframe);
Time_angle_all    = Time_angle_all./double(endframe-startframe);

Table  = table(cat(1,foldernames,foldernames),DisOrAng,[Time_distance_all;Time_angle_all]);
Table2 = table(cat(1,foldernames,foldernames),SAPNumNear,[SAP_num_all;SAP_in_rad_all]);

if(0)
    writetable(Table,'TimeStatistic_tail.csv');
    writetable(Table2,'TimeStatistic_SAP_8cm_norm.csv');
    
    save('TimeStatistic_orient_cumul', 'foldernames','Time_angle_cumul');
    save('SAP_plus_dist.mat', 'dist_of_sap')
    
    disp('saved')
end

disp('end')