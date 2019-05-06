%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

Config_NovAna

AnalysisDay=3; % 3 = first novelty day
fps=30; % 15 rgb; % 25; 30 depth; 

startframe=1;%1;%900;
endframe=(fps*60*10)+startframe;%15000;%18000;
Swindow=40;         % Smooth Window Size %40;
DisThreshold=radius;   % Distance threshold % 10; 60;

mouseSet = 'Iku_photometry2_MoSeq';
%cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' mouseSet])
cd(['/media/alex/DataDrive1/MoSeqData/Iku_photometry2/' mouseSet])
Mice(1).name='Nashville';
Mice(2).name='Omaha';
Mice(3).name='Rochester';
Mice(4).name='Syosset';
Mice(5).name='Universal';
Mice(6).name='Vegas';
Mice(1).novelty='S';
Mice(2).novelty='S';
Mice(3).novelty='S';
Mice(4).novelty='S';
Mice(5).novelty='S';
Mice(6).novelty='S';
Mice(1).datanum=[];
Mice(2).datanum=[];
Mice(3).datanum=[];
Mice(4).datanum=[];
Mice(5).datanum=[];
Mice(6).datanum=[];

%%
%mouseSet = 'CvsS_180831_DLC';
%mouseSet = '7day_preexposure_combine';
mouseSet = 'Iku_photometry_DLC';
%mouseSet = 'MSFP_Test';
%mouseSet = 'FP_Test';

if(strcmp(mouseSet, 'CvsS_180831_DLC'))
    cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/CvsS_180831_DLC/
    Mice(1).name='C4';
    Mice(2).name='C5';
    Mice(3).name='C6';
    Mice(4).name='C7';
    Mice(5).name='S4';
    Mice(6).name='S5';
    Mice(7).name='Mal';
    Mice(8).name='Wash';
    Mice(1).novelty='C';
    Mice(2).novelty='C';
    Mice(3).novelty='C';
    Mice(4).novelty='C';
    Mice(5).novelty='S';
    Mice(6).novelty='S';
    Mice(7).novelty='S';
    Mice(8).novelty='S';
    Mice(1).datanum=71;
    Mice(2).datanum=68;
    Mice(3).datanum=40;
    Mice(4).datanum=37;
    Mice(5).datanum=35;
    Mice(6).datanum=39;
    Mice(7).datanum=46;
    Mice(8).datanum=82;
    %AllActLabels=csvread('CvsS_poke_labels_N1_byHand.csv',1,2);
    
elseif(strcmp(mouseSet, '7day_preexposure_combine'))
    Mice(1).name='C1_Aldehyde';
    Mice(2).name='C2_Ester';
    Mice(3).name='C3_Thiol';
    Mice(4).name='C4_George';
    Mice(5).name='C5_Hermione';
    Mice(6).name='C6_Ron';
    Mice(7).name='S1_Alcohol';
    Mice(8).name='S2_Amine';
    Mice(9).name='S3_Ketone';
    Mice(10).name='S4_Fred';
    Mice(11).name='S5_Harry';
    Mice(12).name='S6_Neville';
    Mice(1).novelty='C';
    Mice(2).novelty='C';
    Mice(3).novelty='C';
    Mice(4).novelty='C';
    Mice(5).novelty='C';
    Mice(6).novelty='C';
    Mice(7).novelty='S';
    Mice(8).novelty='S';
    Mice(9).novelty='S';
    Mice(10).novelty='S';
    Mice(11).novelty='S';
    Mice(12).novelty='S';
    Mice(1).datanum=26;
    Mice(2).datanum=40;
    Mice(3).datanum=43;
    Mice(4).datanum=21;
    Mice(5).datanum=52;
    Mice(6).datanum=67;
    Mice(7).datanum=50;
    Mice(8).datanum=36;
    Mice(9).datanum=14;
    Mice(10).datanum=57;
    Mice(11).datanum=38;
    Mice(12).datanum=42;
    %AllActLabels=csvread('7day_poke_labels_N1_auto');
    
elseif(strcmp(mouseSet, 'Iku_photometry_DLC'))
    cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry_DLC/
    Mice(1).name='Francisco';
    Mice(2).name='Gardner';
    Mice(3).name='Houston';
    Mice(4).name='Ithaca';
    Mice(5).name='Juneau';
    Mice(6).name='Kennebunk';
    Mice(7).name='Miami';
    Mice(1).novelty='S';
    Mice(2).novelty='S';
    Mice(3).novelty='S';
    Mice(4).novelty='S';
    Mice(5).novelty='S';
    Mice(6).novelty='S';
    Mice(7).novelty='S';
    Mice(1).datanum=50;
    Mice(2).datanum=41;
    Mice(3).datanum=27;
    Mice(4).datanum=37;
    Mice(5).datanum=18;
    Mice(6).datanum=41;
    Mice(7).datanum=36;
    %AllActLabels=csvread('IkuFP_poke_labels_N1_auto',1,0);
elseif(strcmp(mouseSet, 'MSFP_Test'))
    cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/MSFP_Test/
    Mice(1).name='MSFP';
    Mice(1).novelty='S';
    Mice(1).datanum=16;%34; %16
    
    cd ./MSFP/Analyzed_Data
    load('Arena_Obj_Pos.mat') %%%
    cd ../..
end

%AllActLabels(:,4)=(AllActLabels(:,3)-AllActLabels(:,1))./fps;

% Mice(1).name='C1_Akbar';
% Mice(2).name='C2_Emperor';
% Mice(3).name='C3_Piett';
% Mice(4).name='S1_Anakin';
% Mice(5).name='S2_Jabba';
% Mice(6).name='S3_Wedge';
% Mice(1).novelty='C';
% Mice(2).novelty='C';
% Mice(3).novelty='C';
% Mice(4).novelty='S';
% Mice(5).novelty='S';
% Mice(6).novelty='S';
% Mice(1).datanum=32;
% Mice(2).datanum=23;
% Mice(3).datanum=30;
% Mice(4).datanum=24;
% Mice(5).datanum=17;
% Mice(6).datanum=29;
% AllActLabels=csvread('ROTJManualLabels.csv',1,2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation (Richard's code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for miceiter=1%:length(Mice)
    cd(Mice(miceiter).name);
    %cd('Analyzed_Data'); %%%
    %cd('head'); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pathname = cd; %%%
    PathRoot=[pathname '/']; %%%
    filelist=dir([PathRoot,'19*']);%'*rgb.mat']); %%%
    flen = length(filelist); %%%
    %load(filelist(AnalysisDay+0).name); %1 %%%
    cd(filelist(AnalysisDay+0).name)
    currSessionName = dir('session*');
    cd(currSessionName.name)
    %byHand = load('RewardResponse_Nashville_190425R.csv');
    cd ./proc/Analyzed_Data
    
    Labels_MoSeq = load(['MoSeqPos_' Mice(miceiter).name '_' filelist(AnalysisDay+0).name]);
    xPos = Labels_MoSeq(:,1);
    yPos = Labels_MoSeq(:,2);
    Labels17 = Labels_MoSeq(:,3);
    %Labels17 = Labels(:,17);
    load('Arena_Obj_Pos.mat')
    
    % Distance
    Xtime=startframe:endframe;
    Dis=Labels17(startframe:endframe,1); %%%
    SDis=smoothdata(Labels17(startframe:endframe,1),'rloess',Swindow); %%%
%     cd /media/alex/DataDrive1/MoSeqData/MSFP_Test/DataSet180922
%     load('MoSeqDataFrame.mat')
%     uuid_subset = ismember(MoSeqDataFrame.uuid, MoSeqDataFrame.session_uuid(1,:), 'rows');
%     session_uuid_start = [find(uuid_subset,1,'first'); find(~uuid_subset,1,'first')];
%     xPos = MoSeqDataFrame.centroid_x_px(session_uuid_start(2):end);
%     yPos = MoSeqDataFrame.centroid_y_px(session_uuid_start(2):end);
%     
%     Dis_whole=sqrt((obj_center(2,1)-xPos).^2+(obj_center(2,2)-yPos).^2);
%     SDis_whole=smoothdata(Dis_whole,'rloess',Swindow);
%     Dis=Dis_whole(startframe:endframe);
%     SDis=smoothdata(Dis_whole(startframe:endframe),'rloess',Swindow);

    
    SDis_whole=smoothdata(Labels17(:,1),'rloess',Swindow); %%%
    DisMin=islocalmin(SDis_whole);
    DisMax=islocalmax(SDis_whole);
    
    findDisMax = find(DisMax);

    % Derivative of Distance
    DisVelocity=diff(Labels17(startframe:endframe,1)');
%     DisVelocity=diff(Dis);
    DisVelocity=[0 DisVelocity];
    SDisVelocity=smoothdata(DisVelocity,'rloess',Swindow);

    Mice(miceiter).PokingLabels=[];
    Mice(miceiter).BoutStart=[];
    
    windowiter=startframe;
    while windowiter<endframe
        % mouse crosses into radius around object
        if SDis_whole(windowiter)<DisThreshold
            
            minfound=0;
            for wpointer=1:length(Labels17)
                % within bout, find local minimum (poke)
                if minfound==0
                    if DisMin(windowiter+wpointer)==1
                        disp(['minfound: ' num2str(windowiter)])
                        minfound=1;
                        Mice(miceiter).PokingLabels=[Mice(miceiter).PokingLabels windowiter+wpointer];
                        
                        % find beginning of bout (last local maximum before poke)
                        prevBoutStartAll = find(findDisMax<(windowiter+wpointer),1,'last');
                        prevBoutStart = findDisMax(prevBoutStartAll);
                        Mice(miceiter).BoutStart=[Mice(miceiter).BoutStart prevBoutStart];
                    end
                end
                
                % start over when mouse leaves radius
                if SDis_whole(windowiter+wpointer)>DisThreshold
                    break
                end
            end
            windowiter=windowiter+wpointer;
        end
        windowiter=windowiter+1;
    end



    DisPlot=figure(1);
    plot(Xtime,SDis)
    hold on
%     plot(Xtime/fps,Labels(startframe:endframe,18),'LineWidth',1.5)
%     plot(Xtime/fps,Labels(startframe:endframe,19),'LineWidth',1.5)
%     plot(Xtime/fps,Labels(startframe:endframe,20),'LineWidth',1.5)
    scatter(Mice(miceiter).PokingLabels,SDis_whole(Mice(miceiter).PokingLabels))
    scatter(Mice(miceiter).BoutStart,SDis_whole(Mice(miceiter).BoutStart))
    line([startframe endframe], [DisThreshold DisThreshold], 'color', 'k')
%     line([poke_byHand-900 poke_byHand-900]', repmat([0 350]', 1, length(poke_byHand)), 'color', 'k')
    set(DisPlot, 'Position', [44 296 1871 505])
    
    
    correctPokingLabels = Mice(1).PokingLabels(Mice(1).BoutStart(1)<Mice(1).PokingLabels);
    boutPlot = figure(2); 
    hold on
    for int = 1:length(Mice(1).BoutStart)
        plot(xPos(Mice(1).BoutStart(int):correctPokingLabels(int)), ...
             yPos(Mice(1).BoutStart(int):correctPokingLabels(int)))%, 'k')
    end
    plot(obj_center(1), obj_center(2), 'r*')
    
    th = 0:pi/50:2*pi;
    x  = obj_center(1,1);
    y  = obj_center(1,2);
    xunit = radius * cos(th) + x;
    yunit = radius * sin(th) + y;
    plot(xunit, yunit)
    
    axis square
    set(gca, 'YDir', 'reverse')
    
%     DisMin_cut = DisMin(startframe:endframe);
%     DisMax_cut = DisMax(startframe:endframe);
%    DisPlot2=figure(2);
%    plot(Xtime,SDis, 'k-')
%    hold on
%    plot(Xtime(DisMin_cut), SDis(DisMin_cut), 'r.')
%    plot(Xtime(DisMax_cut), SDis(DisMax_cut), 'g.')
%     plot(Xtime/fps,Labels(startframe:endframe,18),'LineWidth',1.5)
%
%     DisPatchY=[0 70 70 0];
%     for actiter=1:length(Mice(miceiter).act(:,1))
%         hold on
%         patch([Mice(miceiter).act(actiter,1) Mice(miceiter).act(actiter,1) Mice(miceiter).act(actiter,3) Mice(miceiter).act(actiter,3)]/fps,DisPatchY,[0.6 0.4 0.9],...
%         'FaceAlpha',0.3, 'EdgeColor','none')
%     end
%     xlim([375 445])
%     ylabel('Distance (cm)')
%     xlabel('time (s)')

    
%     VelPlot=figure(2);
%     plot(Xtime/fps,SDisVelocity*fps);
%     hold on
%     plot(Xtime/fps,zeros(1,length(Xtime)),'LineWidth',1.5);
%     DisPatchY=[-30 50 50 -30];
% 
%     for actiter=1:length(Mice(miceiter).act(:,1))
%         hold on
%         patch([Mice(miceiter).act(actiter,1) Mice(miceiter).act(actiter,1) Mice(miceiter).act(actiter,3) Mice(miceiter).act(actiter,3)]/fps,DisPatchY,[0.6 0.4 0.9],...
%         'FaceAlpha',0.3, 'EdgeColor','none')
%     end
%     xlim([375 445])
%     ylabel('Derivative of Distance cm/s')
%     xlabel('time (s)')

    cd ../..

    %cd .. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
disp('end')

%% Saving

csvwrite(['NoveltyResponse_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '_poke'], ...
    Mice(miceiter).PokingLabels)
csvwrite(['NoveltyResponse_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '_boutStart'], ...
    Mice(miceiter).BoutStart)
saveas(DisPlot, ['FindingPokeDTW_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '.tif'])
saveas(boutPlot, ['FindingPokeDTW_boutTrace_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), '.tif'])
disp('saved')

%save(['NoveltyResponse_Day' num2str(AnalysisDay)])
%poking_labels = extractfield(Mice,'PokingLabels');
%csvwrite([num2str(mouseSet), '_poke_labels_Day', num2str(AnalysisDay), '_auto'], poking_labels)

%% add back in datanum, index, and act fields
load('NoveltyResponse_Day3.mat')

Mice(1).datanum=41;
Mice(2).datanum=29;
Mice(3).datanum=21;
Mice(4).datanum=45;
Mice(5).datanum=19;
Mice(6).datanum=26;
Mice(7).datanum=20;
Mice(8).datanum=30;
AllActLabels=csvread('CvsS_poke_labels_N1_auto')';
    
Mice(1).index=1;
for miceiter=2:length(Mice)
    Mice(miceiter).index=Mice(miceiter-1).index+Mice(miceiter-1).datanum;
    endindex=Mice(miceiter).index+Mice(miceiter).datanum;
end

if endindex ~= length(AllActLabels)+1
   error('index calculation error');
end


for miceiter=1:length(Mice)
   Mice(miceiter).act=AllActLabels(Mice(miceiter).index:Mice(miceiter).index+Mice(miceiter).datanum-1,:);
end

save('NoveltyResponse_Day3.mat')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization (sanity check)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

% load MoSeq position information
cd /media/alex/DataDrive1/MoSeqData/MSFP_Test/DataSet180922
load('MoSeqDataFrame.mat')

% load depth_ts and rgb_ts for conversion
cd session_20180924185040/
delimiter  = ' ';
formatSpec = '%*q%f%[^\n\r]';

filename = 'depth_ts.txt';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
depthts = [dataArray{1}];

filename = 'rgb_ts.txt';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
rgbts = [dataArray{1}];

% load DLC tracking and object position
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/MSFP_Test/MSFP/Analyzed_Data
load('MSFP_Test_180924_rgb.mat')
load('Arena_Obj_Pos.mat')
ppc = 42/6.3;

% extract distance from object for DLC and MoSeq
uuid_subset        = ismember(MoSeqDataFrame.uuid, MoSeqDataFrame.session_uuid(1,:), 'rows');
session_uuid_start = [find(uuid_subset,1,'first'); find(~uuid_subset,1,'first')];

xPos = MoSeqDataFrame.centroid_x_mm(session_uuid_start(2):end);
yPos = MoSeqDataFrame.centroid_y_mm(session_uuid_start(2):end);

xyPos_MoSeq = sqrt((obj_center(2,1)-xPos).^2+(obj_center(2,2)-yPos).^2);
xyPos_DLC   = Labels(:,17);

% zscore to better compare y-values from each array
xPos(isnan(xPos))=0;
yPos(isnan(yPos))=0;
xPos_zscore = zscore(xPos);
yPos_zscore = zscore(yPos);
xLabels_zscore = zscore(Labels(:,14));
yLabels_zscore = zscore(Labels(:,15));

xPos_rescale = rescale(xPos, min(Labels(1:15000,14))-5, max(Labels(1:15000,14))+12);
yPos_rescale = rescale(yPos, min(Labels(1:15000,15))-5, max(Labels(1:15000,15)));

xyPos_MoSeq(isnan(xyPos_MoSeq))=0;
xyPos_MoSeq_zscore = zscore(xyPos_MoSeq);
xyPos_DLC_zscore = zscore(-xyPos_DLC); % -1


% translate automatically identified pokes from 25fps (DLC) to 30fps (MoSeq)
poke_byHand = csvread('NoveltyResponse_byHand_korleki.csv'); % depth_ts units
poke_auto   = load('NoveltyResponse_MSFP_N1')';              % DLC output units
rescaling   = length(xyPos_DLC)/length(depthts);             % (DLC output units / depth actual units)
poke_auto   = round(poke_auto/rescaling);                    % depth_ts units
%rescaling   = length(xyPos_MoSeq)/length(xyPos_DLC);
%x_new = (1:length(xyPos_DLC))*rescaling;

poke_adj = zeros(length(poke_auto),1);
for i = 1:length(poke_auto)
    [d, ix] = min(abs(rgbts(poke_auto(i))-depthts));
    poke_adj(i,1) = ix;
end

% define circle around obj center
th = 0:pi/50:2*pi;
xunit = radius * cos(th) + obj_center(2,1);
yunit = radius * sin(th) + obj_center(2,2);

%%
close all

figure(1)
hold on
%plot(xPos_zscore(1:18000),yPos_zscore(1:18000), 'k-')
%plot(xLabels_zscore(1:15000), yLabels_zscore(1:15000), 'r-')
%plot(xPos_zscore(poke_byHand), yPos_zscore(poke_byHand), 'g*')
%plot(xPos_zscore(poke_adj), yPos_zscore(poke_adj), 'c*')

plot(Labels(1:15000,14), Labels(1:15000,15), 'r-')
plot(xPos_rescale(1:18000), yPos_rescale(1:18000), 'k-')
plot(Labels(poke_auto,14), Labels(poke_auto,15), 'g*')

plot(xunit, yunit, 'r')
legend({'DLC', 'MoSeq', 'poke auto'})
set(gca,'ydir','reverse')


axis square


fig2 = figure(2);
set(fig2, 'Position', [46 26 1865 450])
hold on
%plot(xyPos_MoSeq, 'k-')
%plot(xyPos_MoSeq_zscore, 'k-')
%plot(x_new, xyPos_DLC_zscore, 'r.-')
%plot(poke_byHand, xyPos_MoSeq_zscore(poke_byHand), 'g*')
%plot(poke_adj, xyPos_DLC_zscore(poke_adj), 'c*')
% line([poke_byHand poke_byHand]', repmat([5 45], length(poke_byHand), 1)', ... %[-2.5 2]
%     'color', 'g')

plot(1:length(xyPos_DLC), xyPos_DLC, 'r-')
plot(poke_auto, xyPos_DLC(poke_auto), 'g*')
%line([poke_auto poke_auto]', repmat([5 45], length(poke_auto), 1)', ... %[-2.5 2]
%    'color', 'g')


xlim([0 20000])
legend({'MoSeq', 'DLC', 'poke_byHand', 'poke_auto'})


%% convert DLC output units to rgb_ts units
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry_DLC/

for miceiter = 1:length(Mice)
    cd(['./' Mice(miceiter).name '/Analyzed_Data'])
    filename = dir('NoveltyResponse_*_Day3');
    poke_auto = load(filename.name);
    
    cd ..
    rgbts_list = dir('*rgb_ts.txt');
    rgbts = load(rgbts_list(3).name);
    rgbtsLen = length(rgbts);
    
    csv_list = dir('*csv');
    DLC_output = csvread(csv_list(3).name, 3, 0);
    DLC_out_len = length(DLC_output);
    
    rescaling = rgbtsLen/DLC_out_len;
    poke_adj = round(poke_auto*rescaling);
    
    cd Analyzed_Data
    csvwrite(['NoveltyResponse_', Mice(miceiter).name, '_Day', num2str(AnalysisDay), ...
        '_rgbAdj'], poke_adj)
    
    cd ../..
end
