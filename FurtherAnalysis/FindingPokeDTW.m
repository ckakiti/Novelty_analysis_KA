%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

AnalysisDay=2; % 3 = first novelty day
fps=15; % 15 rgb; % 25; 30 depth; 

startframe=900;%750;%900;
endframe=18000+startframe;%15000;%18000;
Swindow=40;         % Smooth Window Size
DisThreshold=10;       % Distance threshold

%mouseSet = 'CvsS_180831_DLC';
%mouseSet = '7day_preexposure_combine';
%mouseSet = 'Iku_photometry_DLC';
mouseSet = 'MSFP_Test';

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
% Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for miceiter=1:length(Mice)
    cd(Mice(miceiter).name);
    cd('Analyzed_Data'); 
    %cd('head'); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pathname = cd;
    PathRoot=[pathname '/'];
    filelist=dir([PathRoot,'*rgb.mat']);
    flen = length(filelist);

    load(filelist(AnalysisDay+0).name); %1

    % Distance
    Xtime=startframe:endframe;
    Dis=Labels(startframe:endframe,17);
    SDis=smoothdata(Labels(startframe:endframe,17),'rloess',Swindow);

    SDis_whole=smoothdata(Labels(:,17),'rloess',Swindow);
    DisMin=islocalmin(SDis_whole);

    % Derivative of Distance
    DisVelocity=diff(Labels(startframe:endframe,17)');
    DisVelocity=[0 DisVelocity];
    SDisVelocity=smoothdata(DisVelocity,'rloess',Swindow);

    Mice(miceiter).PokingLabels=[];

    windowiter=startframe;
    while windowiter<endframe
        if SDis_whole(windowiter)<DisThreshold
            minfound=0;
            for wpointer=1:length(Labels)
                if minfound==0
                    if DisMin(windowiter+wpointer)==1
                        minfound=1;
                        Mice(miceiter).PokingLabels=[Mice(miceiter).PokingLabels windowiter+wpointer];
                    end
                end
                if SDis_whole(windowiter+wpointer)>DisThreshold
                    break
                end
            end
            windowiter=windowiter+wpointer;
        end
        windowiter=windowiter+1;
    end



    DisPlot=figure(1);
    plot(Xtime/fps,SDis)
    hold on
    plot(Xtime/fps,Labels(startframe:endframe,18),'LineWidth',1.5)
    hold on
    plot(Xtime/fps,Labels(startframe:endframe,19),'LineWidth',1.5)
    hold on
    plot(Xtime/fps,Labels(startframe:endframe,20),'LineWidth',1.5)
    hold on
    scatter(Mice(miceiter).PokingLabels/fps,SDis_whole(Mice(miceiter).PokingLabels))
    
    DisMin_cut = DisMin(startframe:endframe);
    DisPlot2=figure(2);
    plot(Xtime/fps,SDis, 'k-')
    hold on
    plot(Xtime(DisMin_cut)/fps, SDis(DisMin_cut), 'r.')
    plot(Xtime/fps,Labels(startframe:endframe,18),'LineWidth',1.5)
%
%     pause
%     saveas(DisPlot, ['FindingPokeDTW_fig_Francisco_Hab2.tif'])
    close all

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

    csvwrite(['NoveltyResponse_', Mice(miceiter).name, '_Day', num2str(AnalysisDay)], ...
        Mice(miceiter).PokingLabels)

    cd ..
    cd ..
    %cd .. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
disp('end')

%% Saving
save(['NoveltyResponse_Day' num2str(AnalysisDay)])

poking_labels = extractfield(Mice,'PokingLabels');
csvwrite([num2str(mouseSet), '_poke_labels_Day', num2str(AnalysisDay), '_auto'], poking_labels)


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

%% convert rgb_ts into depth_ts (for fiber photometry analysis)
poke = load('NoveltyResponse_MSFP_N1');

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

pokeAdjust = zeros(1,length(poke));
for i = 1:length(poke)
    [d, ix] = min(abs(rgbts(poke(i))-depthts));
    pokeAdjust(1,i) = ix;
end
