clear
clc
close all

% Config_NovAna_MoSeq_LegoRubber;
% Config_NovAna_combine3;
% Config_NovAna_Ghana;
Config_NovAna_ATLA;

currSet = 'Aviary_DLC'; %%%
currMouse = 'Robin';
% currDate  = 'allFiles'; %190702
timeShift = 901-900; %6000;%901-900; %For Strawberry/Kiwi: 2488/1725;

cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' currSet '/' currMouse])% '/' currDate]) %%
% cd(['/Users/cakiti/Dropbox (Uchida Lab)/Korleki Akiti/ForMelissa/'])
%cd(['/media/alex/DataDrive1/MoSeqData/CvsS_20180831_MoSeq/' currSet '/' currMouse])% '/' currDate]) %%%

%currSessionName = dir('session*'); %%%
%cd([currSessionName.name '/proc']) %%%
path=cd;
PathRoot=[path '/'];

% vn = dir('*.mp4'); %%% '*.avi'
% vn = vn(1).name; %%%
% flen=1; %%%

filelist=dir([PathRoot,'*.csv']); %csv %mp4
flen = length(filelist); %%
fn = filelist(flen).name; %%
vn=[fn(1:end-4-length(networkname_format)) videoname_format(end-3:end)]; %%%%%%%%%%% +1 -4 -9 -12 %%
Labels = csvread(fn,3,0); %%

video=VideoReader(vn);

% for timeStep = 1:100
%     video.CurrentTime = (timeShift-1)/video.FrameRate; %1;
%     frame=readFrame(video);
%     imshow(frame)
%     display(num2str(timeShift+timeStep))
%     pause
%     timeShift=timeShift+1;
% end
% %
video.CurrentTime = (timeShift-1)/video.FrameRate; %1;
frame=readFrame(video);
imshow(frame)
number_of_obj = input('How many objects? Enter integer # 1-5: ');

if(number_of_obj==1)
    arena=zeros(flen,4);
    obj=zeros(flen,4);
    obj_center=zeros(flen,2);
else
    arena=zeros(flen,4);
    obj=cell(flen,number_of_obj);
    obj_center=cell(flen,number_of_obj);
end

for fi =flen:-1:1
    fn = filelist(fi).name; %%
    vn=[fn(1:end-4-length(networkname_format)) videoname_format(end-3:end)]; %%%%%%%%%%%%% +1 -4 -9 -12 %%
    Labels = csvread(fn,3,0); %%
    
    video=VideoReader(vn);
    video.CurrentTime = (timeShift-1)/video.FrameRate; %1;
    frame=readFrame(video);
    
    if fi == flen
        arena_choice=0;%input('Use deflault arena? 1/0: ');
        if arena_choice == 1
            cur_arena = [133.4 31.2 483.4 391.6];
        elseif arena_choice == 0
            fig1 = figure(1);
            set(fig1, 'Position', [545 1 1200 950]);
            cur_arena=Labelrect(frame,'Please Select Arena');
            
            close all
        else
            error('invalid input');
        end
    end
    arena(fi,:)=cur_arena;
    
    
    fig1 = figure(1);
    set(fig1, 'Position', [545 1 1200 950]);
    if(number_of_obj == 1)
        cur_obj  = Labelrect(frame,['Session ' num2str(fi) '; Please Select Object 1']);
        close all
        
        obj(fi,:) = cur_obj;
        cur_obj_center   = 0.5.*[cur_obj(1)+cur_obj(3),cur_obj(2)+cur_obj(4)];
        obj_center(fi,:) = cur_obj_center;
    else
        
        for objNum = 1:number_of_obj
            clf
            fig1 = figure(1);
            set(fig1, 'Position', [545 1 1200 950]);
            cur_obj{objNum} = Labelrect(frame, ...
                ['Please Select Object ' num2str(objNum)]);
            
            obj(fi,objNum)   = cur_obj(objNum);
            cur_obj_center   = 0.5.*[cur_obj{objNum}(1,1) + cur_obj{objNum}(1,3), ...
                cur_obj{objNum}(1,2) + cur_obj{objNum}(1,4)];
            obj_center(fi,objNum) = {cur_obj_center};
        end
        
        close all
        
    end

end

video.CurrentTime = (timeShift-1)/video.FrameRate; %1;
frame=readFrame(video);
LED_yn = input('Select LED position? 0/1: ');
if(LED_yn==1)
    clf
    fig1 = figure(1);
    set(fig1, 'Position', [545 1 1200 950]);
    LEDpos = Labelrect(frame, 'Please Select LED position');
end
close all

%***********************************************************
% Save
%***********************************************************
% pause

%mkdir Analyzed_Data
cd Analyzed_Data_1obj_7cm_nose

clearvars -except arena obj obj_center currSet currMouse LEDpos
save('Arena_Obj_Pos.mat')

cd ..
cd /home/alex/Programs/Novelty_analysis_KA

% %
fourObj = input('Create 4obj file? 0/1: ');
if(fourObj==1)
% modify Arena_Obj_Pos.mat (when only 1 obj in arena, add 3 pseudo-positions at each corner)

cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/' ...
    currSet '/' currMouse '/Analyzed_Data_1obj'])

if isfile('Arena_Obj_Pos_4obj.mat')
    error('this code has already been run')
end

clear
clc
clf

Config_NovAna;

copyfile Arena_Obj_Pos.mat Arena_Obj_Pos_4obj.mat
load('Arena_Obj_Pos.mat')

path=cd;
% PathRoot=path(1:end-13);%[path '/'];
PathRoot=path(1:end-18);
filelist=dir([PathRoot,'*.csv']);
flen = length(filelist);
number_of_obj = 4;

new_obj_center=cell(flen,number_of_obj);
for i = 1:flen
    new_obj_center(i,1) = {obj_center(i,:)};
end

arena_corners = [arena(1,1:2); [arena(1,3) arena(1,2)]; ...
                 [arena(1,3:4)]; [arena(1,1) arena(1,4)]];

for newCoord = 1:flen
    currX = obj_center(newCoord,1);
    currY = obj_center(newCoord,2);
    
    dist_from_corner = [abs(arena_corners(1,1)-currX) abs(arena_corners(1,2)-currY)];
    
    new_obj_center(newCoord,2) = {[arena_corners(2,1)-dist_from_corner(1) currY]};
    new_obj_center(newCoord,3) = {[arena_corners(3,1)-dist_from_corner(1) ...
                                   arena_corners(3,2)-dist_from_corner(2)]};
    new_obj_center(newCoord,4) = {[currX arena_corners(4,2)-dist_from_corner(2)]};
end
             
% check to make sure coordinates are correct
colors = {'r' 'g' 'b' 'c'};
figure(1)
rectangle('Position', [arena(1,1) arena(1,2) arena(1,3)-arena(1,1) arena(1,4)-arena(1,2)])
set(gca, 'ydir', 'reverse')
hold on

for currCorner = 1:4
    plot(arena_corners(currCorner,1), arena_corners(currCorner,2), ...
        [colors{currCorner} '*'])
    
    for currFile = 1:flen
        plot(new_obj_center{currFile,currCorner}(1), new_obj_center{currFile,currCorner}(2), ...
            [colors{currCorner} '*'])
    end
end

obj_center = new_obj_center;

clear new_obj_center
save('Arena_Obj_Pos_4obj.mat')
end

%***********************************************************
% Functions
%***********************************************************

%Select a rectangular ROI
%ROI = [x_upperleft,y_upperleft,x_lowerright,y_lowerright]
function ROI=Labelrect(frame_name,frame_title)
imshow(frame_name,'InitialMagnification',300);
title(frame_title,'FontSize',15);
mouse=imrect;
pos=getPosition(mouse);% x1 y1 w h
ROI=[pos(1) pos(2) pos(1)+pos(3) pos(2)+pos(4)]; 
end
