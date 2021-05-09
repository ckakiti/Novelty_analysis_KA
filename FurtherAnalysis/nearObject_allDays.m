% Extract "mouse is near object" vector from .mat file
%  (generated by Analysis.m)
clear
close all
clc

basefolder = '/media/alex/DataDrive1/MoSeqData/Dataset_20190723/';

Mice_Index_path=[basefolder 'MiceIndex/MiceIndex_Chess.m'];
run(Mice_Index_path);

currSet = 'StandardSetup_combine';

cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/', ...
    currSet, '/'])
folderpath = cd;
folderd = dir(folderpath);
isub = [folderd(:).isdir];
foldernames = {folderd(isub).name}'; 
foldernames(ismember(foldernames,{'.','..','temp'})) = []; 
folderlen=length(foldernames);

for miceiter = 1:length(Mice)
    currMouse = Mice(miceiter).name;
    mouse_idx = find(contains(foldernames,currMouse));
    cd(foldernames{mouse_idx})
    
    % filter out blank sessions (no MSid)
    anyEmpty = find(cellfun(@isempty, {Mice(miceiter).ExpDay.MSid}));
    for emptyiter = 1:length(anyEmpty)
        disp('file missing')
        Mice(miceiter).ExpDay(anyEmpty(emptyiter)) = [];
    end
        
    allRGBts = dir('*rgb_ts');
    for fileiter = 1:length(allRGBts)
        cd 'Analyzed_Data_1obj_8cm_nose'
        currFiles = dir('*rgb_Converted.mat');%'*Converted.mat');
        currLabels = load(currFiles(fileiter).name, 'Labels');
        temp = currLabels.Labels;
        temp2 = temp(:,21);
        temp3 = crossing(temp2-0.5);
        
        cd ..

        delimiter = ' ';
        formatSpec = '%f%[^\n\r]';%'%*q%f%[^\n\r]';
        filename = allRGBts(fileiter).name;
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, ...
            'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
        fclose(fileID);
        rgbts = [dataArray{1}];
        
        Mice(miceiter).ExpDay(fileiter).NearObj_ts = rgbts(temp3);
    end
    
    cd ..
    
end

        % if using MoSeq position to calculate dist from obj
%         centroidInd = find(ismember(MoSeqDataFrame.uuid, ...
%             Mice(miceiter).ExpDay(dayiter).MSid,'rows'));
%         currX = MoSeqDataFrame.centroid_x_mm(centroidInd);
%         currY = MoSeqDataFrame.centroid_y_mm(centroidInd);
%         
%         close all
%         figure(1)
%         hold on
%         plot(currX, currY, 'k')
%         plot(obj_center(1), obj_center(2), 'b*')
%         
%         modify = input('change object center? 0/1: ');
%         if(modify)
%             objCenterX = input('x value?: ');
%             objCenterY = input('y value?: ');
%             plot(objCenterX, objCenterY, 'r*')
%             pause
%             
%             dist = sqrt((objCenterX-currX).^2+(objCenterY-currY).^2);
%         else
%             dist = sqrt((obj_center(1)-currX).^2+(obj_center(2)-currY).^2);
%         end
%         
%         %nearObj = dist<

%%
cd([basefolder '/MiceIndex'])
save('NearObj_ts_***', 'Mice')

%% combine NearObj files from multiple sets of mice
clear
clc

filefolder = '/media/alex/DataDrive1/MoSeqData/Dataset_20190723/MiceIndex';
cd(filefolder)

NearObj_all = [];
nearObj_files = dir('NearObj_ts*');
for neariter = 1:length(nearObj_files)
    load(nearObj_files(neariter).name)
    NearObj_all = [NearObj_all Mice];
end
Mice = NearObj_all;

if(0)
    cd(filefolder)
    save('NearObj_***','Mice')
end