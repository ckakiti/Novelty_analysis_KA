%% read fiber photometry .dat file (from Mitsuko)
clear
clc
close all

currMouse = 'Nashville';
currDate  = '190425';
cd(['/home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry2_DLC/'...
    currMouse '/' currDate])

% load LED on data as detected by code
LED_on_name = dir('LED_on');
load(LED_on_name.name)
LED_on_max = LED_on_max';

% load FP data (including TTL pulse ~= LED)
file = strcat('imaging.dat');
file_ID = fopen(file,'r');
CC_cart_analog = fread(file_ID, inf, 'double', 0, 'b');

A = reshape(CC_cart_analog, 2, []);
B = reshape(A,2,3,[]);

cart_GCaMP = B(:,1,:);
cart_GCaMP = reshape(cart_GCaMP,[],1);

cart_tdTom = B(:,2,:);
cart_tdTom = reshape(cart_tdTom,[],1);

TTL = B(:,3,:);
TTL = reshape(TTL,[],1);

TTL_on = crossing(TTL,[],1);
TTL_on = (TTL_on(1:2:end)).';
TTL_ts = TTL_on;

% sanity plotting
close all
hold on
plot(cart_GCaMP, 'g')
plot(cart_tdTom, 'r')
line([TTL_on TTL_on]', repmat([2.5 5]', 1, length(TTL_on)), 'color', 'k')

%% convert pokes (rgb frames) -> FP units (scaled with detected LED_on)
% poke_rgb = dir('NoveltyResponse*poke');
% poke_rgb = load(poke_rgb.name)';
cd ./Analyzed_Data_1obj
rgb_mat = dir('*0000.mat');
rgb_mat = load(rgb_mat.name);
poke_rgb = (1:length(rgb_mat.Labels))';
poke_rgb(poke_rgb<=LED_on_max(1))=[];
poke_rgb(poke_rgb>=LED_on_max(end))=[];

pos_within_TTL = zeros(length(poke_rgb),1);
for curr_poke = 1:length(poke_rgb)
    left_bound = find((poke_rgb(curr_poke)-LED_on_max)>0, 1, 'last');
    LED_dist   = LED_on_max(left_bound+1)-LED_on_max(left_bound);
    TTL_dist   = TTL_on(left_bound+1)-TTL_on(left_bound);
    
    within_frac = (poke_rgb(curr_poke)-LED_on_max(left_bound))/LED_dist;
    within_TTL  = within_frac*TTL_dist;
    pos_within_TTL(curr_poke,1) = round(within_TTL)+TTL_on(left_bound);
end

% if(0)
    % sanity check
    figure(2)
    hold on
    plot(poke_rgb/poke_rgb(end), 'r*')
    plot(pos_within_TTL/pos_within_TTL(end), 'b*')
    
    figure(3)
    plot(cart_GCaMP, 'k')
    line([TTL_on TTL_on]', repmat([3 5], length(TTL_on), 1)')
% end
%% save
cd ..
csvwrite([currMouse, '_' currDate '_rgb_ts'], [poke_rgb pos_within_TTL]) 
% previously '_poke_corrected'

save([currMouse, '_' currDate '_FP'], ...
   'cart_GCaMP', 'cart_tdTom', 'TTL_on', 'pos_within_TTL')

close all