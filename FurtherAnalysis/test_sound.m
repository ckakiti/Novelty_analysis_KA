% Test sound output with matlab
clear
clc
close all

cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos/Iku_photometry_DLC/Miami/temp
load('Miami_190125_nidaq.mat')

% create sound vector
x = 1:10000;
a = 1; % amplitude/loudness
b = 0.5; % frequency/pitch
c = 0; % shift
y = a*sin(b*x+c);


% simultaneously plot wave and play sound
fig1 = figure(1);
set(fig1, 'Position', [600 600 1200 450])

plot(y)

axis([x1 x1+dx y1 y2])

%sound(y)