clear all; close all; clc;
addpath('../');

% Example 1: chirp signal
N  = 2048;
fs = 1024;
t  = (0:N-1)/fs;
x  = chirp(t,0,t(end),fs/2,'quadratic')';
% add a bit of noise
x  = x + 0.5*randn(size(x));

% create two subplots
hf = figure;
hf.Position(3:4) = [1.5 1.8].*hf.Position(3:4);
ha1 = subplot('Position',[0.1 0.75 0.8 0.2]);
ha2 = subplot('Position',[0.1 0.1 0.8 0.55]);

% plot signal
axes(ha1);
plot(t,x);
set(gca,'FontSize', 12);
xlabel('Time','FontSize',16);
ylabel('Amplitude','FontSize',16);

% CWT: default parameters
[WT,t_wt,f_wt,scales,coi] = wt(x,1024);
axes(ha2); wt_plot(t,f_wt,WT,coi,hf);

% convolution based wavelets
opts = struct('type','conv');
[WT,t_wt,f_wt,scales,coi] = wt(x,1024,opts);
axes(ha2); wt_plot(t,f_wt,WT,coi,hf);

% FFT-based wavelets, sampling in scales (non-linear frequencies)
opts = struct('type','fft','sampling','scales');
[WT,t_wt,f_wt,scales,coi] = wt(x,1024,opts);
axes(ha2); cla;

wt_plot(t,f_wt,WT,coi,hf);

% [WT,t_wt,f_wt,scales,coi] = wt(x,1024,'name','morl','sampling','freq','norm',1);