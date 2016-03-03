%% Examples
% An example script showing use cases of the function.
% The examples are shown on two types of signals: chirp signal and sum of three
% sinusoids with different frequencies. This script demonstrates most of the
% options and their use, as well as provided plot function.
clear all; clc;
addpath('../');

%% Example 1: chirp signal
close all;
% create signal
N  = 2048;          % number of points
fs = 1024;          % sampling frequency
t  = (0:N-1)/fs;    % time vector
x  = chirp(t,0,t(end),fs/2,'quadratic')';
% add a bit of noise
x  = x + 0.5*randn(size(x));

% plot signal
hFigSig = figure;
hFigSig.Position(3:4) = [1.5 0.6].*hFigSig.Position(3:4);
plot(t,x);
set(gca,'FontSize', 12);
xlabel('Time','FontSize',16);
ylabel('Amplitude','FontSize',16);

% CWT: default parameters
[WT,t_wt,f_wt,scales,coi] = wt(x);
hFigWav = figure;
hFigWav.Position(3:4) = [1.5 1.5].*hFigWav.Position(3:4);
plot_wt(t_wt,f_wt,WT,coi,hFigWav);
title('Default options','FontSize',16);

% speicify sampling frequency and choose convolution-based method
% you can specify options as name-value pairs
[WT,t_wt,f_wt,scales,coi] = wt(x,'type','conv','fs',fs,'plot',1);
title('Convolution-based CWT','FontSize',16);

% convolution-based CWT, choose another wavelet (for FFT-based CWT only Morlet
% wavelet is available for the moment)
[WT,t_wt,f_wt,scales,coi] = wt(x,'type','conv','wname','mexh','fs',fs,'plot',1);
title('Convolution-based CWT, using mexican hat wavelet','FontSize',16);

% FFT-based wavelets, sampling in scales (non-linear frequencies)
% or you can specify options as structure
opts = struct('type','fft','fs',fs,'sampling','scales','plot',1);
[WT,t_wt,f_wt,scales,coi] = wt(x,opts);
title('FFT-based CWT, non-linear frequencies','FontSize',16);

% the same options, increase number of scales
fprintf('Number of scales before: %d\n', length(scales));
opts = setfield(opts,'nscales',80);
[WT,t_wt,f_wt,scales,coi] = wt(x,opts);
title('Increased number of scales','FontSize',16);
fprintf('Number of scales now: %d\n', length(scales));

%% Example 2: sum of sinusoids
close all;
% create signal
N  = 1024;
fs = 1024;
t  = (0:N-1)/fs;
f1 = 50; f2 = 150; f3 = 250;
x  = sin(2*pi*f1*t) + sin(2*pi*f2*t) + sin(2*pi*f3*t);
% add a bit of noise
x  = x + 0.3*randn(size(x));

% plot signal
hFigSig = figure;
hFigSig.Position(3:4) = [1.5 0.6].*hFigSig.Position(3:4);
plot(t,x);
set(gca,'FontSize', 12);
xlabel('Time','FontSize',16);
ylabel('Amplitude','FontSize',16);

% FFT-based wavelets, sampling in frequency
opts = struct('type','fft','fs',fs,'sampling','freq','plot',1);
[WT,t_wt,f_wt,scales,coi] = wt(x,opts);
title('FFT-based CWT, linear frequencies','FontSize',16);

% specify min frequency
opts = setfield(opts,'fmin',30);
[WT,t_wt,f_wt,scales,coi] = wt(x,opts);
title('Min. frequency = 30 Hz','FontSize',16);

% specify mac frequency
opts = setfield(opts,'fmax',350);
[WT,t_wt,f_wt,scales,coi] = wt(x,opts);
title('Min. frequency = 30 Hz, max. frequency = 350 Hz','FontSize',16);

% change frequency step
fprintf('Frequency step before: %.2f Hz\n', f_wt(2)-f_wt(1));
opts = setfield(opts,'fstep',3);
[WT,t_wt,f_wt,scales,coi] = wt(x,opts);
title('Min. frequency = 30 Hz, max. frequency = 350 Hz, \Deltaf = 3 Hz',...
    'FontSize',16);
fprintf('Frequency step now: %.2f Hz\n', f_wt(2)-f_wt(1));

% the same boundaries but with non-linear frequency
opts.sampling = 'scales';
[WT,t_wt,f_wt,scales,coi] = wt(x,opts);
title('Non-linear freq., min. frequency = 30 Hz, max. frequency = 350 Hz',...
    'FontSize',16);

% normalize by scale
opts = setfield(opts,'norm',1);
[WT,t_wt,f_wt,scales,coi] = wt(x,opts);
title('Normalize by scale', 'FontSize',16);

% decrease overlap between scales
opts = setfield(opts,'chi',65);
[WT,t_wt,f_wt,scales,coi] = wt(x,opts);
title('Overlap between scales is 65%', 'FontSize',16);
