function [tfr,f,t,coi,scales] = wt(sig,fs,wname,varargin)
% Calculates time-frequency representations (TFR) using continuous 1-D wavelet
% transform (CWT). Also computes cone of influence (COI), a region where
% coefficients of CWT are affects be edge effect.
%
% This function allows two ways of computing CWT:
% - based on convolution (standard);
% - based on FFT (more computationally efficient).
%
% The choice of scales can be done also using two approaches:
% - construct linearly spaced *frequency* vector and then convert it to scales;
% - geometrical spacing of *scales* as described in [1].
%
% COI is computed according to the method used for calculation of CWT
% (convolution- or FFT-based):
% - convolution-based: COI is defined as regions which are large than wavelet
% support at each scale.
% - FFT-based: COI is defined as e-folding time for the autocorrelation of
% wavelet power at each scale [2].
%
%
% SYNTAX
% [tfr,f,t,coi,scales] = wt(sig,fs,wname,opt)
%
% INPUT
% sig    - signal to analyse;
% fs     - sampling frequency;
% wname  - name of the wavelet function;
% opt    - options for wavelet transform, can be defined either as structure or
%          as name-value pairs.
% 
% OUTPUT
% tfr    - matrix with scalogram values (energy per coefficient of CWT), where
%          rows correspond to frequency, columns to time.
% f      - frequency vector;
% t      - time vector;
% coi    - cone of influence of specified wavelet function.  Coefficients within
%          the cone should be treated with care. Returns matrix of size Nsc\*2,
%          where Nsc is number of scales, Columns correspond to left and right
%          borders of COI;
% scales - vector of scales.
%
% OPTIONS
% type     : 'fft' (default) | 'conv'
%   Which type of CWT to use: convolution of FFT-based.
% sampling : 'freq' (default) | 'scales'
%   Which sampling to use to construct vector of scales: linear sampling of
%   frequencies or geometrical sampling of scales.
% fmax     : float (default: fmax = fs/2)
%   CWT will be computed up to this frequency (should be less or equal to
%   Nyquist frequency).
% fmin     : float
%   Minimal frequency. If nothing is provided, then it is chosen so only 50% of
%   coefficients are affected by the COI.
% fstep    : float (default: fstep = fmin)
%   If sampling was chosen as 'freq', specifies frequency resolution.
% nscales  : int
%   If sampling was chosen as 'scales', specifies desired number of scales. By
%   default nscales is chosen to give 85% overlap between wavelet basis, see [1]
%   for details.
% chi      : int
%   Persantage of overlap between wavelet basis (default is 85%), see [1].
% F        : vector of floats
%   Can be used to create vector of scales, by conversion instead of creating
%   vector inside the function (substitute 'fstep' option).
% norm     : boolean (default: False)
%   If True, then resulting each row of CWT will be normalized by the
%   corresponding scale.
% plot     : boolean (default: False)
%   If True plots CWT with contour plots, values and scales/frequencies are
%   linearly scaled and COI is displayed with hatched regions.
%
% REFERENCES:
% 1. Jordan, D., Miksad, R. W. & Powers, E. J. Implementation of the 
%    continuous wavelet transform for digital time series analysis. Review 
%    of Scientific Instruments 68, 1484 (1997).
% 2. Torrence, C. & Compo, G. P. A Practical Guide to Wavelet Analysis. 
%    Bulletin of the American Meteorological Society 79, 61?78 (1998).
% 
%
% Copyright Mariia Fedotenkova, 2015, INRIA Nancy.
% Licensed for use under GNU General Public License, Version 2.  See LICENSE for
% details.


%-------------------------- Load and check inputs -------------------------

N  = length(sig);   % number of points in the signal
dt = 1/fs;          % sampling time
t  = (0:N-1)*dt;    % time vector

% read options (if any) into structure and set defaults
if nargin < 4
    opt = struct();
elseif nargin >= 4
    try
        opt = struct(varargin{:});
    catch
        error('Specify properties as one or more name-value pairs.');
    end
end
% set default values for missing options
if ~isfield(opt, 'type'), opt.type = 'fft'; end
if ~isfield(opt, 'fmax'), opt.fmax = fs/2; end
if ~isfield(opt, 'chi'),  opt.chi  = 85; end

% sampling of frequencies/scales (freq. default)
if ~isfield(opt, 'sampling') || isfield(opt,'F')
    opt.sampling = 'freq';
elseif strcmpi(opt.sampling,'scales') && isfield(opt,'F')
    warning(['Cannot use both scales and frequency vector. ',...
        'Switching to sampling in frequencies.']);
    opt.sampling = 'freq';
end

%--------------------------- Wavelet parameters ---------------------------
[w0,factor,bound] = wt_params(opt.type,wname,dt);

%------------------------ Frequencies/scales vector -----------------------
% Three ways to construct scales vector:
% - by converting the provided frequency vector to vector of scales;
% - based on linear sampling of frequencies;
% - by geometrical sampling of scales.

if strcmpi(opt.sampling,'freq')
    % sampling in frequencies
    [f,scales] = wt_create_freq(N,bound,factor,opt);
elseif strcmpi(opt.sampling,'scales')
    % sampling in scales
    [f,scales] = wt_create_scales(N,bound,factor,w0,opt);
else
    % wrong type of sampling
    error('Unknown method of constructing scales vector: %s', opt.sampling);
end

%---------------- Calculate wavelet transform coefficients ----------------
if strcmpi(opt.type,'fft')      % FFT-based
    % for now only morlet wavelet is supported with FFT-based method
    wname = 'morl';
    cwtStruct = cwtft({sig,dt},'scales',scales,'wavelet',wname,...
        'padmode','zpd');
    coefs = cwtStruct.cfs;
elseif strcmpi(opt.type,'conv')
    coefs  = cwt(sig, scales, wname);
end

% get scalogram (power)
S = abs(coefs.*conj(coefs));

% normalize by the scales
if isfield(opt,'norm') && opt.norm
    S = bsxfun(@rdivide,S,scales(:));
end

% normalize to get PSD
tfr = S./N/w0*2*pi;
% tfr = 100*S./sum(S(:));

%---------------------------- Cone of influence ---------------------------
% left and right edges of cone of influence
coi = wt_get_coi(scales, N, bound);

% --------------------------------- Plots ---------------------------------
if isfield(opt,'plot') && opt.plot
    wt_plot(t,f,tfr,coi)
%     figure;
%     fpos = get(gcf,'Position');
%     set(gcf,'Position',[fpos(1:2)/4 fpos(3:4)*2]);
%     
%     % plot wavelets
%     contourf(t,f,tfr,20,'edgecolor','none'); set(gca,'YDir','normal');
%     if verLessThan('matlab','8.2')
%         if exist('brewermap','file') == 2
%             colormap(brewermap([],'WhRd'));
%         else
%             colormap(1-hot);
%         end
%         cc = 'k';
%     else
%         cc = 'w';
%     end
%     colorbar;
%     hold on;
%     
%     % plot COI
%     hPatch = patch([L 0 0]*dt,[f f(end) f(1)],min(tfr(:)),'FaceColor','k','EdgeColor',...
%         cc,'LineWidth',1.3);
%     hatchfill(hPatch, 'cross', 45, 10,cc);
%     hPatch = patch([R N N]*dt,[f f(end) f(1)],min(tfr(:)),'FaceColor','k','EdgeColor',...
%         cc,'LineWidth',1.3);
%     hatchfill(hPatch, 'cross', 45, 10,cc);
%     
%     % axes labels
%     xlabel('Time', 'Fontsize', 14)
%     ylabel('Pseudo-frequency', 'Fontsize', 14)
%     
% %     tightfig; 
end

end

