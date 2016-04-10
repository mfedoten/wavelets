function [tfr,f,t,coi,scales] = wt(sig,varargin)
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
% [tfr,t,f,scales,coi] = wt(sig,opt)
%
% INPUT
% sig    - signal to analyse;
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
% wname    : string
%   Contains name of the wavelet to be used. For FFT-based only analytical
%   Morlet ('morl') is supported. For convolution-based CWT wavelet name can be
%   any supported by MATLAB's 'cwt'-function (see 'waveinfo' for more info).
% fs       : float
%   Sampling frequency. If not specified the default value is chosen to be fs=1,
%   so the resulting time vector is equal to time indices.
% type     : 'fft' (default) | 'conv'
%   Which type of CWT to use: convolution of FFT-based.
% sampling : 'freq' (default) | 'scales'
%   Which sampling to use to construct vector of scales: linear sampling of
%   frequencies or geometrical sampling of scales.
% fmax     : float (default: fmax = fs/2)
%   CWT will be computed up to this frequency (should be less or equal to
%   Nyquist frequency).
% fmin     : float
%   Minimal frequency. If nothing is provided, then it is chosen such as only 
%   50% ofcoefficients are affected by the COI.
% fstep    : float (default: fstep = fmin)
%   If sampling was chosen as 'freq', specifies frequency resolution. If empty
%   the default frequency step is chosen to be equal min. frequency.
% nscales  : int
%   If sampling was chosen as 'scales', specifies desired number of scales.
% chi      : int (default is 85%)
%   If sampling was chosen as 'scales', specifies the persantage of overlap 
%   between wavelet basis, see [1] for details.
% F        : vector of floats
%   An arbitrary vector of frequcncies, is used to create vector of scales by
%   converting it directly.
% norm     : boolean (default: False)
%   If True, then each row of resulting CWT will be normalized by the
%   corresponding scale.
% plot     : boolean (default: False)
%   If True, plots CWT and COI.
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

% read options (if any) into structure and set defaults
if nargin < 2
    opt = struct();
elseif nargin >= 2
    try
        opt = struct(varargin{:});
    catch
        error('Specify properties as one or more name-value pairs.');
    end
end
% set default values for missing options
if ~isfield(opt, 'type'), opt.type = 'fft'; end
if ~isfield(opt, 'wname') && strcmpi(opt.type, 'fft')
    opt.wname = 'morl';
elseif ~isfield(opt, 'wname') && strcmpi(opt.type, 'conv')
    opt.wname = 'cmor1.5-0.95';
end
if ~isfield(opt, 'fs'), opt.fs = 1; end
if ~isfield(opt, 'fmax'), opt.fmax = opt.fs/2; end
if ~isfield(opt, 'chi'),  opt.chi  = 85; end

% sampling of frequencies/scales (freq. default)
if ~isfield(opt, 'sampling') || isfield(opt,'F')
    opt.sampling = 'freq';
elseif strcmpi(opt.sampling,'scales') && isfield(opt,'F')
    warning(['Cannot use both scales and frequency vector. ',...
        'Switching to sampling in frequencies.']);
    opt.sampling = 'freq';
end

N  = length(sig);   % number of points in the signal
dt = 1/opt.fs;      % sampling time
t  = (0:N-1)*dt;    % time vector

%--------------------------- Wavelet parameters ---------------------------
[w0,factor,bound] = wt_params(opt.type,opt.wname,dt);

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
    opt.wname = 'morl';
    cwtStruct = cwtft({sig,dt},'scales',scales,'wavelet',opt.wname,...
        'padmode','zpd');
    coefs = cwtStruct.cfs;
elseif strcmpi(opt.type,'conv')
    coefs  = cwt(sig, scales, opt.wname);
end

% get scalogram (power)
S = abs(coefs.*conj(coefs));

% normalize by the scales
if isfield(opt,'norm') && opt.norm
    S = bsxfun(@rdivide,S,scales(:));
end

% normalize to get PSD
tfr = S./N/w0*2*pi;
% tfr = S./opt.fs;

%---------------------------- Cone of influence ---------------------------
% left and right edges of cone of influence
coi = wt_get_coi(scales, N, bound);

% --------------------------------- Plots ---------------------------------
if isfield(opt,'plot') && opt.plot
    plot_wt(t,f,tfr,coi);
end

end

