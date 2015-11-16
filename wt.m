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
% fstep    : float (default: fstep = fmin)
%   If sampling was chosen as 'freq', specifies frequency resolution.
% nscales  : int
%   If sampling was chosen as 'scales', specifies desired number of scales. By
%   default nscales is chosen to give 65% overlap between wavelet basis, see [1]
%   for details.
% F        : vector of floats
%   Can be used to create vector of scales, by conversion instead of creating
%   vector inside the function (substitute 'fstep' option).
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

% read options (if any) into structure
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
if ~isfield(opt, 'sampling'), opt.sampling = 'freq'; end
if ~isfield(opt, 'fmax'), opt.fmax = fs/2; end


%--------------------------- Wavelet parameters ---------------------------

% define if we use convolution- or FFT-based (better) CWT?
if strcmpi(opt.type,'fft')
    % FFT-based wavelet transform
    if ~strcmpi(wname,'morl')
        % for now only Morlet wavelet is possible here:
        wname = 'morl';
        warning(['Sorry, this wavelet is not supported right now. ',...
            'Using analytical Morlet']);
    end
    % wavelet center frequency
    w0 = 6;
    % Fourier factor
    factor = 4*pi/(w0+sqrt(2+w0^2));
    % borders of COI
    bound = sqrt(2)/dt;
    
elseif strcmpi(opt.type,'conv')
    % convolution-based wavelet transform
    % center frequency
    w0 = 2*pi*centfrq(wname);
    % scale /freq. transform factor
    factor = dt/centfrq(wname);
    % borders of COI
    bound = wavsupport(wname);
    bound = bound(2);

else
    % unknown type -> error
    error('Unknown type of CWT: %s. Acceptable values are: "fft" and "conv".',...
        opt.sampling);
end


%------------------------ Frequencies/scales vector -----------------------

% How to construct scales vector: by geometrical sampling of scales or
% based on linear sampling of frequencies?
if strcmpi(opt.sampling,'freq') 
    % frequencies
    if isfield(opt, 'F')  % if freq. vector is provided just use it
        f = opt.F;
    else                  % if no, build new freq vector
        % minimal frequency is chosen so that 1/2 of WT coefficients would
        % be affected by COI see Section VI in [1] for details.
        smax = N/(4*bound);
        fmin = 1/(smax*factor);
        
        % chose the frequncy step 
        if isfield(opt,'fstep') && ~isempty(opt.fstep)
            df = opt.fstep;
        else
            df = fmin;
        end
        
        % construct freq. vector
        f  = df:df:opt.fmax;
        % make spacing between scales "smoother"
        if fmin < df
            f = [fmin f];
        elseif fmin > df
            warning(['Frequency step df=%.3f is too small. ',...
               'Consider to switch to default frequency step df=%.3f.'],...
               df,fmin);
        end
    end
    % convert frequency to scales
    scales = 1./(factor*f);
    
elseif strcmpi(opt.sampling,'scales')
    % scales
    % chose the smallest scale based on the highest frequency
    s0 = 1/(factor*opt.fmax);
    % the biggest scale: see Section VI in [1] for details.
    smax = N/(4*bound);
    if isfield(opt, 'nscales') && ~isempty(opt.nscales)
        % if we are given the number of scales
        I = opt.nscales;
        % derive ratio
        K = exp((log(smax) - log(s0)) / (I - 1));
    else
        % overlap (%) between wavelet basis, set it to an arbitrary value
        chi = .65;
        % proportionality constant
        wd = -sqrt(-2*log(chi));
        % constant ratio sc(i+1)/sc(i)
        K = w0/(w0 + wd);
        % derive number of scales
        I = round((log(smax) - log(s0))/log(K) + 1);
    end
    % build scales vector
    scales = zeros(1,I);
    scales(1) = s0;
    for i=2:I
        scales(i) = scales(i-1)*K;
    end
    scales = fliplr(scales);
    % convert scales to frequencies
    f = 1./(factor*scales);
    
else
    % wrong type
    error('Unknown method of constructing scales vector: %s', opt.sampling);
end


%---------------- Calculate wavelet transform coefficients ----------------

if strcmpi(opt.type,'fft')      % FFT-based
    cwtStruct = cwtft({sig,dt},'scales',scales,'wavelet',wname,...
        'padmode','zpd');
    coefs = cwtStruct.cfs;
elseif strcmpi(opt.type,'conv')
    coefs  = cwt(sig, scales, wname);
end
% get scalogram / power spectrum
S = abs(coefs.*conj(coefs));
% scalogram: % of energy at each scale
% tfr = 100*S./sum(S(:));
% tfr    = S./N;
tfr = S./N/w0*2*pi;


%---------------------------- Cone of influence ---------------------------

% left and right edges of cone of influence
border = ceil(bound * scales);
L  = min(floor(N/2), border);
R = max(ceil(N/2), N-border);
coi = [L(:), R(:)];


% --------------------------------- Plots ---------------------------------
if isfield(opt,'plot') && opt.plot
    
    figure;
    fpos = get(gcf,'Position');
    set(gcf,'Position',[fpos(1:2)/4 fpos(3:4)*2]);
    
    % plot wavelets
    contourf(t,f,tfr,20,'edgecolor','none'); set(gca,'YDir','normal');
    if verLessThan('matlab','8.2')
        if exist('brewermap','file') == 2
            colormap(brewermap([],'WhRd'));
        else
            colormap(1-hot);
        end
        cc = 'k';
    else
        cc = 'w';
    end
    colorbar;
    hold on;
    
    % plot COI
    hPatch = patch([L 0 0]*dt,[f f(end) f(1)],min(tfr(:)),'FaceColor','k','EdgeColor',...
        cc,'LineWidth',1.3);
    hatchfill(hPatch, 'cross', 45, 10,cc);
    hPatch = patch([R N N]*dt,[f f(end) f(1)],min(tfr(:)),'FaceColor','k','EdgeColor',...
        cc,'LineWidth',1.3);
    hatchfill(hPatch, 'cross', 45, 10,cc);
    
    % axes labels
    xlabel('Time', 'Fontsize', 14)
    ylabel('Pseudo-frequency', 'Fontsize', 14)
    
%     tightfig; 
end

end
