function [f,scales] = wt_create_scales(N,bound,factor,w0,opt)

% If min freq. is specified, derive max scale from it, otherwise choose it such
% that 1/2 of WT coefficients would be affected by COI see Section VI in [1] for
% details.
if isfield(opt, 'fmin') && ~isempty(opt.fmin)
    smax = 1/(factor*opt.fmin);
else
    smax = N/(4*bound);
end

% Chose the smallest scale based on the highest frequency
s0 = 1/(factor*opt.fmax);

% construct scale vector
if isfield(opt, 'nscales') && ~isempty(opt.nscales)
    % based on the given number of scales
    I = opt.nscales;
    % derive ratio
    K = exp((log(smax) - log(s0)) / (I - 1));
else
    % based on the given % of overlap between wavelet basis
    chi = opt.chi/100;
    % proportionality constant
    wd = -sqrt(-2*log(chi));
    % constant ratio sc(i+1)/sc(i)
    K = w0/(w0 + wd);
    % derive number of scales
    I = round((log(smax) - log(s0))/log(K) + 1);
end

% build scales vector recursively
scales = zeros(1,I);
scales(1) = s0;
for i=2:I
    scales(i) = scales(i-1)*K;
end

% convert scales to frequencies
scales = fliplr(scales);
f = 1./(factor*scales);

