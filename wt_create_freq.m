function [f,scales] = wt_create_freq(N,bound,factor,opt)

% Either construct a new frequency vector or use the provided one (if any)
if isfield(opt, 'F') && ~isempty(opt.F)
    f = opt.F;
else
    % First find smallest freqyency, which can be resolved -> largest scale and
    % define temporal freq. step, which is equal min. frequency. Minimal frequency
    % is chosen so that 1/2 of WT coefficients would be affected by COI see Section
    % VI in [1] for details.
    smax = N/(4*bound);
    fmin = 1/(smax*factor);
    df   = fmin;
    
    % if the min freq. is specified use it
    if isfield(opt, 'fmin')
        fmin = opt.fmin;
    end
    % use the provided frequncy step
    if isfield(opt,'fstep') && ~isempty(opt.fstep)
        df = opt.fstep;
    end
    
    % construct freq. vector
    f  = fmin:df:opt.fmax;
end

% convert frequency to scales
scales = 1./(factor*f);
