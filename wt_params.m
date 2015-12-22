function [w0,factor,bound] = wt_params(wt_type,wname,dt)

% define if we use convolution- or FFT-based (better) CWT
if strcmpi(wt_type,'fft')                          % FFT-based wavelet transform
    if ~strcmpi(wname,'morl')
        % for now only Morlet wavelet is possible here:
        warning(['Sorry, this wavelet is not supported right now. ',...
            'Using analytical Morlet']);
    end
    % wavelet center frequency
    w0 = 6;
    % Fourier factor
    factor = 4*pi/(w0+sqrt(2+w0^2));
    % borders of COI
    bound = sqrt(2)/dt;
elseif strcmpi(wt_type,'conv')             % convolution-based wavelet transform
    % center frequency
    w0 = 2*pi*centfrq(wname);
    % scale /freq. transform factor
    factor = dt/centfrq(wname);
    % borders of COI
    bound = wavsupport(wname);
    bound = bound(2);
else                                                     % unknown type -> error
    error('Unknown type of CWT: %s. Acceptable values are: "fft" or "conv".',...
        opt.sampling);
end
