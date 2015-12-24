function coi = wt_get_coi(scales,N,bound)
% Returns left and right edges of cone of influence (COI): CWT coefficients
% which are affected by edge effect. Returned value is two-column marix, where
% the first column is the left egde and the second one is the right edge of COI.

border = ceil(bound * scales);
L = min(floor(N/2), border);
R = max(ceil(N/2), N-border);
coi = [L(:), R(:)];