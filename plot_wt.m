function plot_wt(t,f,tfr,coi,hax,flim,fontmin,colmap)
% Plots wavelets. Depending on the scale vector either uses imagesc for linearly
% sampled frequencsies of pcolor for non-linear case. Creates new figure or uses
% provided figure (last input as figure handle). Can also increase font size a
% bit or use provided font for ticks (fontmin input). coi input is used to plot
% COI, it should be two-column matrix, with the first column being the left edge
% of COI and the second one being the right edge.
% INPUT:
% t       : time vector
% f       : frequency vector
% tfr     : CWT coefficients matrix
% coi     : cone of influence (optional); two-column matrix, where the first 
%           column is the right edge of COI and the second one is the left edge;
% flim    : frequency limits, given as two-element vctor: [fmin fmax] OR as a
%           scalar value: fmax.
% hax     : axis handle (optional), if provided plots in this axis, otherwise
%           creates new figure;
% fontmin : minimal font size used in tick and colorbar labels (optional);
% colmap  : colormap (optional); if empty uses default Matlab colormap.


% Use provided axes, if asked
if nargin > 4 && ~isempty(hax) && ishandle(hax)
    hf = get(hax,'Parent');
else
    hf = figure('Units','Centimeters');
    fpos = get(gcf,'Position');
    fpos = [0.6*fpos(1:2) 1.2*fpos(3) 1.2*fpos(3)];
    set(hf,'Position',fpos);
    hax = axes('Units','Centimeters','Position',[0.1*fpos(3:4) 0.8*fpos(3:4)]);
end
axes(hax);

% set Y-limits; plot until max freq., if specified
if nargin > 5 && ~isempty(flim)
    if isscalar(flim)
        yl = [f(1) flim];
    elseif length(flim)==2
        yl = [flim(1) flim(2)];
    end
else
    yl = [f(1) f(end)];
end

% Check if frequency vector is linear or not
if abs(mean(diff(diff(f)))) < eps
    % if linear use image
    imagesc(t,f,tfr); axis xy;
    YLabelStr = 'Frequency (Hz)';
    set(hf,'Render','painters');
else
    % if non-linear take logarithm
    f = log2(f);
    yl = log2(yl);
    % check once again if it's linear or not
    if abs(mean(diff(diff(f)))) < eps
        imagesc(t,f,tfr);  axis xy;
        set(hf,'Render','painters');
    else    % if it's still non-linear, use pcolor
        pcolor(t,f,tfr); shading flat;
    end
    YLabelStr = 'Frequency (log_2)';
end
ylim(yl);

% cone of influence (plot only if exists)
if nargin > 3 && ~isempty(coi)
    L = coi(:,1); R = coi(:,2);
    ccoi = 'w';                 % color to plot COI
    disp_coi(gca,t,L,R,f,ccoi,min(tfr(:)));
end

% change font sizes
if nargin > 6 && ~isempty(fontmin)
    % smallest font is for axes ticks
    set(hax,'FontSize',fontmin);
else
    fontmin = get(hax,'FontSize');
end
fs_labels = fontmin + 2;

% Anotate the plots
ylabel(YLabelStr, 'FontSize', fs_labels);
xlabel('Time (s)', 'FontSize', fs_labels);

% set desired colormap
if nargin > 7 && ~isempty(colmap)
    colormap(hax,colmap);
end

% Colorbar
pos = get(gca,'Position');
cb = colorbar; pause(0.5);
if verLessThan('matlab','8.4')
    set(cb, 'TickLength', [0 0], 'FontSize', fontmin);
    poscb = get(cb, 'Position');
    set(cb, 'Position', [poscb(1)+poscb(3) poscb(2) poscb(3) poscb(4)]);
else
    set(cb, 'TickLength', 0, 'Box', 'off');
    cb.Ruler.Axle.Visible = 'off';
    cb.Ruler.SecondaryLabel.HorizontalAlignment = 'left';
    cb.FontSize = fontmin;
end
set(gca,'Position',pos);

% Decent tick lengths
set(gca,'ticklength',.5*get(gca,'ticklength'));

end

function disp_coi(hAx,tReal,L,R,f,cc,minLvl)
dt = tReal(2) - tReal(1);
axes(hAx);
hold on;
% add ones due to imagesc properties
hPatch = patch([L(1)*dt; L*dt; 0; 0], [f(1)-1 f f(end)+1 f(1)-1],...
    minLvl, 'FaceColor', 'none');
hatchfill(hPatch, 'cross', 45, 10, cc);
plot([L(1);L]*dt,[f(1)-1 f],cc,'linewidth',1.5);
hPatch = patch([R(1); R; length(tReal); length(tReal)]*dt,...
    [f(1)-1 f f(end)+1 f(1)-1], minLvl, 'FaceColor', 'none');
hatchfill(hPatch, 'cross', 45, 10, cc);
plot([R(1); R]*dt,[f(1)-1 f],cc,'linewidth',1.5);
end