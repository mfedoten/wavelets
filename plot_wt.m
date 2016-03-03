function plot_wt(t,f,tfr,coi,hfig,fontmin)
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
% hfig    : figure handle (optional), if provided plots in this figure, 
%           otherwise creates new figure;
% fontmin : minimal font size used in tick and colorbar labels (optional).


% Set font sizes
% get default font
fs_default = get(0,'DefaultAxesFontSize');
% smallest font for axes
if nargin > 5 && fontmin
    set(0,'DefaultAxesFontSize',fontmin);
else
    set(0,'DefaultAxesFontSize',fs_default+2);
end
% axes ticks
fs_ticks  = get(0,'DefaultAxesFontSize');
% axes labels
fs_labels = fs_ticks + 4;

% Use provided figure, if asked
if nargin > 4 && ishandle(hfig)
    hf = hfig;
else
    hf = figure;
    pos = get(hf,'Position');
    % increase the size a bit
    set(hf,'Position',[pos(1) pos(2) 1.5*pos(3) 1.5*pos(4)]);
end
set(hf,'Render','painters');

% Check if frequency vector is linear or not
if sum(diff(diff(f))) < eps
    % if linear use image
    imagesc(t,f,tfr); axis xy;
    YLabelStr = 'Frequency';
else
    % if non-linear take logarithm
    f = log2(f);    
    % check once again if it's linear or not
    if sum(diff(diff(f))) < 10^-15
        imagesc(t,f,tfr);  axis xy;
    else    % if it's still non-linear, use pcolor
        pcolor(t,f,tfr); shading flat;
    end
    YLabelStr = 'Frequency (log_2)';
end

% Decent tick lengths
set(gca,'ticklength',.5*get(gca,'ticklength'));

% Anotate the plots
ylabel(YLabelStr, 'FontSize', fs_labels);
xlabel('Time', 'FontSize', fs_labels);

% Colorbar
pos = get(gca,'Position');
cb = colorbar; pause(0.5);
if verLessThan('matlab','8.4')
    set(cb, 'TickLength', [0 0], 'FontSize', fs_ticks);
    poscb = get(cb, 'Position');
    set(cb, 'Position', [poscb(1)+poscb(3) poscb(2) poscb(3) poscb(4)]);
else
    set(cb, 'TickLength', 0, 'Box', 'off');
    cb.Ruler.Axle.Visible = 'off';
    cb.Ruler.SecondaryLabel.HorizontalAlignment = 'left';
    cb.FontSize = fs_ticks;
end
set(gca,'Position',pos);

% cone of influence (plot only if exists)
if nargin > 3 && ~isempty(coi)
    L = coi(:,1); R = coi(:,2);
    ccoi = 'w';                 % color to plot COI
    disp_coi(gca,t,L,R,f,ccoi,min(tfr(:)));
end

% set default font size back
set(0,'DefaultAxesFontSize',fs_default);

end

function disp_coi(hAx,tReal,L,R,f,cc,minLvl)
dt = tReal(2) - tReal(1);
axes(hAx);
hold on;
% add ones due to imagesc properties
hPatch = patch([L(1)*dt; L*dt; 0; 0], [-1 f f(end)+1 -1],...
    minLvl, 'FaceColor', 'none');
hatchfill(hPatch, 'cross', 45, 10, cc);
plot([L(1);L]*dt,[-1 f],cc,'linewidth',1.5);
hPatch = patch([R(1); R; length(tReal); length(tReal)]*dt,...
    [-1 f f(end)+1 -1], minLvl, 'FaceColor', 'none');
hatchfill(hPatch, 'cross', 45, 10, cc);
plot([R(1); R]*dt,[-1 f],cc,'linewidth',1.5);

% % for contourf or pcolor plots
% hPatch = patch([L; 0; 0]/fs,[f f(end) 0],minLvl,'FaceColor','none');
% hatchfill(hPatch, 'cross', 45, 10, cc);
% plot(L/fs, f, cc, 'linewidth', 2.5);
% hPatch = patch([R; length(tReal); length(tReal)]/fs,[f f(end) 0],minLvl,...
%     'FaceColor','none');
% hatchfill(hPatch, 'cross', 45, 10, cc);
% plot(R/fs, f, cc, 'linewidth', 2.5);
end