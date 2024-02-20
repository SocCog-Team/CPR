function CPR_plotDyadDist(T1, T2, lab)

% Params
fs                  = 16;
lw                  = 2;
nBins               = 30;
col                 = {[0 0 0];[1 0 0]};
alph                = .75;
cnt                 = 0;

% Extract data
snr                 = unique(T1.rdp_coh);

if strcmp(lab, 'eccentricity')
    d1              = cellfun(@nanmean,T1.js_str);
    d2              = cellfun(@nanmean,T2.js_str);
elseif strcmp(lab, 'accuracy')
    for iState = 1:size(T1,1)
        js_dev1        	= rad2deg(circ_dist(deg2rad(T1.js_dir{iState}),deg2rad(T1.rdp_dir(iState))));  % Minimum RDP-Joystick difference
        d1(iState)      = nanmean(abs(1 - abs(js_dev1) / 180)); % Joystick accuracy
        js_dev2        	= rad2deg(circ_dist(deg2rad(T2.js_dir{iState}),deg2rad(T2.rdp_dir(iState))));  
        d2(iState)      = nanmean(abs(1 - abs(js_dev2) / 180));                          
    end
end

figure
hold on

for iCoh = 1:2:length(snr)*2
    clear dat1 dat2
    
    cnt             = cnt+1;
    dat1          	= d1(T1.rdp_coh == snr(cnt));
    dat2         	= d2(T2.rdp_coh == snr(cnt));
    f_auroc(dat1, dat2)
    
    [N,edges]     	= histcounts(dat1,nBins);
%     [N,edges]     	= histcounts(dat1,nBins,'Normalization','pdf');
    center        	= edges(2:end) - (edges(2)-edges(1))/2;
    s1              = stairs(center, (N./max(N))+iCoh, 'Color', [col{1} .75], 'LineWidth', lw);
    
    [N,edges]     	= histcounts(dat2,nBins,'Normalization','pdf');
    [N,edges]     	= histcounts(dat2,nBins,'Normalization','pdf');
    center       	= edges(2:end) - (edges(2)-edges(1))/2;
    s2              = stairs(center, (-N./max(N))+iCoh, 'Color', [col{2} .75], 'LineWidth', lw);
    
    %     p2              = plot(center,(-N./max(N))+iCoh,'Color', [col{2} .75], 'LineWidth', lw);
    %     p2              = plot(center,-N+iCoh,'Color', [col{2} alph], 'LineWidth', lw);
    
    extr(iCoh,:)    = [center(1) center(end)];
end

for iCoh = 1:2:length(snr)*2
    ln          	= line([min(extr(:,1)) max(extr(:,2))],[iCoh iCoh], 'Color', 'k', 'LineStyle', ':');
end

lg                  = legend(unique(T1.ID),unique(T2.ID));
lg.Location         = 'southwest';

ax                  = gca();
ax.XLabel.String    = lab;
ax.YLabel.String    = 'Coherence [%]';
ax.YTick            = 1:2:length(snr)*2;
ax.YTickLabel       = round(snr,2) * 100;
% ax.YLim             = [0 (length(snr)*2)+1];
ax.XLim             = [min(extr(:,1)) max(extr(:,2))];
ax.XLim             = [0 1];
ax.FontSize         = fs;

end