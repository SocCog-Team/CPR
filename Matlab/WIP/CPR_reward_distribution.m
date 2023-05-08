function f = CPR_reward_distribution(t)

addpath /Users/fschneider/Documents/MATLAB/CircStat2012a

cmap        = [0 0 0; gray(256)];
steps       = .05;
bins        = 0:steps:1;
c           = 0;

% f     	= figure('Units', 'normalized', 'Position', [0 0 .75 .5]);
% u{1}   	= uipanel('position',[0,0,.5,1],'BackgroundColor', [1 1 1]);
% u{2}  	= uipanel('position',[.5,0,.5,1],'BackgroundColor', [1 1 1]);

% Extracte experimental data
for iState = 1:size(t,1)
    if t.trg_shown(iState) == false
        continue
    end
    
    for iTarget = 1:length(t.trg_ts{iState})
        trg_pos     = find(t.frme_ts{iState} >= t.trg_ts{iState}(iTarget), 1, 'first');
        
        if isempty(trg_pos)
            continue
        end
        
        c               = c+1;
        js_dev          = rad2deg(circ_dist(deg2rad(t.js_dir{iState}(trg_pos)),deg2rad(t.rdp_dir(iState))));
        trg_acc(c)      = abs(1 - abs(js_dev / 180));
        trg_conf(c)     = t.js_str{iState}(trg_pos);
        trg_coh(c)      = t.rdp_coh(iState);
        trg_hit(c)   	= t.trg_hit{iState}(iTarget);
    end
end

% Initiate figure
f           	    = figure;
sch                 = scatterhist(trg_acc,trg_conf,...
    'Kernel','on', ...
    'Color','kbr', ...
    'Direction','out', ...
    'Marker','none');

xlim([0 1])
ylim([0 1])

% Calculate reward matrix
acc              	= 0:.001:1;
conf                = 0:.001:1;
rew                 = acc' .* conf;

% Visualise matrix
% f                   = figure;
% h                   = surf(acc,conf,rew);
% h.LineStyle         = 'none';
% ax                  = gca;
% ax.XLabel.String    = 'Accuracy';
% ax.YLabel.String    = 'Confidence';
% ax.ZLabel.String    = '% Reward';
% ax.FontSize         = 16;
% colormap(jet(256))

% Determine arc width for each confidence level
for j = 1:length(conf)
    arc(j)          = 180 - (180 * conf(j));
end

% Cap arc width at target width (2dva == 12.7587deg at chosen position)
idx                 = arc < 12.76;
arc(idx)            = 12.76;

% For each confidence level, calculate minimum accuracy required to hit
% the target at given arc width - normalised values
hit_width_acc       = 1 - ((arc/2) / 180);
hit_width_acc(idx)  = 1 - (12.76/2)/180; % arc width fixed

% Remove position that cannot yield reward from reward matrix
for iAcc = 1:length(acc)
    indx            = conf < hit_width_acc(iAcc);
    rew(iAcc,indx)  = nan;
end

% Plot reward matrix
hold on
ax                  = gca;
im                  = imagesc(acc,conf,rew);
ax.XLabel.String    = 'Accuracy';
ax.YLabel.String    = 'Eccentricity';
ax.FontSize         = 22;
ax.XLim             = [0 1];
ax.YLim             = [0 1];
%     ax.Title.String     = titl{iSub};
cb                  = colorbar;
cb.Label.String     = '% Reward';

% Plot target-wise subject performance
snr                 = unique(trg_coh);
% cl                  = linspace(0,1,size(snr,2));
cl = jet(size(snr,2));

for iCoh = 1:size(snr,2)
    cIdx                = trg_coh == snr(iCoh);
    sc                  = scatter(trg_acc(trg_hit&cIdx),trg_conf(trg_hit&cIdx));
    sc.Marker           = '.';
    sc.SizeData         = 75;
%     sc.MarkerFaceColor  = [1 cl(length(snr)+1-iCoh) 0];
%     sc.MarkerEdgeColor  = [1 cl(length(snr)+1-iCoh) 0];
    sc.MarkerFaceColor  = cl(iCoh,:);
    sc.MarkerEdgeColor  = cl(iCoh,:);
    sc.MarkerEdgeAlpha  = .5;
    sc.MarkerFaceAlpha  = .5;
    hdl(iCoh)           = sc;
    
    sc                  = scatter(trg_acc(~trg_hit&cIdx),trg_conf(~trg_hit&cIdx));
    sc.Marker           = '.';
    sc.SizeData         = 75;
%     sc.MarkerFaceColor  = [1 cl(length(snr)+1-iCoh) 0];
%     sc.MarkerEdgeColor  = [1 cl(length(snr)+1-iCoh) 0];
    sc.MarkerFaceColor  = cl(iCoh,:);
    sc.MarkerEdgeColor  = cl(iCoh,:);
    sc.MarkerEdgeAlpha  = .5;
    sc.MarkerFaceAlpha  = .5;
end

ax.XLim             = [0 1];
ax.XTick            = [0:.2:1];
ax.YTick            = [0:.2:1];
colormap(cmap)

% [lg,obj]         	= legend(hdl, {num2str(round(snr',2))});
% lg.Box            	= 'off';
% lg.Title.String  	= 'SNR';
% lg.Position(1:2)  	= [0 .05];
% 
% for k = length(snr)+1:length(snr)*2
%     obj(k).Children.MarkerSize = 20;
% end

end
