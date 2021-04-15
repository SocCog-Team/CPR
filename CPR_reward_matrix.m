% clear all
close all
cd '/Users/fschneider/ownCloud/CPR_data/Pilot_free_viewing/'

% Load experimental data
load('/Users/fschneider/ownCloud/CPR_data/Pilot_free_viewing/pop_tbl.mat')

subjID 	= {'fxs_20210204';'kan_20210204'};
titl    = {'Human: High confidence';'Human: Low confidence'};
f     	= figure('Units', 'normalized', 'Position', [0 0 .75 .5]);
u{1}   	= uipanel('position',[0,0,.5,1],'BackgroundColor', [1 1 1]);
u{2}  	= uipanel('position',[.5,0,.5,1],'BackgroundColor', [1 1 1]);

for iSub = 1:2
    % Extract subject data
    stbl = t(strcmp(t.ID,subjID{iSub}),:);
    
    % Extracte experimental data
    [acc_exp,conf_exp,rew_exp,trg_hit,trg_coh] = getTargetValue(stbl);
    
    % Initiate figure
    f_tmp           	= figure;
    sch                 = scatterhist(acc_exp,conf_exp,...
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
    ax.YLabel.String    = 'Confidence';
    ax.FontSize         = 16;
    ax.XLim             = [0 1];
    ax.YLim             = [0 1];
    ax.Title.String     = titl{iSub};
    cb                  = colorbar;
    cb.Label.String     = '% Reward';
    
    % Plot target-wise subject performance
    snr                 = unique(trg_coh);
    cl                  = linspace(.5,1,size(snr,1));
    for iCoh = 1:size(snr,1)
        cIdx                = trg_coh == snr(iCoh);
        sc                  = scatter(acc_exp(trg_hit&cIdx),conf_exp(trg_hit&cIdx));
        sc.Marker           = '.';
        sc.SizeData         = 50;
        sc.MarkerFaceColor  = [cl(iCoh) cl(size(cl,2)+1-iCoh) 0];
        sc.MarkerEdgeColor  = [cl(iCoh) cl(size(cl,2)+1-iCoh) 0 ];
        sc.MarkerEdgeAlpha  = .5;
        sc.MarkerFaceAlpha  = .5;
        hdl(iCoh)           = sc;
        
        sc                  = scatter(acc_exp(~trg_hit&cIdx),conf_exp(~trg_hit&cIdx));
        sc.Marker           = '.';
        sc.SizeData         = 2;
        sc.MarkerFaceColor  = [cl(iCoh) cl(size(cl,2)+1-iCoh) 0];
        sc.MarkerEdgeColor  = [cl(iCoh) cl(size(cl,2)+1-iCoh) 0];
        sc.MarkerEdgeAlpha  = .5;
        sc.MarkerFaceAlpha  = .5;
    end
    
    ax.XLim             = [0 1];
    ax.XTick            = [0:.2:1];
    ax.YTick            = [0:.2:1];
    
    set(sch,'parent',u{iSub})
    u{iSub}.BorderWidth = 0;
    close(f_tmp)
end

cmap                    = [0 0 0; gray(256)];
colormap(cmap)


[lg,obj]                = legend(hdl, {num2str(round(snr,2))});
lg.Box                  = 'off';
lg.Title.String         = 'SNR';
lg.Position(1:2)        = [0 .05];

for i = 6:10
    obj(i).Children.MarkerSize = 20;
end

% dest_dir = '/Users/fschneider/ownCloud/Shared/SFB_1528_A1_Project - @SFB1528 [ownCloud, FS, ST]/Summary_plots/';
% print(f, [dest_dir 'rewMatrix_subj'], '-r300', '-dpng');

%% Helper function

function [acc,conf,rew,trg_hit,trg_coh] = getTargetValue(t)

indx        = logical(t.trg_shown);
js_str      = t.js_str(indx);
js_dir      = t.js_dir(indx);
rdp_dir     = t.rdp_dir(indx);
js_str_ts   = t.js_str_ts(indx);
trg_ts      = t.trg_ts(indx);
trg_hit     = t.trg_hit(indx);
trg_coh     = t.trl_coh(indx);

for i = 1:length(js_str)
    trg_pos     = find(js_str_ts{i} >= trg_ts(i), 1, 'first');
    js_dev      = abs(rdp_dir(i) - js_dir{i}(trg_pos-1));
    
    acc(i)      = abs(1 - (js_dev / 180));
    conf(i)     = js_str{i}(trg_pos-1);
end

conf(conf > 1)  = 1;
rew             = conf .* acc;
end