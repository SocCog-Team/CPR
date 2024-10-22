% clear all
close all
cd '/Users/fschneider/ownCloud/CPR_data/Pilot_free_viewing/'

% Load experimental data
load('/Users/fschneider/ownCloud/CPR_data/Pilot_free_viewing/pop_tbl.mat')
subjID      = {'fxs_20210204';'kan_20210204';'cla_210430'};
titl        = {'Human: High confidence';'Human: Low confidence'; 'NHP_Claas'};
cmap        = [0 0 0; gray(256)];
steps       = .05;
bins        = 0:steps:1;
dest_dir    = '/Users/fschneider/ownCloud/Shared/SFB_1528_A1_Project - @SFB1528 [ownCloud, FS, ST]/Summary_plots/';

% f     	= figure('Units', 'normalized', 'Position', [0 0 .75 .5]);
% u{1}   	= uipanel('position',[0,0,.5,1],'BackgroundColor', [1 1 1]);
% u{2}  	= uipanel('position',[.5,0,.5,1],'BackgroundColor', [1 1 1]);

for iSub = 1:3
    
    % Extract subject data
    if iSub < 3
        stbl    = t(strcmp(t.ID,subjID{iSub}),:);
    else
        tmp = load('/Users/fschneider/Documents/MWorks/XBI/cla_20210430_tbl.mat');
        stbl = tmp.t;
    end
    
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
%     ax.Title.String     = titl{iSub};
    cb                  = colorbar;
    cb.Label.String     = '% Reward';
    
    % Plot target-wise subject performance
    snr                 = unique(trg_coh);
    cl                  = linspace(0,1,size(snr,1));
    for iCoh = 1:size(snr,1)
        cIdx                = trg_coh == snr(iCoh);
        sc                  = scatter(acc_exp(trg_hit&cIdx),conf_exp(trg_hit&cIdx));
        sc.Marker           = '.';
        sc.SizeData         = 75;
        sc.MarkerFaceColor  = [1 cl(length(snr)+1-iCoh) 0];
        sc.MarkerEdgeColor  = [1 cl(length(snr)+1-iCoh) 0];
        sc.MarkerEdgeAlpha  = .5;
        sc.MarkerFaceAlpha  = .5;
        hdl(iCoh)           = sc;
        
        sc                  = scatter(acc_exp(~trg_hit&cIdx),conf_exp(~trg_hit&cIdx));
        sc.Marker           = '.';
        sc.SizeData         = 40;
%         sc.MarkerFaceColor  = [cl(iCoh) cl(size(cl,2)+1-iCoh) 0];
%         sc.MarkerEdgeColor  = [cl(iCoh) cl(size(cl,2)+1-iCoh) 0];
        sc.MarkerFaceColor  = [1 cl(length(snr)+1-iCoh) 0];
        sc.MarkerEdgeColor  = [1 cl(length(snr)+1-iCoh) 0];
        sc.MarkerEdgeAlpha  = .5;
        sc.MarkerFaceAlpha  = .5;
    end
    
    ax.XLim             = [0 1];
    ax.XTick            = [0:.2:1];
    ax.YTick            = [0:.2:1];
    colormap(cmap)

    [lg,obj]         	= legend(hdl, {num2str(round(snr,2))});
    lg.Box            	= 'off';
    lg.Title.String  	= 'SNR';
    lg.Position(1:2)  	= [0 .05];
    
%     if iSub < 3
%         indx = 6:10;
%     else
%         indx = 5:8;
%     end
    
    for i = 6:10
        obj(i).Children.MarkerSize = 20;
    end

    print(f_tmp, [dest_dir 'rewMatrix_subj' num2str(iSub)], '-r300', '-dpng');

    %     set(sch,'parent',u{iSub})
    %     u{iSub}.BorderWidth = 0;
    %     close(f_tmp)
    
    %%% Proportional confidence per coherence level
    clear bin_coh frc_coh
    for iBin = 1:size(bins,2)-1
        cIdx                    = [];
        cIdx                    = conf_exp > bins(iBin) & conf_exp <= bins(iBin+1);
        bin_coh{iBin}           = trg_coh(cIdx);
        
        for iCoh = 1:size(snr,1)
            cohIdx          	= [];
            cohIdx            	= bin_coh{iBin} == snr(iCoh);
            frc_conf(iBin,iCoh) = length(bin_coh{iBin}(cohIdx)) / length(bin_coh{iBin});
        end
    end
    
    f                           = figure; hold on
    im                          = imagesc(frc_conf);
    ax                          = gca;
    ax.YLim                     = [1 size(bins,2)-1];
    ax.XLim                     = [.5 length(snr)+.5];
    ax.XTick                    = 1:length(snr);
    ax.XTickLabel               = {round(snr,2)};
    ax.YTick                    = [2 size(bins,2)-2];
    ax.YTickLabel               = {'Low','High'};
    ax.YTickLabelRotation       = 90;
    ax.XLabel.String            = 'Coherence level';
    ax.YLabel.String            = 'Confidence';
    ax.FontSize                 = 22;
    cb                          = colorbar;
    cb.Label.String             = 'Proportion';
    colormap(cmap)
    caxis([0 1])
    axis square
    print(f, [dest_dir 'ConfBins_subj' num2str(iSub)], '-r300', '-dpng');
        
    %%% Proportional accuracy per coherence level
    clear bin_coh frc_acc
    for iBin = 1:size(bins,2)-1
        aIdx                    = [];
        aIdx                    = acc_exp > bins(iBin) & acc_exp <= bins(iBin+1);
        bin_coh{iBin}           = trg_coh(aIdx);
        
        for iCoh = 1:size(snr,1)
            cohIdx          	= [];
            cohIdx            	= bin_coh{iBin} == snr(iCoh);
            frc_acc(iBin,iCoh)  = length(bin_coh{iBin}(cohIdx)) / length(bin_coh{iBin});
        end
    end
    
    
    f                           = figure; hold on
    im                          = imagesc(frc_acc);
    ax                          = gca;
    ax.YLim                     = [1 size(bins,2)-1];
    ax.XLim                     = [.5 length(snr)+.5];
    ax.XTick                    = 1:length(snr);
    ax.XTickLabel               = {round(snr,2)};
    ax.YTick                    = [2 size(bins,2)-2];
    ax.YTickLabel               = {'Low','High'};
    ax.YTickLabelRotation       = 90;
    ax.XLabel.String            = 'Coherence level';
    ax.YLabel.String            = 'Accuracy';
    ax.FontSize                 = 22;
    cb                          = colorbar;
    cb.Label.String             = 'Proportion';
    colormap(cmap)
    caxis([0 1])
    axis square
    print(f, [dest_dir 'AccBins_subj' num2str(iSub)], '-r300', '-dpng');
    
    %%% Average accuracy per confidence level
    for iBin = 1:size(bins,2)-1
        cIdx                    = [];
        cIdx                    = conf_exp > bins(iBin) & conf_exp <= bins(iBin+1);
        bin_acc{iBin}           = acc_exp(cIdx);
        bin_coh{iBin}           = trg_coh(cIdx);
        
        for iCoh = 1:size(snr,1)
            cohIdx          	= [];
            cohIdx            	= bin_coh{iBin} == snr(iCoh);
            avg_acc(iBin,iCoh) = median(bin_acc{iBin}(cohIdx));
        end
    end

    f                           = figure; hold on
    im                          = imagesc(avg_acc);
    ax                          = gca;
    ax.YLim                     = [1 size(bins,2)-1];
    ax.XLim                     = [.5 length(snr)+.5];
    ax.XTick                    = 1:length(snr);
    ax.XTickLabel               = {round(snr,2)};
    ax.YTick                    = [2 size(bins,2)-2];
    ax.YTickLabel               = {'Low','High'};
    ax.YTickLabelRotation       = 90;
    ax.XLabel.String            = 'Coherence level';
    ax.YLabel.String            = 'Confidence';
    ax.FontSize                 = 22;
    cb                          = colorbar;
    cb.Label.String             = 'Avg Accuracy';
    colormap(cmap) 
    caxis([0 1])
    axis square
    print(f, [dest_dir 'AvgAccBins_subj' num2str(iSub)], '-r300', '-dpng');
    
end

%% Helper function

function [acc,conf,rew,trg_hit,trg_coh] = getTargetValue(t)

indx        = logical(t.trg_shown);
js_str      = t.js_str(indx);
js_dir      = t.js_dir(indx);
rdp_dir     = t.rdp_dir(indx);
js_str_ts   = t.js_str_ts(indx);
trg_ts      = t.trg_ts(indx);
trg_hit     = t.trg_hit(indx);
try
    trg_coh = t.trl_coh(indx);
catch
    trg_coh = t.ss_coh(indx);
end

for i = 1:length(js_str)
    trg_pos     = find(js_str_ts{i} >= trg_ts(i), 1, 'first');
    js_dev      = abs(rdp_dir(i) - js_dir{i}(trg_pos-1));
    
    acc(i)      = abs(1 - (js_dev / 180));
    conf(i)     = js_str{i}(trg_pos-1);
end

conf(conf > 1)  = 1;
rew             = conf .* acc;
end